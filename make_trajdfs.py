import os
import multiprocessing
from multiprocessing import Pool
import uproot
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import argparse
import random

# local modules
from modules.TrajectoryMCSFitter import *
from modules.load_pandas import *
from modules.RangeE import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s:%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('logs/make_trajdfs.log')
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.propagate = False

# turn off FutureWarning and RuntimeWarning
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

CPU_COUNT = multiprocessing.cpu_count()


if CPU_COUNT == 0:
    CPU_COUNT = os.cpu_count()

def branches(plane):

    hitplane = "hits{}".format(plane)
    tpbranches = [hitplane+".rr", hitplane+".pitch", hitplane+".tp.x", hitplane+".tp.y", hitplane+".tp.z",
                hitplane+".dir.x", hitplane+".dir.y", hitplane+".dir.z"]  
    spbranches = [hitplane+".h.sp.x", hitplane+".h.sp.y", hitplane+".h.sp.z",
                    hitplane+".h.time", hitplane+".h.wire", hitplane+".h.tpc", hitplane+".h.hasSP"]
    htruthbranches = ["truth.p.true"+hitplane+".wire", "truth.p.true"+hitplane+".time", "truth.p.true"+hitplane+".pitch",
                        "truth.p.true"+hitplane+".rr", "truth.p.true"+hitplane+".plane", "truth.p.true"+hitplane+".tpc",
                        "truth.p.true"+hitplane+".p.x", "truth.p.true"+hitplane+".p.y", "truth.p.true"+hitplane+".p.z",]

    trkbranches = ["selected", "length", "dir.x", "dir.y", "dir.z",
                    "start.x", "start.y", "start.z", "end.x", "end.y", "end.z"]
    trktruthbranches = ["truth.p.start.x", "truth.p.start.y", "truth.p.start.z", "truth.p.end.x", "truth.p.end.y", "truth.p.end.z",
                    "truth.p.startE", "truth.p.genE", "truth.p.contained"]

    hitbranches = tpbranches + spbranches 
    trkbranches = trkbranches + trktruthbranches

    return hitbranches, trkbranches


def get_stopping_trajs(fidx):
    thisfile = files[fidx]
    events = uproot.open(thisfile+":caloskim")["TrackCaloSkim"]
    hitbranches, trkbranches = branches(plane)
    pointsdf = loadbranches(events["hits{}".format(plane)], hitbranches)["hits{}".format(plane)]
    trackdf = loadbranches(events, trkbranches)

    # is reco direction flipped?
    truth_len = np.linalg.norm(trackdf.truth.p.start - trackdf.truth.p.end, axis=1)
    truth_dir = (trackdf.truth.p.end - trackdf.truth.p.start) / truth_len[:, np.newaxis]
    reco_dir = trackdf.dir
    dir_dot = reco_dir.x * truth_dir.x + reco_dir.y * truth_dir.y + reco_dir.z * truth_dir.z
    trackdf["flipped"] = dir_dot 
    # trackdf.loc[(dir_dot < 0), "flipped"] = True
    # trackdf.loc[(dir_dot >= 0), "flipped"] = False
    # trackdf.loc[(dir_dot.isna()), "flipped"] = np.nan

    # collect columns needed for Trajectory into a dataframe with a single level of columns
    tp_df = pointsdf.tp.copy()
    tp_df = tp_df.droplevel(1, axis=1)
    tp_df["rr"] = pointsdf["rr"]
    tp_df["pitch"] = pointsdf["pitch"]
    tp_df["tpc"] = pointsdf[("h","tpc","")]
    tp_df["wire"] = pointsdf[("h","wire","")]
    tp_df["time"] = pointsdf[("h","time","")]

    sp_df = pointsdf.h.sp.copy()
    sp_df["rr"] = pointsdf["rr"]
    sp_df["pitch"] = pointsdf["pitch"]
    sp_df["tpc"] = pointsdf[("h","tpc","")]
    sp_df["wire"] = pointsdf[("h","wire","")]
    sp_df["time"] = pointsdf[("h","time","")]

    this_tp_dfs = []
    this_sp_dfs = []
    for tidx in trackdf[trackdf.selected == 0].index: # stopping tracks
        # skip if track is too short
        if trackdf.loc[tidx].length.values[0] < 15:
            continue

        # make trajectories
        tp_traj = Trajectory(tp_df.loc[tidx])
        sp_traj = Trajectory(sp_df.loc[tidx])

        # get segments needed for MCS fit
        prep_tp = fitMcsPrepare(tp_traj)
        prep_sp = fitMcsPrepare(sp_traj)

        # add useful track-level info
        # truth p
        try:
            truth_p = np.sqrt(trackdf.loc[tidx].truth.p.startE.values[0]**2 - mass(PID)**2)
            prep_tp["truth_p"] = truth_p
            prep_sp["truth_p"] = truth_p
        except TypeError:
            logger.info("no truth p info available for track {} in file {}, skipping...".format(tidx, thisfile))
            continue

        # p from reco length
        try:
            prep_tp["trk_length"] = tp_traj.length
            range_p = muonRangeP(tp_traj.length)/1e3 # [GeV]
            prep_tp["range_p"] = range_p

            prep_sp["trk_length"] = sp_traj.length
            range_p = muonRangeP(sp_traj.length)/1e3 # [GeV]
            prep_sp["range_p"] = range_p
        except ValueError:
            # TODO: interpolation fails in some cases, not sure if this is scipy-specific
            logger.info("range estimation failed for track {} in file {}, skipping...".format(tidx, thisfile))
            continue

        # add flipped info
        prep_tp["flipped"] = trackdf.loc[tidx].flipped.values[0]
        prep_sp["flipped"] = trackdf.loc[tidx].flipped.values[0]

        # track direction
        prep_tp["dirx"] = trackdf.loc[tidx].dir.x.values[0]
        prep_tp["diry"] = trackdf.loc[tidx].dir.y.values[0]
        prep_tp["dirz"] = trackdf.loc[tidx].dir.z.values[0]

        prep_sp["dirx"] = trackdf.loc[tidx].dir.x.values[0]
        prep_sp["diry"] = trackdf.loc[tidx].dir.y.values[0]
        prep_sp["dirz"] = trackdf.loc[tidx].dir.z.values[0]

        # add file, track indices, and breakpoint indices
        prep_tp["fidx"] = fidx
        prep_tp["tidx"] = tidx
        prep_tp.set_index(["fidx", "tidx", "breakpoint"], inplace=True)

        prep_sp["fidx"] = fidx
        prep_sp["tidx"] = tidx
        prep_sp.set_index(["fidx", "tidx", "breakpoint"], inplace=True)

        this_tp_dfs.append(prep_tp)
        this_sp_dfs.append(prep_sp)

    return this_tp_dfs, this_sp_dfs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run MCS analysis on calib ntuple files')
    parser.add_argument('-i', '--input', help='list of input files', default=-1)
    parser.add_argument('-o', '--output', help='output file', default=-1)
    parser.add_argument('-p', '--plane', help='which plane', default=2)
    args = parser.parse_args() 

    input_files = args.input
    output_file = args.output
    plane = args.plane

    #  "hist_gen_g4_detsim_reco1_reco2-4f5ce275-882b-7f85-5f6b-c800b59d4ce4.root"
    #  "hist_gen_g4_detsim_reco1_reco2-0b687a4e-6dcd-5d03-bcc4-e700896ef4d7.root"
    files = []
    with open("hist_gen_g4_detsim_reco1_reco2-0b687a4e-6dcd-5d03-bcc4-e700896ef4d7.root", "rb") as f:
        for line in f:
            files.append(line.strip())
    # a few prints to help me debug whether the files are being read in correctly
    print("Number of files: ", len(files))
    print("Range of files: ", range(len(files)))
    if len(files) != 0:
        print("files list is FILLED")
    else:
        print("files list is EMPTY")


    tp_dfs = []
    sp_dfs = []

    nfiles = len(files)
    with Pool(processes=min(20,CPU_COUNT)) as pool:
        for ret in tqdm(pool.imap(get_stopping_trajs, range(nfiles)), total=nfiles):
            this_tp_dfs, this_sp_dfs = ret
            tp_dfs.extend(this_tp_dfs)
            sp_dfs.extend(this_sp_dfs)

    tp_df = pd.concat(tp_dfs)
    sp_df = pd.concat(sp_dfs)

    tp_df.to_hdf(output_file, key="tp_df", mode="w")
    sp_df.to_hdf(output_file, key="sp_df", mode="a")


