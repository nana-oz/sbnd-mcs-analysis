import os
import multiprocessing
from multiprocessing import Pool

import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse

# local modules
from modules.TrajectoryMCSFitter import *

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s:%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('logs/make_mcsdfs.log')
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


def get_mcs_fit(idx):
    # Get TP and SP MCS fits
    fwd_tp, bwd_tp = fitMcs(tp_df.loc[idx], PID, momdep, angres, cutseg)
    fwd_sp, bwd_sp = fitMcs(sp_df.loc[idx], PID, momdep, angres, cutseg)
    
    # print("tp, fwd", fwd_tp['p'], "bwd", bwd_tp['p'], 
    #     "truth", tp_df.loc[idx].truth_p.unique()[0], 
    #     "range", tp_df.loc[idx].range_p.unique()[0])
    
    # Store results
    ret= {
        'fidx': idx[0],
        'tidx': idx[1], 
        'truth_p': tp_df.loc[idx].truth_p.unique()[0],
        'tp_fwd_p': fwd_tp['p'], 'tp_bwd_p': bwd_tp['p'],
        'tp_fwd_p_unc': fwd_tp['pUnc'], 'tp_bwd_p_unc': bwd_tp['pUnc'],
        'tp_nsegs': len(tp_df.loc[idx]),
        'tp_range_p': tp_df.loc[idx].range_p.unique()[0],
        'sp_fwd_p': fwd_sp['p'], 'sp_bwd_p': bwd_sp['p'],
        'sp_fwd_p_unc': fwd_sp['pUnc'], 'sp_bwd_p_unc': bwd_sp['pUnc'],
        'sp_nsegs': len(sp_df.loc[idx]),
        'sp_range_p': sp_df.loc[idx].range_p.unique()[0],
        'dirx': tp_df.loc[idx].dirx.unique()[0],
        'diry': tp_df.loc[idx].diry.unique()[0],
        'dirz': tp_df.loc[idx].dirz.unique()[0],
        'flipped': tp_df.loc[idx].flipped.unique()[0]
    }

    return ret
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run MCS analysis on calib ntuple files')
    parser.add_argument('-i', '--input', help='list of input files', default=-1)
    parser.add_argument('-o', '--output', help='output file', default=-1)
    parser.add_argument('-m', '--momdep', help='momentum dependence', default=True)
    parser.add_argument('-a', '--angres', help='angular resolution', default=5)
    parser.add_argument('-c', '--cutseg', help='number of segments to cut at track end', default=0)

    args = parser.parse_args() 

    input_df = args.input
    output_file = args.output
    cutseg = args.cutseg
    momdep = args.momdep
    angres = args.angres

    tp_df = pd.read_hdf(input_df, key="tp_df")
    sp_df = pd.read_hdf(input_df, key="sp_df")

    tp_idxs = tp_df.groupby(level=[0,1]).max().index
    sp_idxs = sp_df.groupby(level=[0,1]).max().index
    assert np.all(tp_idxs == sp_idxs)

    results = []
    with Pool(processes=min(20,CPU_COUNT)) as pool:
        for ret in tqdm(pool.imap(get_mcs_fit, tp_idxs), total=len(tp_idxs)):
            results.append(ret)

    # Create dataframe from results
    mcs_df = pd.DataFrame(results).set_index(['fidx', 'tidx'])

    mcs_df.to_hdf(output_file, key="mcs_df")