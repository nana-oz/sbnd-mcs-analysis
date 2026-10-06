# Double Gaussian function from MicroBooNE paper
$f(\theta; E_{\mu}) =  \frac{A(E_{\mu})}{\sigma_1(E_{\mu})\sqrt{2\pi}}e^{\theta^2/2{\sigma}^2_1(E_{\mu})} + \frac{1 − A(E_{\mu})}{\sigma_2 (E_{\mu})\sqrt{2\pi}}e^{\theta^2/2σ^2_2(E_{\mu})}$
- $\theta$ is the measured scattering angle
- $E_{\mu}$ is the energy of the muon
- $\sigma_1(E_{\mu}) = \sigma_{\text{pred}}(E_\mu) = \sqrt{\frac{2}{3}\kappa^2(E_\mu)\sigma_H^2(E_\mu) + \sigma_{\text{res}}^2}$
- $\sigma_H \approx \frac{S_2}{p c \beta} = \frac{S_2}{E_\mu \left( 1 - \frac{m^2}{E_\mu^2} \right)}$, where $S_2 = 13.6$ MeV and $m \approx 105.7$ MeV (muon mass)

# What the original `mcs_MCSimStudy.ipynb` contains

## Prototype
### Single Gaussian
Single Gaussian is plotted using `plt.plot(x, pdf, color='red', linewidth=2, label='Fitted Gaussian')`, where:

- $x$, `x = np.linspace(-500, 500, nbins)`. In the code, the variable $x$ corresponds to the scattering angle $\theta$ from the paper.
  
- Probability density function, `pdf = stats.norm.pdf(x, mu, sigma) * norm` <br>
  $\text{PDF}(x) = \frac{1}{\sigma \sqrt{2 \pi}} e^{-\frac{1}{2}(\frac{x-\mu}{\sigma})^2}$ where initial guesses are:
  - `area` = 0.78
  - `mu` = 0.452 ... mean (or center peak) of the scattering angle distribution.
  -  `sigma` = 95
 
  and the scaling factor `norm` is number of data points ($N$) * the bin width ($\Delta x$)


### Double Gaussian
Double Gaussian is plotted using ` plt.plot(x, double_pdf, color='orange', linewidth=2, label=f'Fitted Double Gaussian (χ²/ndf = {round(chi2_ndf,2)})', zorder=4)`, where:

- Probability density function, `double_pdf = pdf + pdf_tail` where,
  - `pdf = scale *fraction * np.exp(-0.5 * ((x - mu) / sigma) ** 2)` ... Probability for core gaus, $\text{Scale} * \text{Fraction} * e^{-0.5((x-\mu)/\sigma)^2}$ with initial guesses:
    - `scale` = 10000.
    - `fraction` = 0.02 ... $A(E_{\mu})$ from the paper (how much each Gaussian contributes to the total PDF area)
    - `mu` = 0.
    - `sigma` = 1000. ... $\sigma_1(E_{\mu})$ from the paper
  - `pdf_tail = scale*(1-fraction) * np.exp(-0.5 * ((x - mu_tail) / sigma_tail) ** 2)` ... Probability for tail gaus with initial guesses:
    - `mu_tail` = 0.
    - `sigma_tail` = 14100. ... $\sigma_2(E_{\mu})$
    - `1-fraction` = (1-0.02) ... $1-A(E_{\mu})$ from the paper (how much each Gaussian contributes to the total PDF area)

Double Gaussian fit is performed with a function `dgaus_fitResult = double_gaus_fit_5param(muon_gen2, 'dtheta_yz_prime', xrange=[-500,500], energyBounds=[0.7, 1.5], name=r"$\theta_{yz}^'$",nbins=100,core=True,tail=True)`
- `energyBounds=[0.7, 1.5]` ... $E_{\mu}$


### Root's Optimizer for Fitting
#### Overview
1. Hardcode initial values (initial guesses)
2. Convert Data and Functions to C++ ROOT
3. Optimize parameters using `.Fit()`

#### Detail
`h_dtheta_xz.Fit(func1, "S")`. `.Fit()` is used for parameter fitting. 
- Ref: https://root.cern/manual/fitting/
- It adjusts 6 parameters (`scale`, `fraction`, `mu`, `sigma`, `mu_tail`, `sigma_tail`) at once to lower $\chi^2$.
- Each variable corresponds to: (p0 = `scale`, p1 = `fraction`, p2 = `mu`, p3 = `sigma`, p4 = `mu_tail`, p5 = `sigma_tail`). 

Function `th1_from_series` is defined to optimize the tuning process. Details are in Appendix.


## Prototype --> Stage 0
Function `fit_sigma_res(series, sig_H, initSubList, fit_xmin=-500, fit_xmax=500, nbins=1000, rebin=None, sigLims=None)` performs a ROOT double-Gaussian fit on angular residual data and subtracts Multiple Coulomb Scattering (σH​) to extract the detector's intrinsic baseline resolution (σres​).
- This function is used in later stages

- `sig_res1 = ((sigma_1_0**2) - ((2/3)*(sig_H**2)))**0.5` ... $\sigma_{\text{res1}} = \sqrt{\sigma_1^2 - \frac{2}{3}\sigma^2_H}$
- `sig_res2 = ((sigma_2_0**2) - ((2/3)*(sig_H**2)))**0.5`
and
- `sig_res1_err = ((sigma_1_0*sigma_1_0_err)/(((sigma_1_0**2)-(2/3)*(sig_H**2))**0.5))**0.5`
- `sig_res2_err = ((sigma_2_0*sigma_2_0_err)/(((sigma_2_0**2)-(2/3)*(sig_H**2))**0.5))**0.5`

then
- `sigma_1_0` and `sigma_1_0_err` are determined by (the same process as the Prototype):
  - setting initial guess
  - then fitting it with `.Fit()` function 

In the paper:
- $\sigma_{\text{res}} = \sqrt{\sigma_{\text{pred}}^2 - \frac{2}{3}\kappa^2(E_\mu)\sigma^2_H}$


## Stage 0 (Detector Calibration)
### Stage 0 -- Version 0 (not complete)
- Goal: Extract $\sigma_{\text{res}}$ in $\sigma_{\text{pred}}^2(E_\mu) = \frac{2}{3}\kappa^2(E_\mu)\sigma_H^2(E_\mu) + \sigma_{\text{res}}^2$.
- Energy selection: High energy ($0.7\text{--}1.5\text{ GeV}$) where scattering is minimal.
- Outcome: Intrinsic Detector Resolution ($\sigma_{\text{res}}$)

- Detector Angular Resolution ($\sigma_{\text{res}}$): Hardware noise, wire spacing, and reconstruction smearing. It is constant and independent of track energy.
- Multiple Coulomb Scattering ($\sigma_H$): Pure physics of muons interacting with Liquid Argon nuclei. It scales inversely with energy ($1/E_\mu$).


Function `double_gaus_subplots` is set in the stage 0. Inside the function:
- `sig_H = 13.6/(E_mean*(1-((mu_mass_MeV**2)/(E_mean**2))))` = $\sigma_H = \frac{13.6}{E_{\mu}\left(1-\frac{m^2}{E^2_{\mu}}\right)}$
  - `13.6` ... $S_2 = 13.6 MeV$
  - `mu_mass` ... muon mass in MeV. (`mu_mass_MeV = 105.7`)
  - `E_mean` ... is calculated with `E = (((p**2.)+(mu_mass**2.))**0.5)` where `mu_mass = .1057 #GeV/c^2`, `p` is muon momenta.


### Stage 0 -- Version 3
Function `double_gaus_stageZero_v3(df, column, energyBounds=[0.1, 0.9], maxBin=[0.7, 1.5], sigLims=None, eStep=0.1, xranges=None, yranges=None, initParams=[], nbins=1000, fitRange=None, errdf=False, core=True, tail=False, rebin=None)`. 
- Calculation of `sig_H` is the same as v0. `sig_H = 13.6 / (E_mean * (1 - ((mu_mass_MeV**2) / (E_mean**2))))`.

Step:
- Separate cases into $yz'$ projection ($\theta'_{yz}$) and $xz'$ projection.
  - $yz'$ is more sensitive to $\vert{}v_x\vert{} \to 0$. --> needs |v_x| binning
- 


## Stage 1 (Tail Tuning & Physics)
- Goal: Validate scattering predictions ($\sigma_H$) across energy.
- Energy selection: Low energy ($0.1\text{--}0.9\text{ GeV}$) where scattering dominates.
- Outcome: Low-Energy Model Validation & Tuning ($\chi^2/\text{ndf}$)

Function `sigma1_fit(series,sig_res_arr,sig_H,meanEnergy,nbins=1000,fitRanges=None,initSubList=None,rebin=None)`.


# Appendix 

### Explanation of how the Root's Optimizer for Fitting work (memo):

Function Definition & Parameters
```
def th1_from_series(series, name=r"$d\theta_{xz}$", title="dQ/dx;dQ/dx [A.U./cm];Counts",
                    nbins=100, xmin=-250.0, xmax=250.0, weights=None, edges=None):
```
- `series`: The input Pandas Series containing raw numerical values (e.g., scattering angles).
- `name`: ROOT's internal object identifier.
- `title`: ROOT title string formatted as "Main Title; X-axis Label; Y-axis Label".
- `nbins`, `xmin`, `xmax`: Standard fixed-width binning options.
- `weights`: Optional array of event weights.
- `edges`: Optional explicit bin boundary array for non-uniform bin widths.

Data Cleaning (NaN Handling)
- `vals = series.to_numpy(dtype=float)` Converts the Pandas Series into a clean 1D NumPy array of standard C-compatible floating-point numbers.
- `vals = vals[~np.isnan(vals)]` Removes all missing or invalid (NaN) entries from the array using a boolean filter mask.

```
if weights is not None:
  w = np.asarray(weights, dtype=float)
  w = w[~np.isnan(series.to_numpy())]
else:
  w = None
```
If weights are provided:
- Converts weights to a NumPy array of floats.
- Applies the exact same NaN mask derived from the original series to ensure the weight array length stays aligned with vals.
Sets w to None if no weights were passed, signaling that the data is unweighted.

Bin Edge Calculations
```
if edges is None:
  edges = np.linspace(xmin, xmax, nbins + 1)
else:
  edges = np.asarray(edges, dtype=float)
  nbins = len(edges) - 1
```
- If no custom bin edges were specified, generates nbins + 1 evenly spaced boundary values between xmin and xmax. (For 100 bins, you need 101 edge markers).
- If custom bin edges were provided:
  - Converts edges into a NumPy array.
  - Updates nbins to equal the number of bins implied by those boundaries (len(edges) - 1).

Fast Binning with NumPy
- `counts, edges = np.histogram(vals, bins=edges, weights=w)` Uses NumPy's vectorized C implementation to sort the raw data into the specified bins and sum up the counts (and optional weights). This step executes orders of magnitude faster in Python than filling a ROOT histogram point-by-point in a loop (h.Fill(x)).

ROOT Histogram Construction
- `h = ROOT.TH1F(name, title, nbins, array('d', edges))` Instantiates a 1D single-precision ROOT histogram (TH1F). Passing array('d', edges) (a Python array of double-precision floats) allows ROOT to support variable-width binning structures seamlessly.
- `h.Sumw2()` Tells ROOT to store and calculate the sum of squared weights ($\sum w^2$) for every bin. This ensures ROOT tracks uncertainties correctly if the histogram is scaled, normalized, or weighted later.

Populating Bins & Setting Uncertainties
```    
    for i, c in enumerate(counts, start=1):
        h.SetBinContent(i, float(c))
        # Set sqrt(N) if unweighted; else keep Sumw2 weights
        if w is None:
            h.SetBinError(i, np.sqrt(c))
```
- Loops through each bin count c from the NumPy result:
  - start=1: Important because ROOT histograms use 1-based bin indexing (Bin 1 to Bin N). Bin 0 is reserved for Underflow and Bin N+1 for Overflow.
- Sets the height of ROOT bin i equal to the corresponding count $c$ computed by NumPy.
- If the dataset is unweighted, manually assigns Poisson statistical uncertainty ($\sigma_i = \sqrt{N_i}$) to bin i. ROOT's MINUIT optimizer uses these error values to weight each bin when evaluating $\chi^2$.

Return Output
- `return h` Returns the fully configured and populated ROOT.TH1F histogram, ready to be passed directly to h.Fit(...).
  


