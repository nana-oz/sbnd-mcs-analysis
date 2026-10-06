## Double Gaussian function from MicroBooNE paper
$f(\theta; E_{\mu}) =  \frac{A(E_{\mu})}{\sigma_1(E_{\mu})\sqrt{2\pi}}e^{\theta^2/2{\sigma}^2_1(E_{\mu})} + \frac{1 − A(E_{\mu})}{\sigma_2 (E_{\mu})\sqrt{2\pi}}e^{\theta^2/2σ^2_2(E_{\mu})}$
- $\theta$ is the measured scattering angle
- $E_{\mu}$ is the energy of the muon

## What the original `mcs_MCSimStudy.ipynb` contains
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
    - '1-fraction` = (1-0.02) ... $1-A(E_{\mu})$ from the paper (how much each Gaussian contributes to the total PDF area)

Double Gaussian fit is performed with a function `dgaus_fitResult = double_gaus_fit_5param(muon_gen2, 'dtheta_yz_prime', xrange=[-500,500], energyBounds=[0.7, 1.5], name=r"$\theta_{yz}^'$",nbins=100,core=True,tail=True)`
- `energyBounds=[0.7, 1.5]` ... $E_{\mu}$


## Root's Optimizer for Fitting
Function `th1_from_series` is defined to optimize the tuning parameters. Details are in Appendix.

## Appendix 

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
  


