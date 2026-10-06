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
  - `pdf = scale *fraction * np.exp(-0.5 * ((x - mu) / sigma) ** 2)` ... Probability for core gaus, $\text{Scale} * \text{Fraction} * e^{-0.5((x-\mu)/\sigma)^2}$
    - `scale` = 10000.
    - `fraction` = 0.02 ... $A(E_{\mu})$ from the paper (how much each Gaussian contributes to the total PDF area)
    - `mu` = 0.
    - `sigma` = 1000. ... $\sigma_1(E_{\mu})$ from the paper
  - `pdf_tail = scale*(1-fraction) * np.exp(-0.5 * ((x - mu_tail) / sigma_tail) ** 2)` ... Probability for tail gaus
    - `mu_tail` = 0.
    - `sigma_tail` = 14100. ... $\sigma_2(E_{\mu})$
    - '1-fraction` = (1-0.02) ... $1-A(E_{\mu})$ from the paper (how much each Gaussian contributes to the total PDF area)

Double Gaussian fit is performed with a function `dgaus_fitResult = double_gaus_fit_5param(muon_gen2, 'dtheta_yz_prime', xrange=[-500,500], energyBounds=[0.7, 1.5], name=r"$\theta_{yz}^'$",nbins=100,core=True,tail=True)`
- `energyBounds=[0.7, 1.5]` ... $E_{\mu}$


### Root's Optimizer for Fitting
Function `th1_from_series` is defined to optimize the tuning parameters (I think).


### Comment/Questions
- In the code, several values are hardcoded, so I am not sure where they come from. But I also think they are initial guess for fitting. 
- I assume that extra $\frac{1}{\sqrt{2\pi}}$ in the paper is absorbed into `scale` parameter in the code. I am not sure the reason/purpose of it.

