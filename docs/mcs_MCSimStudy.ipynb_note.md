### Single Gaussian

Single Gaussian is plotted using `plt.plot(x, pdf, color='red', linewidth=2, label='Fitted Gaussian')`, where:

- $x$
  `x = np.linspace(-500, 500, nbins)`
  
- Probability density function
  `pdf = stats.norm.pdf(x, mu, sigma) * norm #Probability Density Function`
  $\text{PDF}(x) = \frac{1}{\sigma \sqrt{2 \pi}} e^{-\frac{1}{2}(\frac{x-\mu}{\sigma})^2}$



