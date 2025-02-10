
# pldis Package

<!-- badges: start -->
[![R-CMD-check](https://github.com/njerezlillo/pldis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/njerezlillo/pldis/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
![Lifecycle: experimental](https://img.shields.io/badge/Lifecycle-Experimental-orange)
<!-- badges: end -->

This package offers a collection of tools for fitting the discrete power-law model, utilizing both frequentist and Bayesian approaches for parameter estimation and model evaluation. The package includes:

- Log-likelihood: Calculates the log-likelihood function for the discrete power-law model.
- Jeffreys prior: Defines a non-informative prior for Bayesian estimation.
- Posterior distribution: Computes the posterior distribution of the model parameters.
- Fitting: Provides methods for estimating model parameters, with the option to apply either Bayesian or frequentist approaches.

## Progress status

- [x] Set up the package structure  
- [x] Write functions  
- [x] Document functions
- [x] Write examples for each function in the package
- [x] Check the documentation
- [x] Publish on GitHub  
- [ ] Complete the "Example" section on GitHub
- [ ] Distribute on CRAN

## Installation

You can install the package using :

``` r
# install.packages("devtools")
devtools::install_github("njerezlillo/pldis")
```

## Example

This section presents a concise example illustrating the use of the main functions of the package.

We begin by defining a dataset of 13 observations:

``` r
library(pwpldis)

x <- c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 10)
```

Next, we estimate the discrete power-law model for the dataset `x` using the `fit_pldis` function with both approaches to compare the results. The argument `xm = 1` specifies that the lower bound is fixed at 1. Consequently, the function returns the estimated scaling parameter under both estimation methods:

``` r
fit_pldis(x, xm = 1, bayesian = FALSE)
fit_pldis(x, xm = 1, bayesian = TRUE)
```

In many contexts, however, the power-law model does not apply to the entire dataset but rather holds from a certain lower bound onward. Therefore, we aim to determine the threshold that ensures the best fit between the observed data and the power-law model. Specifically, this threshold is chosen so that the cumulative distribution function of the observed data is as similar as possible to the cumulative distribution function of the fitted power-law model. This alignment is assessed using the Kolmogorov-Smirnov statistic  (see Clauset *et. al.* 2007 for details).

<p align="center">
  <img src="KS.png" alt="">
</p>

Our function allows this threshold to be estimated automatically by omitting the `xm` argument. In this case, the function selects the lower bound that best fits the model while estimating the scaling parameter using either Bayesian or frequentist methods:

``` r
fit_pldis(x, bayesian = FALSE)
fit_pldis(x, bayesian = TRUE)
```

## Citation

To cite `pldis` package in publications, please use the following format:

Jerez-Lillo N (2025). *pldis: Fitting Discrete Power-Law Model*. R package version 1.0.0, [https://github.com/njerezlillo/pldis](https://github.com/njerezlillo/pldis).

For LaTeX users, the corresponding BibTeX entry is:

```bibtex
@Manual{
  title = {pldis: Fitting Discrete Power-Law Model},
  author = {Nixon Jerez-Lillo},
  year = {2025},
  note = {R package version 1.0.0},
  url = {https://github.com/njerezlillo/pldis},
}
```

## References  

[**Beyond the Power Law: Estimation, Goodness-of-Fit, and a Semiparametric Extension in Complex Networks**](https://arxiv.org/abs/2311.11200)  
*Jerez-Lillo N., Rodrigues F. A., Ferreira P. H., Ramos P. L.*  
arXiv preprint arXiv:2311.11200 

[**Power-law distributions in empirical data**](https://doi.org/10.1137/070710111)  
*Clauset A., Shalizi C. R., Newman M. E. J.*  
SIAM Review, 51(4), 661–703 (2009)
