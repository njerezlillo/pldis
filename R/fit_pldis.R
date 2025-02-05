#' Fitting the Discrete Power-Law Model
#'
#' Estimates the parameters of a discrete power-law model, including the lower
#' bound \eqn{x_{\min}} and the scaling parameter \eqn{\alpha}, using both
#' frequentist and Bayesian approaches.
#'
#' @param x A numeric vector of observed data values.
#' @param xm A numeric value specifying the lower bound \eqn{x_{\min}}.
#' The default is NULL. If not specified, the function will automatically
#' estimate \eqn{x_{\min}} from the data.
#' @param bayesian A logical value indicating whether to use a Bayesian approach
#' to estimate \eqn{x_{\min}}. If `TRUE`, the function employs Bayesian estimation
#' to refine the maximum likelihood estimate. If `FALSE`, the frequentist method is used.
#'
#' @details
#' This function estimates the lower bound \eqn{x_{\min}} of a discrete power-law
#' model using the Kolmogorov-Smirnov (KS) statistic, defined by:
#' \deqn{K = \max |F_n(x)-F(x; \hat{\boldsymbol{\theta}})|.}
#'
#' The KS statistic compares the empirical cumulative distribution function
#' of the observed data with the theoretical cumulative distribution function
#' of the discrete power-law model.
#'
#' This approach aims to find the value of \eqn{x_{\min}} that makes the
#' distribution of the data and the power-law model fit as closely as possible.
#' In other words, it chooses \eqn{x_{\min}} so that the cumulative distribution
#' function of the observed data, restricted to values greater than or equal to
#' \eqn{x_{\min}}, is as similar as possible to the CDF of the power-law model
#' fitted to the same range of data (see Clauset et. al. 2007 for details).
#'
#' When `bayesian = TRUE`, a Bayesian framework is applied to estimate both
#' \eqn{x_{\min}} and the scaling parameter \eqn{\alpha}.
#'
#' @return
#' - `fit_xmin_pldis()`: Returns a list containing the estimated lower bound
#'   of the discrete power-law model based on the data.
#' - `fit_pldis()`: Returns a list containing both estimated parameters.
#'
#' @references
#' Jerez-Lillo, N., Rodrigues, F. A., Ferreira, P. H., & Ramos, P. L. (2025).
#' Beyond the Power Law: Estimation, Goodness-of-Fit, and a Semiparametric
#' Extension in Complex Networks. arXiv preprint arXiv:2311.11200. Available at:
#' \url{https://arxiv.org/abs/2311.11200}
#'
#' @examples
#' # Example data
#' x <- c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 10)
#'
#' # Estimation of xmin using the fit_xmin_pldis function:
#' # First, we estimate xmin using the frequentist approach (bayesian = F)
#' fit_xmin_pldis(x, bayesian = F)
#' # Then, we estimate xmin using the Bayesian approach (bayesian = T)
#' fit_xmin_pldis(x, bayesian = T)
#'
#' # Estimation of the full discrete power-law model parameters using fit_pldis:
#' # First, we fix xmin to 1 and use the frequentist approach to estimate
#' # alpha (bayesian = F)
#' fit_pldis(x, xm = 1, bayesian = F)
#' # Next, we fix xmin to 1 and use the Bayesian approach to estimate
#' # both alpha and xmin (bayesian = T)
#' fit_pldis(x, xm = 1, bayesian = T)
#'
#' # If no xmin value is provided, the function will automatically estimate xmin:
#' # Estimation of the full discrete power-law model with the frequentist approach
#' fit_pldis(x, bayesian = F)
#' # Estimation with the Bayesian approach
#' fit_pldis(x, bayesian = T)
#'
#' @seealso [post_pldis] [loglik_pldis]
#'
#'
#' @importFrom graphics hist
#'
#' @export
fit_xmin_pldis <- function (x, bayesian = T) {
  xmins <- sort(unique(x))
  xmins <- xmins[-length(xmins)]

  vec <-  seq(1.1, 6.5, .01)
  zvec <- zeta(vec)
  xmax <- max(x)
  dat <- matrix(0, nrow = length(xmins), ncol = 2)
  z <- x
  vecprior <- 0

  for (xm in 1:length(xmins)) {
    xmin <- xmins[xm]
    z    <- z[z >= xmin]
    n    <- length(z)

    if (xmin == 1) {
      zdiff <- rep(0,length(vec))
    } else {
      zdiff <- apply(rep(t(1:(xmin-1)),length(vec))^-t(kronecker(t(array(1,xmin-1)),vec)),2,sum)
    }

    x_grid <- xmin:xmax
    n_x_grid <- length(x_grid)

    if (bayesian) {
      Jprior <- Vectorize(function(t) Jprior_pldis(t, xmin), "t")
      vecprior <- Jprior(vec)
    }

    L <- vecprior - vec * sum(log(z)) - n * log(zvec - zdiff)
    I <- which.max(L)

    # compute KS statistic
    fit <- cumsum(((x_grid ^ -vec[I])) / (zvec[I] - sum((1:(xmin - 1)) ^ -vec[I])))
    cdi <- cumsum(hist(z, c(min(z) - 1, (xmin + .5):xmax, max(z) + 1), plot = F)$counts / n)
    dat[xm,] <- c(max(abs(fit - cdi)), vec[I])
  }

  I     <- which.min(dat[, 1])
  xmin  <- xmins[I]
  alpha <- dat[I, 2]

  return(xmin)
}

#' @rdname fit_xmin_pldis
#' @export
fit_pldis <- function (x, xm = NULL, bayesian = T) {
  if (is.null(xm)) xm <- fit_xmin_pldis(x, bayesian)
  alphas <- seq(1.1, 6.5, .01)

  if (bayesian) {
    l <- Vectorize(function(z) post_pldis(z, xm, x[x >= xm]), "z")
  } else {
    l <- Vectorize(function(z) loglik_pldis(z, xm, x[x >= xm]), "z")
  }

  return(list(xmin = xm, alpha = alphas[which.max(l(alphas))]))
}
