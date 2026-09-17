#' Fitting the Discrete Power-Law Model
#'
#' Estimates the parameters of a discrete power-law model, including the lower
#' bound \eqn{x_{\min}} and the scaling parameter \eqn{\alpha}, using either
#' the ordinary maximum likelihood estimator (MLE) or the discrete
#' Fisher-penalized MLE (Firth-type bias reduction).
#'
#' @param x A numerical vector of observed data.
#' @param xm A numeric value specifying the lower bound \eqn{x_{\min}}.
#' The default is `NULL`. If not specified, the function will automatically
#' estimate \eqn{x_{\min}} from the data.
#' @param penalized A logical value indicating whether to estimate \eqn{\alpha}
#' (and, when `xm` is not supplied, to select \eqn{x_{\min}}) using the
#' discrete Fisher-penalized MLE, \eqn{\hat\alpha_{Fd}} (Firth-type bias
#' reduction; see [loglik_fd_pldis]). The default is `TRUE`. If `FALSE`, the
#' ordinary (unpenalized) MLE is used instead.
#'
#' @details
#' This function estimates the parameters of a discrete power-law model in
#' two steps:
#'
#' First, the lower bound \eqn{x_{\min}} is estimated using the Kolmogorov-Smirnov
#' (KS) statistic, defined as:
#' \deqn{K = \max |F_n(x) - F(x; \hat{\boldsymbol{\theta}})|.}
#'
#' This method identifies the value of \eqn{x_{\min}} that best aligns the
#' empirical data distribution with the power-law model. Specifically, it selects
#' \eqn{x_{\min}} such that the cumulative distribution function (CDF) of the
#' observed data, restricted to values \eqn{x \geq x_{\min}}, closely matches the
#' CDF of the fitted power-law model over the same range (see Clauset et al., 2007
#' for details). When `penalized = TRUE`, the scaling exponent evaluated at
#' each candidate threshold is the discrete Fisher-penalized MLE rather than
#' the ordinary MLE.
#'
#' Second, the scaling parameter \eqn{\alpha} is estimated at the selected
#' (or supplied) \eqn{x_{\min}}, again using either the ordinary MLE or the
#' discrete Fisher-penalized MLE according to `penalized`.
#'
#' @return
#' - `fit_xmin_pldis()`: Returns the estimated lower bound based on the data.
#' - `fit_pldis()`: Returns a list containing both estimated parameters.
#'
#' @references
#' Jerez-Lillo, N., Rodrigues, F. A., Ferreira, P. H., & Ramos, P. L. (2025).
#' Bias Reduction and Goodness-of-Fit in Discrete Power-Law Models.
#' Submitted for publication.
#'
#' @examples
#' # Example data
#' x <- c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 10)
#'
#' # Estimation of xmin:
#' # First, we estimate xmin using the ordinary MLE (penalized = FALSE)
#' fit_xmin_pldis(x, penalized = FALSE)
#' # Then, we estimate xmin using the discrete Fisher-penalized MLE
#' # (penalized = TRUE)
#' fit_xmin_pldis(x, penalized = TRUE)
#'
#' # Estimation of the full discrete power-law model parameters:
#' # First, we fix xmin to 1 and use the ordinary MLE to estimate
#' # alpha (penalized = FALSE)
#' fit_pldis(x, xm = 1, penalized = FALSE)
#' # Next, we fix xmin to 1 and use the discrete Fisher-penalized MLE to
#' # estimate alpha (penalized = TRUE)
#' fit_pldis(x, xm = 1, penalized = TRUE)
#'
#' # If no xmin value is provided, the function automatically estimates it:
#' # Estimation of the full discrete power-law model with the ordinary MLE
#' fit_pldis(x, penalized = FALSE)
#' # Estimation with the discrete Fisher-penalized MLE
#' fit_pldis(x, penalized = TRUE)
#'
#' @seealso [loglik_fd_pldis] [loglik_pldis]
#'
#' @importFrom graphics hist
#'
#' @export
fit_xmin_pldis <- function (x, penalized = TRUE) {
  xmins <- sort(unique(x))
  xmins <- xmins[-length(xmins)]

  vec <-  seq(1.1, 6.5, .01)
  zvec <- zeta(vec)
  xmax <- max(x)
  dat <- matrix(0, nrow = length(xmins), ncol = 2)
  z <- x
  vecpenalty <- 0

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

    if (penalized) {
      penfun <- Vectorize(function(t) fisher_penalty_pldis(t, xmin), "t")
      vecpenalty <- penfun(vec)
    }

    L <- vecpenalty - vec * sum(log(z)) - n * log(zvec - zdiff)
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
fit_pldis <- function (x, xm = NULL, penalized = TRUE) {
  if (is.null(xm)) xm <- fit_xmin_pldis(x, penalized = penalized)
  alphas <- seq(1.1, 6.5, .01)

  if (penalized) {
    l <- Vectorize(function(z) loglik_fd_pldis(z, xm, x[x >= xm]), "z")
  } else {
    l <- Vectorize(function(z) loglik_pldis(z, xm, x[x >= xm]), "z")
  }

  return(list(xmin = xm, alpha = alphas[which.max(l(alphas))]))
}
