#' Observed-Information Standard Error of the Scaling Exponent
#'
#' Computes the standard error of an estimator of the scaling exponent from
#' the observed information.
#'
#' @param f A one-argument function of \eqn{\alpha} returning the
#' log-likelihood (or Fisher-penalized log-likelihood) to be differentiated,
#' e.g. \code{function(z) loglik_pldis(z, xmin, x)} or
#' \code{function(z) loglik_fd_pldis(z, xmin, x)}.
#' @param ahat A numeric value at which to evaluate the observed information,
#' typically the corresponding estimate of \eqn{\alpha} (e.g. the MLE or
#' \eqn{\hat\alpha_{Fd}}).
#'
#' @details
#' The observed information is minus the second derivative of \code{f} at
#' \code{ahat}, computed numerically via \code{numDeriv::hessian}. The
#' resulting standard error is:
#' \deqn{\text{SE}(\hat\alpha) = \left[-f''(\hat\alpha)\right]^{-1/2}.}
#'
#' @return A numeric value representing the observed-information standard
#' error of the scaling exponent.
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
#' ahat <- fit_pldis(x, xm = 1, penalized = FALSE)$alpha
#' se_observed(function(z) loglik_pldis(z, 1, x), ahat)
#'
#' ahat_fd <- fit_pldis(x, xm = 1, penalized = TRUE)$alpha
#' se_observed(function(z) loglik_fd_pldis(z, 1, x), ahat_fd)
#'
#' @seealso [se_expected] [loglik_pldis] [loglik_fd_pldis]
#'
#' @importFrom numDeriv hessian
#'
#' @export
se_observed <- function(f, ahat) 1 / sqrt(-numDeriv::hessian(f, ahat)[1, 1])
