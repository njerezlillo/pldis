#' Log-Likelihood of the Discrete Power-Law Model
#'
#' Computes the log-likelihood of a discrete power-law model given a lower bound,
#' scaling parameter, and observed dataset.
#'
#' @param alpha A numeric value representing the scaling parameter \eqn{\alpha},
#' which must be greater than 1.
#' @param xmin A numeric value specifying the lower bound \eqn{x_{\min}}.
#' @param x A numeric vector of observed data values.
#'
#' @details The computation relies on the Hurwitz zeta function from the
#' VGAM package.
#'
#' @return A numeric value representing the log-likelihood of the data
#' given the specified parameters.
#'
#' @examples
#' # Example data
#' x <- c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 10)
#'
#' # Evaluating the log-likelihood for different parameter values
#' loglik_pldis(2.0, 1, x)
#' loglik_pldis(2.5, 1, x)
#' loglik_pldis(2.5, 2, x)
#' loglik_pldis(2.5, 2, x)
#'
#' @importFrom VGAM zeta
#'
#' @export
loglik_pldis <- function (alpha, xmin, x) {
  n <- length(x)
  l <- -n * log(zeta(alpha, shift = xmin)) - alpha * sum(log(x))
  return(l)
}
