#' Discrete Fisher-Penalized Log-Likelihood of the Power-Law Model
#'
#' Computes the discrete Fisher-penalized log-likelihood for the power-law
#' model, obtained by adding the Firth-type bias-reduction penalty derived
#' from the Fisher information to the ordinary log-likelihood.
#'
#' @param alpha A numeric value representing the scaling parameter \eqn{\alpha},
#' which must be greater than 1.
#' @param xmin A numeric value specifying the lower bound \eqn{x_{\min}}.
#' @param x A numerical vector of observed data.
#'
#' @details
#' This objective combines the log-likelihood of the discrete power-law
#' model with the Fisher-type penalty for \eqn{\alpha}. Given a dataset
#' \eqn{x} and a scaling parameter \eqn{\alpha}, it is defined as:
#' \deqn{\mathcal{L}_{Fd}(\alpha) = \mathcal{L}(\alpha) + p_1(\alpha),}
#' where \eqn{\mathcal{L}(\alpha)} is the ordinary log-likelihood and
#' \eqn{p_1(\alpha)} is the discrete Fisher-type penalty.
#'
#' This implementation relies on the Hurwitz zeta function, provided by the
#' VGAM package, to accurately compute the required terms.
#'
#' @return A numeric value representing the discrete Fisher-penalized
#' log-likelihood evaluated at the given parameters and observed data.
#'
#' @references
#' Firth, D. (1993). Bias reduction of maximum likelihood estimates.
#' *Biometrika*, 80(1), 27-38.
#'
#' Jerez-Lillo, N., Rodrigues, F. A., Ferreira, P. H., & Ramos, P. L. (2025).
#' Bias Reduction and Goodness-of-Fit in Discrete Power-Law Models.
#' Submitted for publication.
#'
#' @examples
#' # Example data
#' x <- c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 10)
#'
#' # Evaluating the discrete Fisher-penalized log-likelihood for different
#' # parameter values
#' loglik_fd_pldis(2.0, 1, x)
#' loglik_fd_pldis(2.5, 1, x)
#' loglik_fd_pldis(2.0, 2, x)
#' loglik_fd_pldis(2.5, 2, x)
#'
#' @seealso [fisher_penalty_pldis] [loglik_pldis] [loglik_fc_pldis]
#'
#' @export
loglik_fd_pldis <- function (alpha, xmin, x) {
  l <- fisher_penalty_pldis(alpha, xmin) + loglik_pldis(alpha, xmin, x)
  return(l)
}
