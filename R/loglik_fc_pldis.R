#' Continuous-Penalty Log-Likelihood of the Discrete Power-Law Model
#'
#' Computes an alternative Fisher-penalized log-likelihood for the discrete
#' power-law model, using the simpler penalty derived from the Fisher
#' information of the continuous power-law model as an approximation to
#' the discrete one.
#'
#' @param alpha A numeric value representing the scaling parameter \eqn{\alpha},
#' which must be greater than 1.
#' @param xmin A numeric value specifying the lower bound \eqn{x_{\min}}.
#' @param x A numerical vector of observed data.
#'
#' @details
#' This objective combines the log-likelihood of the discrete power-law
#' model with a continuous-penalty approximation for \eqn{\alpha}. Given a
#' dataset \eqn{x} and a scaling parameter \eqn{\alpha}, it is defined as:
#' \deqn{\mathcal{L}_{Fc}(\alpha) = \mathcal{L}(\alpha) + p_2(\alpha),}
#' where \eqn{\mathcal{L}(\alpha)} is the ordinary log-likelihood and
#' \eqn{p_2(\alpha)} is given by:
#' \deqn{p_2(\alpha) = -\log(\alpha - 1)}
#' This continuous approximation is computationally simpler than the exact
#' discrete penalty (see [fisher_penalty_pldis]).
#'
#' This implementation relies on the Hurwitz zeta function, provided by the
#' VGAM package, to accurately compute the required terms.
#'
#' @return A numeric value representing the continuous-penalty log-likelihood
#' evaluated at the given parameters and observed data.
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
#' # Evaluating the continuous-penalty log-likelihood for different
#' # parameter values
#' loglik_fc_pldis(2.0, 1, x)
#' loglik_fc_pldis(2.5, 1, x)
#' loglik_fc_pldis(2.0, 2, x)
#' loglik_fc_pldis(2.5, 2, x)
#'
#' @seealso [loglik_pldis] [loglik_fd_pldis]
#'
#' @export
loglik_fc_pldis <- function (alpha, xmin, x) {
  l <- -log(alpha - 1) + loglik_pldis(alpha, xmin, x)
  return(l)
}
