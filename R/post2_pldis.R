#' (2nd) Posterior Distribution of the Discrete Power-Law Model
#'
#' Computes an alternative posterior distribution for a discrete power-law model
#' using the Jeffreys prior derived from the continuous version of the model.
#'
#' @param alpha A numeric value representing the scaling parameter \eqn{\alpha},
#' which must be greater than 1.
#' @param xmin A numeric value specifying the lower bound \eqn{x_{\min}}.
#' @param x A numerical vector of observed data.
#'
#' @details
#' This posterior distribution is obtained by combining the likelihood function
#' of the discrete power-law model with Jeffreys prior for \eqn{\alpha} derived
#' from the continuous version of the model. Given a dataset \eqn{x} and a
#' scaling parameter \eqn{\alpha}, the posterior is proportional to:
#' \deqn{\pi_2(\alpha \mid x) \propto L(\alpha \mid x) \pi_2(\alpha),}
#' where \eqn{L(\alpha \mid x)} is the likelihood function and \eqn{\pi_2(\alpha)}
#' is given by:
#' \deqn{\pi_2(\alpha) \propto \dfrac{1}{\alpha - 1}}
#'
#' This implementation relies on the Hurwitz zeta function, provided by the
#' VGAM package, to accurately compute the required terms.
#'
#' @return A numeric value representing the alternative posterior distribution
#' evaluated at the given parameters and observed data.
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
#' # Evaluating the 2nd posterior distribution for different parameter values
#' post2_pldis(2.0, 1, x)
#' post2_pldis(2.5, 1, x)
#' post2_pldis(2.0, 2, x)
#' post2_pldis(2.5, 2, x)
#'
#' @seealso [loglik_pldis]
#'
#' @export
post2_pldis <- function (alpha, xmin, x) {
  post2 <- -log(alpha - 1) + loglik_pldis(alpha, xmin, x)
  return(post2)
}
