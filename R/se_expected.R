#' Expected-Information Standard Error of the Scaling Exponent
#'
#' Computes the asymptotic standard error of an estimator of the scaling
#' exponent, obtained from the Fisher information of the discrete
#' power-law model.
#'
#' @param alpha A numeric value representing the scaling parameter \eqn{\alpha},
#' which must be greater than 1.
#' @param xmin A numeric value specifying the lower bound \eqn{x_{\min}}.
#' @param n A numeric value specifying the sample size.
#'
#' @details
#' Let \eqn{I(\alpha)} denote the sample Fisher information of the scaling
#' exponent for the discrete power-law model, given by:
#' \deqn{I\left(\alpha\right) = n\left(\frac{\zeta''(\alpha, x_{\min})}
#' {\zeta(\alpha, x_{\min})} - \left(\frac{\zeta'(\alpha, x_{\min})}
#' {\zeta(\alpha, x_{\min})}\right)^2\right)}
#'
#' @return A numeric value representing the information standard
#' error of the scaling exponent.
#'
#' @references
#' Jerez-Lillo, N., Rodrigues, F. A., Ferreira, P. H., & Ramos, P. L. (2025).
#' Bias Reduction and Goodness-of-Fit in Discrete Power-Law Models.
#' Submitted for publication.
#'
#' @examples
#' se_expected(2.0, 1, 100)
#' se_expected(2.5, 1, 100)
#' se_expected(2.0, 2, 250)
#'
#' @seealso [fisher_penalty_pldis] [se_observed]
#'
#' @export
se_expected <- function(alpha, xmin, n)
  sqrt((exp(fisher_penalty_pldis(alpha, xmin)) ^ 2 * n) ^ -1)
