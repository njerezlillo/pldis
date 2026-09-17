#' Discrete Fisher-Type Penalty for the Power-Law Model
#'
#' Computes the Firth-type bias-reduction penalty for the discrete power-law
#' model, derived from its Fisher information. Added to the ordinary
#' log-likelihood, it defines the discrete Fisher-penalized log-likelihood
#' used to obtain a nearly unbiased estimator of the scaling exponent
#' (see [loglik_fd_pldis]).
#'
#' @param alpha A numeric value representing the scaling parameter \eqn{\alpha},
#' which must be greater than 1.
#' @param xmin A numeric value specifying the lower bound \eqn{x_{\min}}.
#'
#' @details
#' Following Firth (1993), a bias-reducing penalty for a one-parameter model
#' is one half of the log of the Fisher information. For the discrete
#' power-law model this penalty is:
#' \deqn{\pi_1\left(\alpha\right) = \frac{1}{2}\log\left(\frac{\zeta''(\alpha, x_{\min})}
#' {\zeta(\alpha, x_{\min})} - \left(\frac{\zeta'(\alpha, x_{\min})}
#' {\zeta(\alpha, x_{\min})}\right)^2\right)}
#'
#' This implementation relies on the Hurwitz zeta function, provided by the
#' VGAM package, to accurately compute the required terms.
#'
#' @return A numeric value representing the Fisher-type penalty evaluated at
#' the given scaling parameter.
#'
#' @references
#' Firth, D. (1993). Bias reduction of maximum likelihood estimates.
#' *Biometrika*, 80(1), 27-38.
#'
#' Jerez-Lillo, N., Rodrigues, F. A., Ferreira, P. H., & Ramos, P. L. (2026).
#' Bias Reduction and Goodness-of-Fit in Discrete Power-Law Models.
#' Submitted for publication.
#'
#' @examples
#' fisher_penalty_pldis(2.0, 1)
#' fisher_penalty_pldis(2.5, 1)
#' fisher_penalty_pldis(2.0, 2)
#' fisher_penalty_pldis(2.5, 2)
#'
#' @importFrom VGAM zeta
#'
#' @export
fisher_penalty_pldis <- function (alpha, xmin){
  if (xmin > 2) {
    va <- seq(1, (xmin - 1), 1)
    deriv1 <- zeta(alpha, deriv = 1) + sum(log(va) / (va ^ alpha))
    deriv2 <- zeta(alpha, deriv = 2) - sum((log(va) ^ 2) / (va ^ alpha))
  } else {
    deriv1 <- zeta(alpha, deriv = 1)
    deriv2 <- zeta(alpha, deriv = 2)
  }
  penalty <-
    0.5 * log((deriv2 / zeta(alpha, shift = xmin)) -
                ((deriv1 / zeta(alpha, shift = xmin)) ^ 2))
  return(penalty)
}
