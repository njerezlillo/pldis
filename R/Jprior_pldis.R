#' Jeffreys Prior for the Discrete Power-Law Model
#'
#' Computes the Jeffreys prior for a discrete power-law model, which is a
#' non-informative prior commonly used in Bayesian inference to ensure parameter
#' invariance under reparameterization.
#'
#' @param alpha A numeric value representing the scaling parameter \eqn{\alpha},
#' which must be greater than 1.
#' @param xmin A numeric value specifying the lower bound \eqn{x_{\min}}.
#'
#' @details
#' The Jeffreys prior for the discrete power-law model is derived from its
#' Fisher information. Its computation involves the second derivative of the
#' log-likelihood function with respect to \eqn{\alpha}, while keeping
#' \eqn{x_{\min}} fixed. It is given by:
#' \deqn{\pi_1\left(\alpha\right) \propto \sqrt{\frac{\zeta''(\alpha, x_{\min})}
#' {\zeta(\alpha, x_{\min})} - \left(\frac{\zeta'(\alpha, x_{\min})}
#' {\zeta(\alpha, x_{\min})}\right)^2}}
#'
#' This implementation relies on the Hurwitz zeta function, provided by the
#' VGAM package, to accurately compute the required terms.
#'
#' @return A numeric value representing the Jeffreys prior evaluated at
#' the given scaling parameter.
#'
#' @references
#' Clauset, A., Shalizi, C. R., & Newman, M. E. (2009). Power-law distributions
#' in empirical data. *SIAM Review*, 51(4), 661-703.
#'
#' @importFrom VGAM zeta
#'
#' @export
Jprior_pldis <- function (alpha, xmin){
  if (xmin > 2) {
    va <- seq(1, (xmin - 1), 1)
    deriv1 <- zeta(alpha, deriv = 1) + sum(log(va) / (va ^ alpha))
    deriv2 <- zeta(alpha, deriv = 2) - sum((log(va) ^ 2) / (va ^ alpha))
  } else {
    deriv1 <- zeta(alpha, deriv = 1)
    deriv2 <- zeta(alpha, deriv = 2)
  }
  prior <-
    0.5 * log((deriv2 / zeta(alpha, shift = xmin)) -
                ((deriv1 / zeta(alpha, shift = xmin)) ^ 2))
  return(prior)
}
