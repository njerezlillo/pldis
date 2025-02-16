% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglik_pldis.R
\name{loglik_pldis}
\alias{loglik_pldis}
\title{Log-Likelihood of the Discrete Power-Law Model}
\usage{
loglik_pldis(alpha, xmin, x)
}
\arguments{
\item{alpha}{A numeric value representing the scaling parameter \eqn{\alpha},
which must be greater than 1.}

\item{xmin}{A numeric value specifying the lower bound \eqn{x_{\min}}.}

\item{x}{A numerical vector of observed data.}
}
\value{
A numeric value representing the log-likelihood of the data
given the specified parameters.
}
\description{
Computes the log-likelihood of a discrete power-law model given a lower bound,
scaling parameter, and observed dataset.
}
\details{
The log-likelihood function of the discrete power-law model is given by:
\deqn{\mathcal{L}(\alpha, x_{\min}) = -n\log\zeta(\alpha, x_{\min}) -
\alpha\sum_{i=1}^{n}\log(x_i).}

The computation relies on the Hurwitz zeta function from the
VGAM package.
}
\examples{
# Example data
x <- c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 10)

# Evaluating the log-likelihood for different parameter values
loglik_pldis(2.0, 1, x)
loglik_pldis(2.5, 1, x)
loglik_pldis(2.0, 2, x)
loglik_pldis(2.5, 2, x)

}
