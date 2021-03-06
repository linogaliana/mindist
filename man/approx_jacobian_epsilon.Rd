% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/approx_jacobian.R
\name{approx_jacobian_epsilon}
\alias{approx_jacobian_epsilon}
\title{Approximate jacobian matrix for loss function}
\usage{
approx_jacobian_epsilon(theta, model_function, step = 1e-06, ...)
}
\arguments{
\item{theta}{Vector of structural parameters. Assuming a named vector.}

\item{model_function}{Function that should be used
to transform \eqn{\theta} parameter into
moment conditions}

\item{step}{\eqn{h} step to numerically compute
derivative}

\item{...}{Additional arguments}
}
\description{
A numerical estimation for the jacobian matrix
giving the change of loss function dimensions
\eqn{\mathcal{l}(\theta)} for marginal change
in each dimension of \eqn{\theta} vector
}
\details{
Jacobian matrix is not derived from gradient
methods but is numerically approximated using
a small \eqn{h} step (\code{step} argument).

Parallel implementation is proposed
but is not efficient for the moment: it is
usually slower than the sequential approach
}
