
#' Overidentification test for method of simulated moments
#'
#' @param theta_hat Estimated theta parameter
#' @param S Estimated covariance matrix
#' @param model_function \eqn{g(\theta)} transformation
#'  defining moment conditions
#' @inheritParams estimation_theta
#' @param ... Additional arguments that should be passed
#'  to `model_function`
#'
#' Under standard hypotheses, GMM yields
#'  \deqn{n*g(\theta)^T*W*g(\theta)}
#' @return List with following elements
#' @importFrom stats pchisq
#' @export

sargan_test <- function(theta_hat, S,
                        model_function = {function(theta) theta},
                        ...){


  g_theta <- model_function(theta = theta_hat,
                            ...,
                            weights = S,
                            return_moment = TRUE)

  g_theta <- g_theta[['epsilon']]

  J <- nrow(S) *( g_theta %*% S %*% g_theta )

  overid_dim <- nrow(S) - length(theta_hat)

  pvalue <- stats::pchisq(J, df = overid_dim, lower.tail = FALSE)


  return(list("J" = J, "pvalue" = pvalue))
}
