#' Efficient weight matrix in GMM framework
#'
#' Implementation of the GMM optimal matrix \insertCite{hansen1982large}{mindist}
#'  in a simulated method of moments framework
#'
#' @param epsilon Vector of difference between sample and simulated moments
#'
#'
#' @details
#' The objective function is:
#' \deqn{\theta = \arg \min_{\theta} \bigg\{ \mathcal{l}(\theta) = \epsilon(\theta)^\intercal W(\theta) \epsilon(\theta) \bigg\}}
#' The optimal weight matrix is unknown, an estimator is
#' \deqn{\widehat{W} \big(\widehat{\theta_{(1)}}\big) = \mathbb{E}\bigg(\varepsilon \left(\widehat{\theta_{(1)}} \right) \varepsilon \left( \widehat{\theta_{(1)}}\right)^\intercal \bigg)}
#' where \eqn{\theta_{(1)}} is a first step GMM estimator
#'
#'
#' @export

optimal_weight_matrix <- function(epsilon){

  temp <- epsilon %*% matrix(1,
                             nrow = 1,
                             ncol = length(epsilon))

  temp <- diag(length(epsilon)) * temp

  f <- crossprod(temp)

  f <- f/length(epsilon)
  W <- solve(f)

  return(W)
}

