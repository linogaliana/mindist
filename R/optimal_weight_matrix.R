#' Efficient weight matrix in GMM framework
#'
#' Implementation of the GMM optimal matrix \insertCite{hansen1982large}{mindist}
#'  in a simulated method of moments framework
#'
#' @param epsilon Vector of difference between sample and simulated moments
#'
#'
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

