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

  f <- diag(epsilon^2)

  f <- f/length(epsilon)
  W <- solve(f)

  return(W)
}

