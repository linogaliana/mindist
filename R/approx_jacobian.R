#' Approximate jacobian matrix for loss function
#'
#' A numerical estimation for the jacobian matrix
#' giving the change of loss function dimensions
#' \eqn{\mathcal{l}(\theta)} for marginal change
#' in each dimension of \eqn{\theta} vector
#'
#' @details
#' Jacobian matrix is not derived from gradient
#'  methods but is numerically approximated using
#'  a small \eqn{h} step (\code{step} argument).
#'
#' Parallel implementation is proposed
#'  but is not efficient for the moment: it is
#'  usually slower than the sequential approach
#'
#' @inheritParams loss_function
#' @inheritParams loss_over_identified
#' @param model_function Function that should be used
#'  to transform \eqn{\theta} parameter into
#'  moment conditions
#' @param step \eqn{h} step to numerically compute
#'  derivative
#'
#' @return A \eqn{K \times K} matrix where
#'  `K` is the number of dimensions. If \eqn{i} represents a row and
#'  \eqn{j} a column, jacobian matrix element is
#'  \deqn{\frac{\delta \epsilon_i}{\delta \beta_j}}
#' @param ... Additional arguments
#'
#'
#' @export


approx_jacobian_epsilon <- function(
  theta,
  model_function,
  step = 1e-6,
  ...){



  # epsilon_f <- function(est_theta, true_theta = c(4,2)){
  #   x <- rnorm(1000, mean = est_theta[1], sd = est_theta[2])
  #   epsilon <- c(mean(x) - true_theta[1], sd(x) - true_theta[2])
  #   return(epsilon)
  # }


  all_theta <- lapply(seq_along(theta), function(i) list(
    c(theta[i]+step, theta[-i]),
    c(theta[i]-step, theta[-i]))
  )

  all_theta <- unlist(all_theta, recursive = FALSE)
  grad <- lapply(all_theta, extract_moment,
                 model_function = model_function,
                 return_moment = TRUE, ...)

  names(grad) <- do.call(
    c,
    lapply(names(theta), function(g) paste0(g, 1:2))
  )



  H <- lapply(names(theta), function(param) ({
    (grad[[paste0(param,"1")]] - grad[[paste0(param,"2")]])/(2*step)
  })
  )

  names(H) <- names(theta)

  data.table::setDT(H)

  N <- length(grad[[1]])

  # message("Z argument missing, using ones")
  #
  # Z <- matrix(1L, nrow = nrow(H), ncol = length(theta))
  #
  # return(t(Z) %*%  as.matrix(H)/N)

  # return(as.matrix(H)/N)

  return(H)

}
