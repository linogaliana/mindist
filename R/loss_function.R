#' Compute loss function in a flexible way
#'
#' Compute loss associated with a vector of structural parameters
#'  between some data and a transformation of structural parameters
#'
#' @param theta Vector of structural parameters. Assuming a named vector.
#' @param prediction_function Function that transforms \eqn{\theta}
#'  into vector of moments. In Newer and MacFadden, denoted \eqn{g(\theta)}
#' @param ... Additional arguments that should be used to control
#'  `prediction_function` behavior. This function should return a
#'  `data.table` object with a variable denoted `epsilon` giving the
#'  distance that should be minimized
#' @param verbose Logical value indicating whether we want to print informative
#'  messages to monitor progress
#' @param return_moment Logical value indicating whether we want to
#'  return moment matrix (observed, simulation, difference)
#' @return A distance between empirical and simulated moments
#'
#'
#' @return Denoting \eqn{\epsilon(\theta)} the difference between
#'  population and simulated moments, loss function (or
#'  objective function) writes down
#'  \deqn{ L(\theta) = \epsilon(\theta)'  W(\theta) \epsilon(\theta)}
#' \eqn{W(\theta)} is determined by \code{weight} argument
#'
#' \link{estimation_theta} for GMM estimation
#'
#' @export


loss_function <- function(theta,
                          prediction_function = {function(theta) theta},
                          weights = 1L,
                          verbose = FALSE,
                          return_moment = FALSE,
                          ...
){

  loss_scale <- match.arg(loss_scale)
  args <- list(...)

  if ('match_call' %in% names(args)){
    envir <- parent.frame()
    return(envir)
  }


  # COMPUTE MOMENT GIVEN THETA PARAMETER ----------------------------


  df_moment <- prediction_function(theta = theta, ...)

  if (return_moment) return(df_moment)


  if (verbose) cat(
    paste0("\n======= Parameter values: ========== \n ", paste(paste0(names(theta), ": "), theta, collapse = "; "), "\n")
  )



  # COMPUTE LOSS ---------------------------------------

  # RESIDUAL ERROR
  epsilon <- df_moment[['epsilon']]


  if (verbose) print(epsilon)

  # If weights not provided: W is (eps*eps')^{-1}
  if (is.null(weights)){
    W <- optimal_weight_matrix(epsilon)
  } else{
    W <- weights
  }

  if (verbose){
    print(df_moment)
    print("Weight matrix W: ")
    print(W)
  }


# LOSS FUNCTION COMPUTATIONS ------------------

  return(
    as.numeric(
      t(epsilon) %*% W %*% epsilon
    )
  )


}
