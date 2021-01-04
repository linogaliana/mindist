#' Compute loss function in a flexible way
#'
#' Compute loss associated with a vector of structural parameters
#'  between some data and a transformation of structural parameters
#'
#' @param theta Vector of structural parameters. Assuming a named vector.
#' @param prediction_function Function that transforms \eqn{\theta}
#'  into vector of moments. In Newer and MacFadden, denoted \eqn{g(\theta)}
#' @param weights Weight matrix \eqn{W} that should be used
#' @param moments_weights User-defined weights that should be applied
#'  to reweight moments importance. This is a user choice
#'  independent of using an optimal weight matrix
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
#'
#' \link{estimation_theta} for GMM estimation
#'
#' @export


loss_function <- function(theta,
                          prediction_function = {function(theta) theta},
                          weights = 1L,
                          verbose = FALSE,
                          return_moment = FALSE,
                          moments_weights = NULL,
                          ...
){

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

  if (!inherits(df_moment, "data.table")) stop("prediction_function should return a data.table with a column epsilon (e.g. diff between observed and predicted values) that aims to be minimized")


  # COMPUTE LOSS ---------------------------------------

  # RESIDUAL ERROR
  epsilon <- df_moment[['epsilon']]

  if (!is.null(moments_weights)){
    if (is.character(moments_weights)){
      print("moments_weights is provided as character, assuming this is a variable name")
      moments_weights <- df_moment[[moments_weights]]
    }
    # epsilon <- epsilon*moments_weights
  }


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

  eps_weight <- epsilon
  if (!is.null(moments_weights)){
    eps_weight <- sqrt(moments_weights)*eps_weight
  }

  if (length(weights) == 1){
    # with default weights = 1L,
    # equivalent to sse_hat = e_adj' * W_0 * e_adj in Einav et al. replication code in objective_function.m
    return(sum(eps_weight^2)*weights)
  }

  # we mimic Einav et al. replication code (compute_std_errors.m)
  return(
    as.numeric(
      t(eps_weight) %*% W %*% eps_weight
    )
  )



}
