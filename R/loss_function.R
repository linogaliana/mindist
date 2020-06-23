#' Compute loss function between simulated and empirical moments
#'
#' Compute loss associated with a vector of structural parameters
#'  between empirical \eqn{(m_k)_k} and simulated
#'  moments \eqn{(\hat{m}_k)_k}
#' @inheritParams create_moment_data
#' @param theta \eqn{theta} parameter
#' @param prediction_function Function that transforms \eqn{\theta}
#'  into vector of moments
#' @param ... Additional arguments that should be used to control
#'  `prediction_function` behavior
#' @param weights GMM weight matrix that should be used to rescale
#'  loss function
#' @param dim Dimension for the loss function. If equal to 1 (default),
#'  loss function is real valued and \link{optim} can be called, otherwise
#'  optimization should be performed with other function
#' @param loss_scale Should the loss be based on difference in levels (\emph{level}),
#'   in log (\emph{log}) or in euclidian norm (\emph{square}) between empirical and
#'  simulated. Ignored if \code{dim=1}
#' @param verbose Logical value indicating whether we want to print informative
#'  messages to monitor progress
#' @param return_moment Logical value indicating whether we want to
#'  return moment matrix (observed, simulation, difference)
#' @return A distance between empirical and simulated moments,
#'  based on the \code{loss_scale} and \code{dim} parameters
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
                          dim = 1L,
                          loss_scale = c("level","log","square"),
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


  # UNWEIGHTED GMM ESTIMATES =============

  if (length(weights) == 1L){

    if (verbose) print(df_moment)
    W <- diag(length(epsilon))

  } else{
    # WEIGHTED GMM ESTIMATES =============

    print(epsilon)

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
  }

  # LOSS FUNCTION COMPUTATIONS ------------------

  if (dim==1){
    return(
      as.numeric(
        t(epsilon) %*% W %*% epsilon
      )
    )
  } else{

    if (loss_scale == "square"){
      df_moment[,'l' := sqrt((get('moment_simulations')-get('moment_data'))^2)]
    } else if (loss_scale == "level"){
      df_moment[,'l' := get('moment_simulations')-get('moment_data')]
    } else{
      df_moment[,'l' := abs(get('moment_simulations')-get('moment_data'))]
      df_moment[get('Nmoment')==2, c('l') := log(get('l'))]
    }

    return(
      df_moment[['l']]
    )

  }

}
