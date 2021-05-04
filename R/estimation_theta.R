#' Structural estimation for life-cycle model over microsimulated data
#'
#' @description
#' User front-end function for indirect inference on microsimulated data.
#' Estimation is performed using simulated method of moments and
#' is a particular case of extremum estimator. Theoretical
#' framework can be found on
#' \insertCite{gourieroux1993simulation;textual}{mindist}.
#'
#' Weight matrix \eqn{W} is either assumed to be identity
#'  (`approach` = \emph{"one_stage"}) or
#'  estimated with the optimal weight matrix as
#'  proposed by \insertCite{hansen1982large;textual}{mindist}.
#'  In the former case, known as feasible GMM,
#'  (`approach` = \emph{"two_stage"})
#'
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @param theta_0 Initial values for \eqn{\theta} parameter. This can be a
#'  named vector
#' @param prediction_function Function to transform `theta` into sample conditions.
#'  See examples.
#' @param objective_function How to transform output from prediction_function into
#'  an objective function. Included for flexibility but not recommanded to change
#'  that default option
#' @param approach Estimation approach. Either *one_step* or *two_step* (default)
#' @param ... Additional arguments that should be used by
#'  `model_function`
#'
#' @examples \dontrun{
#'
#' n <- 10000L
#' ncol <- 3
#' x <- replicate(ncol, rnorm(n))
#' y <- 2*x[,2] + runif(n)
#'
#' df <- data.frame(y = y, x1 = x[,1], x2 = x[,2],
#'                  x3 = x[,3])
#'
#' moment_function <- function(theta, ...){
#'   return(
#'     data.table::data.table(
#'       'epsilon' = as.numeric(df$x1*(df$y - cbind(1L, df$x1) %*% theta))
#'     )
#'   )
#' }
#' objective_function <- function(theta, weights = 1L, return_moment = FALSE, ...){
#'   if (return_moment) return(moment_function(theta, ...))
#'   if (weights == 1L) return(
#'   t(moment_function(theta, ...)$epsilon) %*%
#'      moment_function(theta, ...)$epsilon
#'   )
#'   return(
#'     t(moment_function(theta, ...)$epsilon) %*%
#'       weights %*%
#'       moment_function(theta, ...)$epsilon
#'     )
#' }
#'
#' msm1 <- estimation_theta(theta_0 = c(0, 0),
#'                          model_function = objective_function,
#'                          approach = "one_step")
#'
#' }
#'
#'
#' @importFrom stats optim
#' @importFrom MASS ginv
#' @export
#'

estimation_theta <- function(theta_0,
                             prediction_function = {function(theta) theta},
                             objective_function = loss_function,
                             # optim_args = list(),
                             approach = c("two_step","one_step"),
                             ...){

  approach <- match.arg(approach)

  # if ('method' %in% names(optim_args)){
  #   optim_method <- as.character(optim_args$method)
  # } else{
  #   optim_method <- "Nelder-Mead"
  # }

  # STEP 1: ESTIMATE W WEIGHT MATRIX -----------------------------------------


  # 1.1: theta_{(1)} vector =====================================

  NM_step1 <- optim(fn = objective_function,
                    par = theta_0,
                    prediction_function = prediction_function,
                    ...,
                    weights = 1L,
                    # method = optim_method,
                    # control = optim_args,
                    return_moment = FALSE)



  # # Create the first step vector of interest (theta_1)
  # theta_1 <- theta_0
  # theta_1[] <- NM_step1$par # We don't change names

  # 1.2: W_{(1)} vector ===========================================


  df_moment <- objective_function(theta = NM_step1$`par`,
                                  prediction_function = prediction_function,
                                  ...,
                                  weights = 1L,
                                  return_moment = TRUE)


  W_1 <- optimal_weight_matrix(df_moment[['epsilon']])



  if (approach == "two_step"){

    # STEP 2: ESTIMATE THETA WITH W MATRIX -------------------------


    # Create the initial point for second step =====================


    # DETERMINE THETA_HAT =====================

    # WE START FROM \theta_{(1)} WITH W(\theta_{(1)})
    NM_step2 <- optim(fn = objective_function,
                      par = NM_step1$`par`,
                      prediction_function = prediction_function,
                      ...,
                      weights = W_1,
                      # method = optim_method,
                      # control = optim_args,
                      return_moment = FALSE)




  } else{

    NM_step2 <- NM_step1

  }

  estimator_theta <- NM_step2$`par`


  # FINAL STEP: ESTIMATOR VARIANCE ------------------------------


  Gamma <- approx_jacobian_epsilon(
    theta = estimator_theta,
    model_function = prediction_function,
    ...,
    step = 1e-6)


  # We must use Gamma = Jacobian/N_moments
  # Gamma <- Gamma/N_moments


  df_moment_optim <- prediction_function(theta = estimator_theta,
                                         ...,
                                         weights = W_1,
                                         return_moment = TRUE)

  # CAPTURE ENVIRONMENT
  envir <- prediction_function(theta = estimator_theta,
                               match_call = TRUE,
                               ...,
                               weights = W_1)


  # message("Z argument missing, using ones")

  # Z <- matrix(1L, nrow = nrow(df_moment_optim), ncol = length(estimator_theta))

  # Delta <- optimal_weight_matrix(df_moment_optim[['epsilon']]*Z)
  # Delta <- optimal_weight_matrix(df_moment_optim[['epsilon']])

  args <- list(...)

  if (is.null(args[['Z']])){
    # Formula (9.100) in Davidson & MacKinnon (2004)
    syst <- t(as.matrix(Gamma)) %*% W_1 %*% as.matrix(Gamma)
    matinv <- tryCatch(
      solve(syst),error = function(e){
        message("as.matrix(Gamma)) %*% W_1 %*% as.matrix(Gamma) not invertible. Using Moore-Penrose inverse")
        MASS::ginv(syst)
      }
    )
  } else{
    Z <- args[['Z']]
    matinv <- solve(t(as.matrix(Gamma)) %*% Z %*% W_1 %*% Z %*% as.matrix(Gamma))
  }





  # Vtheta <- matinv %*% t(as.matrix(Gamma)) %*% W_1 %*% Delta %*% W_1 %*% as.matrix(Gamma) %*% matinv
  # Vtheta <- Vtheta/N_moments

  # (eq. 9.100 in Davidson MacKinnon)
  # Vtheta <- matinv
  Vtheta <- nrow(Gamma)*matinv

  se_estimator <- sqrt(diag(Vtheta)) # to match with (6) in stata help
  names(se_estimator) <- names(estimator_theta)

  out <- list(
    "estimates" = list(
      "theta_hat" = estimator_theta,
      "se_theta_hat" = se_estimator,
      "theta_1" = NM_step1$`par`,
      "W_1" = W_1,
      "jacobian" = Gamma,
      "vcov_theta" = Vtheta
    ),
    "moments" = list(
      "moment_first_step" = df_moment,
      "moment_optimum" = df_moment_optim
    ),
    "NelderMead" = list(
      "NM_step1" = NM_step1,
      "NM_step2" = NM_step2
    ),
    "initialParams" = list("theta_0" = theta_0),
    "envir" = envir
  )

  class(out) <- c("mindist", class(out))

  attr(out, "approach") <- approach

  return(
    out
  )
}
