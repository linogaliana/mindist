#' Structural estimation for life-cycle model over microsimulated data
#'
#' @description
#' User front-end function for indirect inference on microsimulated data.
#' Estimation is performed using simulated method of moments and
#' is a particular case of extremum estimator. Theoretical
#' framework can be found on
#' \insertCite{gourieroux1993simulation;textual}{wealthyR}.
#'
#' Weight matrix \eqn{W} is either assumed to be identity
#'  (`method` = \emph{"one_stage"}) or
#'  estimated with the optimal weight matrix as
#'  proposed by \insertCite{hansen1982large;textual}{wealthyR}.
#'  In the former case, known as feasible GMM,
#'  (`method` = \emph{"two_stage"})
#'
#' @details
#' The idea in GMM framework is to find \eqn{\theta} estimator such that
#' \deqn{\hat\theta=\argmin_\theta
#'   \left(\mathcal{K}^o-\mathcal{K}^s(\theta)\right)^\intercal
#'   W^{-1}\left(\mathcal{K}^o-\mathcal{K}^s(\theta)\right)}
#' In practice, we use Nelder-Mead algorithm (see \link[stats]{optim}) to converge
#' to an estimates for \eqn{\theta}.
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @param theta_0 Initial values for \eqn{\theta} parameter. This can be a
#'  named vector
#' @param model_function Function to transform `theta` into sample conditions.
#'  See examples.
#' @param optim_args Arguments that should be used to control \link[stats]{optim}
#'  routine
#' @param method Estimation method. Either *one_step* or *two_step* (default)
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
#'                          method = "one_step")
#'
#' }
#'
#'
#' @importFrom stats optim
#' @export
#'

estimation_theta <- function(theta_0,
                             model_function = {function(theta) theta},
                             optim_args = list(),
                             method = c("two_step","one_step"),
                             ...){

  method <- match.arg(method)

  # STEP 1: ESTIMATE W WEIGHT MATRIX -----------------------------------------


  # 1.1: theta_{(1)} vector =====================================

  NM_step1 <- optim(fn = model_function,
                    par = theta_0,
                    ...,
                    weights = 1L,
                    method = "Nelder-Mead",
                    control = optim_args,
                    return_moment = FALSE)



  # # Create the first step vector of interest (theta_1)
  # theta_1 <- theta_0
  # theta_1[] <- NM_step1$par # We don't change names

  # 1.2: W_{(1)} vector ===========================================


  df_moment <- model_function(theta = NM_step1$`par`,
                              ...,
                              weights = 1L,
                              return_moment = TRUE)


  W_1 <- optimal_weight_matrix(df_moment[['epsilon']])



  if (method == "two_step"){

    # STEP 2: ESTIMATE THETA WITH W MATRIX -------------------------


    # Create the initial point for second step =====================


    # DETERMINE THETA_HAT =====================

    # WE START FROM \theta_{(1)} WITH W(\theta_{(1)})
    NM_step2 <- optim(fn = model_function,
                      par = NM_step1$`par`,
                      ...,
                      weights = W_1,
                      method = "Nelder-Mead",
                      control = optim_args,
                      return_moment = FALSE)




  } else{

    NM_step2 <- NM_step1

  }

  estimator_theta <- NM_step2$`par`


  # FINAL STEP: ESTIMATOR VARIANCE ------------------------------


  Gamma <- approx_jacobian_epsilon(
    theta = estimator_theta,
    model_function = model_function,
    ...,
    step = 1e-6)


  # We must use Gamma = Jacobian/N_moments
  # Gamma <- Gamma/N_moments


  df_moment_optim <- model_function(theta = estimator_theta,
                                    ...,
                                    weights = W_1,
                                    return_moment = TRUE)

  # CAPTURE ENVIRONMENT
  envir <- model_function(theta = estimator_theta,
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
    matinv <- solve(t(as.matrix(Gamma)) %*% W_1 %*% as.matrix(Gamma))
  } else{
    Z <- args[['Z']]
    matinv <- solve(t(as.matrix(Gamma)) %*% Z %*% W_1 %*% Z %*% as.matrix(Gamma))
  }




  # Vtheta <- matinv %*% t(as.matrix(Gamma)) %*% W_1 %*% Delta %*% W_1 %*% as.matrix(Gamma) %*% matinv
  # Vtheta <- Vtheta/N_moments

  # (eq. 9.100 in Davidson MacKinnon)
  Vtheta <- matinv
  Vtheta <- matinv
  # Vtheta <- matinv*nrow(Gamma)

  se_estimator <- sqrt(diag(Vtheta))/nrow(Gamma) # to match with (9.100) in MacKinnon (because our Vtheta is not 1/n factor)
  names(se_estimator) <- names(estimator_theta)


  return(
    list(
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
  )
}
