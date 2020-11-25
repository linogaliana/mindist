calibration_theta <- function(theta,
                             model_function = {function(theta) theta},
                             # optim_args = list(),
                             approach = c("two_step","one_step"),
                             ...){

  approach <- match.arg(approach)

  df_moment <- model_function(theta = theta,
                              ...,
                              weights = 1L,
                              return_moment = TRUE)


  W_1 <- optimal_weight_matrix(df_moment[['epsilon']])



  estimator_theta <- theta


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
        "W_1" = W_1,
        "jacobian" = Gamma,
        "vcov_theta" = Vtheta
      ),
      "moments" = list(
        "moment_first_step" = df_moment,
        "moment_optimum" = df_moment_optim
      ),
      "envir" = envir
    )
  )
}
