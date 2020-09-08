context("test-estimation_theta")


n <- 1000L
ncol <- 3

mu <- 2
sd <- 2


# PART 1: ESTIMATE PARAMETERS FROM NORMAL DISTRIBUTION --------------


x <- rnorm(n, mu, sd)
theta_observed <- c(mean(x),sd(x))

moment_function <- function(theta, ...){
  x_hat <- rnorm(n, mean = theta[1], sd = theta[2])
  theta_simul <- c(mean(x_hat),sd(x_hat))
  return(
    data.table::data.table(
      'epsilon' = theta_observed - theta_simul
    )
  )
}

objective_function <- function(theta, prediction_function, weights = 1L, return_moment = FALSE, ...){
  if (return_moment) return(moment_function(theta))
  if (length(weights)==1) return( t(moment_function(theta, ...)$epsilon) %*% moment_function(theta, ...)$epsilon)
  eps <- moment_function(theta, ...)$epsilon
  return(t(eps) %*% weights %*% eps)
}

msm1 <- estimation_theta(theta_0 = c("mu" = 0, "sigma" = 0.2),
                         model_function = objective_function,
                         prediction_function = moment_function,
                         approach = "one_step")

msm2 <- estimation_theta(theta_0 = c("mu" = 0, "sigma" = 0.2),
                         model_function = objective_function,
                         prediction_function = moment_function,
                         approach = "two_step")


test_that("Method of simulated moments should be close from theoretical parameters", ({
  expect_equal(c(mu,sd), as.numeric(msm1$NelderMead$`NM_step1`$`par`), tolerance = 10e-1)
}))

test_that("Method of simulated moments should be close from theoretical parameters", ({
  expect_equal(c(mu,sd), as.numeric(msm2$NelderMead$`NM_step1`$`par`), tolerance = 10e-1)
}))



# PART 2: REPLICATE OLS --------------------------------


# x <- replicate(ncol, rnorm(n))
#
# df <- data.frame(x1 = x[,1], x2 = x[,2],
#                  x3 = x[,3])
#
# df$y <- 1 + 2*df$x1 + rnorm(n)
#
#
# prediction_function <- function(theta){
#   return(
#     cbind(1L, df$x1) %*% theta
#   )
# }
#
# # FORMALISM REQUIRED FOR OUR FUNCTIONS
# moment_function <- function(theta){
#   return(
#     data.table::data.table(
#       'epsilon' = as.numeric(df$x1*(df$y - cbind(1L, df$x1) %*% theta))
#     )
#   )
# }
#
#
# objective_function <- function(theta, prediction_function, weights = 1L, return_moment = FALSE, ...){
#   if (return_moment) return(moment_function(theta))
#   if (length(weights)==1) return(
#     moment_function(theta)$epsilon*moment_function(theta)$epsilon
#     )
#   eps <- moment_function(theta)
# #epsilon/as.numeric(cbind(1L, df$x1) %*% theta)
#   return(t(eps) %*% weights %*% eps)
# }
#
#
# # 1: LINEAR REGRESSION ===========
#
# ols <- lm(y~ x1, data = df)
#
# requireNamespace("gmm", quietly = TRUE)
#
#
# # 2a : GMM PACKAGE
# gmm1 <- gmm::gmm(y~ x1, x = ~x1, wmatrix = "ident", data = df,
#                  vcov="TrueFixed", weightsMatrix = diag(2))
# gmm2 <- gmm::gmm(y~ x1, x = ~x1, data = df)
#
#
#
#
# msm1 <- estimation_theta(theta_0 = c("const" = 0.1, "beta1" = 0),
#                          model_function = objective_function,
#                          prediction_function = moment_function,
#                          approach = "one_step")
#
#
# msm2 <- estimation_theta(theta_0 = c("const" = 0, "beta1" = 0),
#                          model_function = objective_function,
#                          prediction_function = moment_function,
#                          approach = "two_step")
#
#
# test_that("[One-step] Method of simulated moments should be close from OLS estimator using E(xu) condition", ({
#   expect_equal(as.numeric(ols$coefficients), as.numeric(msm1$NelderMead$`NM_step1`$`par`), tolerance = 10e-2)
# }))
#
# test_that("[One-step] Method of simulated moments should be close from theoretical GMM estimator", ({
#   expect_equal(as.numeric(gmm1$coefficients), as.numeric(msm1$NelderMead$`NM_step1`$`par`), tolerance = 10e-2)
# }))
#
#
# test_that("[Two-steps] Method of simulated moments should be close from OLS estimator using E(xu) condition", ({
#   expect_equal(as.numeric(ols$coefficients), as.numeric(msm2$estimates$theta_hat), tolerance = 10e-2)
# }))
#
# test_that("[Two-steps] Method of simulated moments should be close from theoretical GMM estimator", ({
#   expect_equal(as.numeric(gmm1$coefficients), as.numeric(msm2$estimates$theta_hat), tolerance = 10e-2)
# }))
#
#
#
#
# se_ols <- summary(ols)$coefficients[,'Std. Error']
# se_gmm1 <- summary(gmm1)$coefficients[,'Std. Error']
# se_gmm2 <- summary(gmm2)$coefficients[,'Std. Error']
# se_msm1 <- msm1$estimates$se_theta_hat
# se_msm2 <- msm2$estimates$se_theta_hat


# test_that("[Two-steps ; s.e.] Method of simulated moments should be close from OLS estimator using E(xu) condition", ({
#   expect_equal(as.numeric(se_ols), as.numeric(se_msm2), tolerance = 10e-2)
# }))
#
# test_that("[Two-steps ; s.e.] Method of simulated moments should be close from theoretical GMM estimator", ({
#   expect_equal(as.numeric(se_gmm2), as.numeric(se_msm2), tolerance = 10e-2)
# }))
