testthat::context("Test print method")



n <- 1000L
ncol <- 3

mu <- 2
sd <- 2


# PART 1: ESTIMATE PARAMETERS FROM NORMAL DISTRIBUTION --------------

## EXPERIMENTATION 1

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

msm1 <- estimation_theta(theta_0 = c("mu" = 0, "sigma" = 0.2),
                         # model_function = objective_function,
                         prediction_function = moment_function,
                         approach = "one_step")
msm2 <- estimation_theta(theta_0 = c("mu" = 0, "sigma" = 0.2),
                         # model_function = objective_function,
                         prediction_function = moment_function,
                         approach = "two_step")


msm1_print <- capture.output(msm1)
msm2_print <- capture.output(msm2)


testthat::test_that("Header are correct",{
  testthat::expect_equal(
    trimws(msm1_print[1:2]),
    c("Minimum distance estimation", "1 step estimation:")
  )
  testthat::expect_equal(
    trimws(msm2_print[1:2]),
    c("Minimum distance estimation", "2 steps estimation:")
  )
})


