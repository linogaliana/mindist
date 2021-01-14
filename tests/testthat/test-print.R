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


testthat::test_that("Headers are correct",{
  testthat::expect_equal(
    trimws(msm1_print[1:2]),
    c("Minimum distance estimation", "1 step estimation:")
  )
  testthat::expect_equal(
    trimws(msm2_print[1:2]),
    c("Minimum distance estimation", "2 steps estimation:")
  )
})

testthat::test_that("Number parameters correct",{

  testthat::expect_equal(
    as.numeric(trimws(
      gsub("Number of parameters:", "",
           msm1_print[grep("Number of parameters:", msm1_print)]
      )
    )),
    length(msm1$estimates$theta_hat)
  )

  testthat::expect_equal(
    as.numeric(trimws(
      gsub("Number of parameters:", "",
           msm1_print[grep("Number of parameters:", msm2_print)]
      )
    )),
    length(msm2$estimates$theta_hat)
  )


})

testthat::test_that("Number moments correct",{

  testthat::expect_equal(
    as.numeric(trimws(
      gsub("Number of moments:", "",
           msm1_print[grep("Number of moments:", msm1_print)]
      )
    )),
    nrow(msm1$moments$moment_optimum)
  )

  testthat::expect_equal(
    as.numeric(trimws(
      gsub("Number of moments:", "",
           msm1_print[grep("Number of moments:", msm2_print)]
      )
    )),
    nrow(msm2$moments$moment_optimum)
  )

  testthat::expect_equal(
    length(msm1_print[msm1_print == "(Exact identification)"]),
    1
  )
  testthat::expect_equal(
    length(msm2_print[msm2_print == "(Exact identification)"]),
    1
  )

})

digits <- max(3L, getOption("digits") - 3L)

testthat::test_that("E(\\epsilon'\\epsilon) correct",{
  testthat::expect_equal(
    sum(grepl(
      sprintf("Step 1: %s",
              format(
                mean(
                  msm1$moments$moment_first_step$epsilon^2
                  ), digits = digits
              )
      ), msm1_print)
    ),
    1
  )
  testthat::expect_equal(
    sum(grepl(
      sprintf("Step 1: %s",
              format(
                mean(
                  msm2$moments$moment_first_step$epsilon^2
                ), digits = digits
              )
      ), msm2_print)
    ),
    1
  )
  testthat::expect_equal(
    sum(grepl(
      sprintf("Step 2: %s",
              format(
                mean(
                  msm2$moments$moment_optimum$epsilon^2
                ), digits = digits
              )
      ), msm2_print)
    ),
    1
  )

})
