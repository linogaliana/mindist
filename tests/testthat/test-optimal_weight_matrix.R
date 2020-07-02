testthat::context("optimal_weight_matrix W computation")


epsilon <- rnorm(10)


testthat::test_that("Dimensions are ok", {
  testthat::expect_equal(
    dim(optimal_weight_matrix(epsilon)),
    c(10,10)
  )
})

testthat::test_that("eps' W eps = dim(eps)^2", {
  testthat::expect_equal(
    as.numeric(
      epsilon %*% optimal_weight_matrix(epsilon) %*% epsilon
      ),
    length(epsilon)^2
  )
})

