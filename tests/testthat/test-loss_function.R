context("test-loss_function")


# play with prediction function --------------------------

testthat::test_that(
  "By default, prediction_function is identity",
  {
    testthat::expect_identical(
      loss_function(2L, return_moment = TRUE),
      2L
    )
    testthat::expect_identical(
      loss_function(seq_len(10), return_moment = TRUE),
      seq_len(10)
    )
  }
)


testthat::test_that(
  "prediction_function can be changed",
  {
    testthat::expect_equal(
      loss_function(2L, return_moment = TRUE, prediction_function = function(theta) return(theta^2)),
      4L
    )
    testthat::expect_equal(
      loss_function(seq_len(10), return_moment = TRUE, prediction_function = function(theta) return(sqrt(theta))),
      sqrt(seq_len(10))
    )
    testthat::expect_equal(
      loss_function(seq_len(10), return_moment = TRUE, prediction_function = function(theta, a) return(a*theta), a = -3),
      -3*seq_len(10)
    )
    testthat::expect_equal(
      loss_function(seq_len(10), return_moment = TRUE, prediction_function = function(theta, a) return(data.frame(x = a*theta, y = a)), a = -3),
      data.frame(x = -3*seq_len(10), y = -3)
    )
  }
)


testthat::test_that(
  "error when prediction_function does not depend on theta",
  testthat::expect_error(
    loss_function(2L, return_moment = TRUE, prediction_function = function(x) return(x^2))
  )
)

testthat::test_that(
  "error when prediction_function does not return a data.table with epsilon",
  testthat::expect_error(
    loss_function(2L, prediction_function = function(theta) return(theta^2))
  )
)


# simplest setup: no weight or moment_weights -------------------

testthat::test_that(
  "In that case, l(theta) = mean(prediction_function(theta)^2)",
  testthat::expect_equal(
    loss_function(c(2L,3L), prediction_function = function(theta) return(data.table::data.table(epsilon = theta))),
    sum(c(2L,3L)^2)/4 #4: N^2
  )
)

testthat::test_that(
  "Even with more complex setup, l(theta) = sum(prediction_function(theta)^2)",{
    testthat::expect_equal(
      loss_function(c(2L,3L), prediction_function = function(theta) return(data.table::data.table(epsilon = theta + 2))),
      sum((c(2L,3L) + 2)^2)/4
    )
    testthat::expect_equal(
      loss_function(c(2L,3L), prediction_function = function(theta, a) return(data.table::data.table(epsilon = theta + a*log(abs(theta)))), a = 2),
      sum((c(2L,3L) + 2*log(c(2L,3L)))^2)/4
    )
  }
)


# add weight argument ---------------------

# weight is scalar
testthat::test_that(
  "scalar: just a multiplicative effect",
  {
    testthat::expect_equal(
      loss_function(c(2L,3L), prediction_function = function(theta) return(data.table::data.table(epsilon = theta)), weights = 2L),
      sum(2*c(2L,3L)^2)/4
    )
    testthat::expect_equal(
      loss_function(c(2L,3L), prediction_function = function(theta) return(data.table::data.table(epsilon = theta)), weights = log(2)),
      sum(log(2)*c(2L,3L)^2)/4
    )
  }
)

# weight is a matrix
testthat::test_that(
  "matrix: \\eps' %*% W %*% \\eps",
  {
    testthat::expect_equal(
      loss_function(c(2L,3L), prediction_function = function(theta) return(data.table::data.table(epsilon = theta)), weights = diag(2)),
      sum(c(2L,3L)^2)/4
    )
    testthat::expect_equal(
      loss_function(c(2L,3L), prediction_function = function(theta) return(data.table::data.table(epsilon = theta)), weights = diag(rep(2,2))),
      sum(2*c(2L,3L)^2)/4
    )
    testthat::expect_equal(
      loss_function(c(2L,3L), prediction_function = function(theta) return(data.table::data.table(epsilon = theta)), weights = diag(c(1,2))),
      as.numeric(c(2,3) %*% diag(c(1,2)) %*% c(2,3))/4
    )
  }
)

# weight is NULL: W = (eps'eps)^{-1} -> l(theta) = length(theta)^2
testthat::test_that(
  "matrix: l(theta) = length(theta)^2",
  {
    testthat::expect_equal(
      loss_function(seq_len(10), prediction_function = function(theta) return(data.table::data.table(epsilon = theta)), weights = NULL),
      100/10^2
    )
    # prediction_function does not matter, arithmetic equality
    testthat::expect_equal(
      loss_function(seq_len(10), prediction_function = function(theta) return(data.table::data.table(epsilon = theta + log(theta))), weights = NULL),
      100/10^2
    )
    testthat::expect_equal(
      loss_function(seq_len(13), prediction_function = function(theta) return(data.table::data.table(epsilon = theta + log(theta))), weights = NULL),
      13^2/13^2
    )
  }
)


# add weight_moments argument ---------------------

# weight_moments are multiplicative
testthat::test_that(
  "[without weights] l(theta) = sum(weight_moments*prediction_function(theta)^2)",
  {
    testthat::expect_equal(
      loss_function(c(2L,3L), prediction_function = function(theta) return(data.table::data.table(epsilon = theta)),
                    moments_weights = 4L
                    ),
      sum(4L*c(2L,3L)^2)/2^2
    )

    testthat::expect_equal(
      loss_function(c(2L,3L), prediction_function = function(theta) return(data.table::data.table(epsilon = theta)),
                    moments_weights = c(1L,4L)
                    ),
      sum(c(1L,4L)*c(2L,3L)^2)/2^2
    )
  }
)

testthat::test_that(
  "[without weights] moments_weighting_formula control how weighting occurs",
  {
    testthat::expect_equal(
      loss_function(c(2L,3L), prediction_function = function(theta) return(data.table::data.table(epsilon = theta)),
                    moments_weights = 4L,
                    moments_weighting_formula = as.formula("w ~ I(moments_weights^2)")),
      sum(4L^2*c(2L,3L)^2)/2^2
    )

    testthat::expect_equal(
      loss_function(c(2L,3L), prediction_function = function(theta) return(data.table::data.table(epsilon = theta)),
                    moments_weights = c(1L,4L),
                    moments_weighting_formula = as.formula("w ~ I(log(moments_weights))")),
      sum(log(c(1L,4L))*c(2L,3L)^2)/2^2
    )
  }
)


testthat::test_that(
  "[with weights = NULL] no longer l(theta) = length(theta)^2",
  testthat::expect_equal(
    loss_function(seq_len(10),
                  prediction_function = function(theta) return(data.table::data.table(epsilon = theta)),
                  weights = NULL,
                  moments_weights = 1/seq_len(10)),
    as.numeric(t(sqrt(1/seq_len(10)) * seq_len(10)) %*% optimal_weight_matrix(seq_len(10)) %*% (sqrt(1/seq_len(10)) * seq_len(10)))/10^2
  )
)


testthat::test_that(
  "[with weights as matrix]",
  {
    testthat::expect_equal(
      loss_function(seq_len(10),
                    prediction_function = function(theta) return(data.table::data.table(epsilon = theta)),
                    weights = diag(10),
                    moments_weights = 1/seq_len(10)),
      as.numeric(t(sqrt(1/seq_len(10)) * seq_len(10)) %*% diag(10) %*% (sqrt(1/seq_len(10)) * seq_len(10)))/10^2
    )
    testthat::expect_equal(
      loss_function(seq_len(10),
                    prediction_function = function(theta) return(data.table::data.table(epsilon = theta)),
                    weights = diag(seq_len(10)),
                    moments_weights = 1/seq_len(10)),
      as.numeric(t(sqrt(1/seq_len(10)) * seq_len(10)) %*% diag(seq_len(10)) %*% (sqrt(1/seq_len(10)) * seq_len(10)))/10^2
    )
    testthat::expect_equal(
      loss_function(seq_len(10),
                    prediction_function = function(theta) return(data.table::data.table(epsilon = theta)),
                    weights = diag(seq_len(10)),
                    moments_weights = seq_len(10),
                    moments_weighting_formula = as.formula("w ~ I(log(moments_weights))")
                    ),
      as.numeric(t(sqrt(log(seq_len(10))) * seq_len(10)) %*% diag(seq_len(10)) %*% (sqrt(log(seq_len(10))) * seq_len(10)))/10^2
    )
  }
)

