library(fmriparametric)
library(testthat)

context(".ridge_linear_solve input validation")

# Non-finite values in X should trigger an error

test_that("non-finite values in X cause error", {
  X <- matrix(c(1, NA, 3, 4), nrow = 2)
  Y <- matrix(rnorm(4), nrow = 2)
  expect_error(
    fmriparametric:::.ridge_linear_solve(X, Y, 0),
    "X contains non-finite values"
  )
})

# Non-finite values in Y should trigger an error

test_that("non-finite values in Y cause error", {
  X <- matrix(rnorm(4), nrow = 2)
  Y <- matrix(c(1, Inf, 3, 4), nrow = 2)
  expect_error(
    fmriparametric:::.ridge_linear_solve(X, Y, 0),
    "Y contains non-finite values"
  )
})
