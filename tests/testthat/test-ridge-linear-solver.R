
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
})
  


X <- matrix(rnorm(20), nrow = 5)
Y <- matrix(rnorm(20), nrow = 5)

# valid call to ensure baseline
expect_silent(fmriparametric:::.ridge_linear_solve(X, Y, lambda_ridge = 0))

test_that("negative lambda triggers error", {
  expect_error(
    fmriparametric:::.ridge_linear_solve(X, Y, lambda_ridge = -1),
    "lambda_ridge"
  )
})

test_that("NA lambda triggers error", {
  expect_error(
    fmriparametric:::.ridge_linear_solve(X, Y, lambda_ridge = NA_real_),
    "lambda_ridge"

  )
})
