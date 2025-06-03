library(testthat)
library(fmriparametric)

context(".ridge_linear_solve input validation")

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
