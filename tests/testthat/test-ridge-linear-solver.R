library(fmriparametric)
library(testthat)

test_that(".ridge_linear_solve errors when X and Y have different rows", {
  X <- matrix(rnorm(6), nrow = 2, ncol = 3)
  Y <- matrix(rnorm(9), nrow = 3, ncol = 3)
  expect_error(fmriparametric:::.ridge_linear_solve(X, Y),
               "same number of rows")
})

test_that("non-finite values in X cause error", {
  X <- matrix(c(1, NA, 3, 4), nrow = 2)
  Y <- matrix(rnorm(4), nrow = 2)
  expect_error(
    fmriparametric:::.ridge_linear_solve(X, Y, 0),
    "X contains non-finite values"
  )
})

test_that("non-finite values in Y cause error", {
  X <- matrix(rnorm(4), nrow = 2)
  Y <- matrix(c(1, Inf, 3, 4), nrow = 2)
  expect_error(
    fmriparametric:::.ridge_linear_solve(X, Y, 0),
    "Y contains non-finite values"
  )
})

test_that("valid call works", {
  X <- matrix(rnorm(20), nrow = 5)
  Y <- matrix(rnorm(20), nrow = 5)
  expect_silent(fmriparametric:::.ridge_linear_solve(X, Y, lambda_ridge = 0))
})

test_that("negative lambda triggers error", {
  X <- matrix(rnorm(20), nrow = 5)
  Y <- matrix(rnorm(20), nrow = 5)
  expect_error(
    fmriparametric:::.ridge_linear_solve(X, Y, lambda_ridge = -1),
    "lambda_ridge"
  )
})

test_that("NA lambda triggers error", {
  X <- matrix(rnorm(20), nrow = 5)
  Y <- matrix(rnorm(20), nrow = 5)
  expect_error(
    fmriparametric:::.ridge_linear_solve(X, Y, lambda_ridge = NA_real_),
    "lambda_ridge"
  )
})