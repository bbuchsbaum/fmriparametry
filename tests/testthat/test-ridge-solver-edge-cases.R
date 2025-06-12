# Edge case tests for ridge regression solver
library(fmriparametric)
library(testthat)

test_that(".ridge_linear_solve handles singular matrices", {
  # Case 1: Rank deficient X
  n <- 20
  p <- 5
  X_rank_def <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  X_rank_def <- cbind(X_rank_def, X_rank_def[, 1:2])  # Duplicate columns
  Y <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  
  # Should not error with ridge regularization
  result <- fmriparametric:::.ridge_linear_solve(X_rank_def, Y, lambda = 0.1)
  
  expect_equal(dim(result), c(5, 2))
  expect_true(all(is.finite(result)))
  
  # Case 2: Near-zero lambda with rank deficient X
  result_small_lambda <- fmriparametric:::.ridge_linear_solve(
    X_rank_def, Y, lambda = 1e-10
  )
  
  expect_true(all(is.finite(result_small_lambda)))
})

test_that(".ridge_linear_solve handles extreme scales", {
  n <- 30
  p <- 4
  
  # Case 1: Very small scale data
  X_small <- matrix(rnorm(n * p, sd = 1e-8), nrow = n, ncol = p)
  Y_small <- matrix(rnorm(n * 2, sd = 1e-8), nrow = n, ncol = 2)
  
  result_small <- fmriparametric:::.ridge_linear_solve(
    X_small, Y_small, lambda = 1e-16
  )
  
  expect_true(all(is.finite(result_small)))
  expect_true(all(abs(result_small) < 1e6))  # Coefficients shouldn't explode
  
  # Case 2: Very large scale data
  X_large <- matrix(rnorm(n * p, sd = 1e6), nrow = n, ncol = p)
  Y_large <- matrix(rnorm(n * 2, sd = 1e6), nrow = n, ncol = 2)
  
  result_large <- fmriparametric:::.ridge_linear_solve(
    X_large, Y_large, lambda = 1e12
  )
  
  expect_true(all(is.finite(result_large)))
})

test_that(".ridge_linear_solve matches R's ridge implementation", {
  set.seed(42)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n), nrow = n, ncol = 1)
  lambda <- 0.5
  
  # Our implementation
  our_result <- fmriparametric:::.ridge_linear_solve(X, Y, lambda)
  
  # Standard ridge formula: (X'X + lambda*I)^{-1} X'Y
  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)
  I_p <- diag(p)
  expected <- solve(XtX + lambda * I_p, XtY)
  
  expect_equal(our_result, expected, tolerance = 1e-10)
})

test_that(".ridge_linear_solve handles zero matrices appropriately", {
  n <- 20
  p <- 5
  
  # Case 1: Zero X matrix
  X_zero <- matrix(0, nrow = n, ncol = p)
  Y <- matrix(rnorm(n), nrow = n, ncol = 1)
  
  result_zero_X <- fmriparametric:::.ridge_linear_solve(
    X_zero, Y, lambda = 0.1
  )
  
  # With zero X, coefficients should be zero
  expect_true(all(result_zero_X == 0))
  
  # Case 2: Zero Y matrix
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y_zero <- matrix(0, nrow = n, ncol = 1)
  
  result_zero_Y <- fmriparametric:::.ridge_linear_solve(
    X, Y_zero, lambda = 0.1
  )
  
  # With zero Y, coefficients should be zero
  expect_true(all(result_zero_Y == 0))
})

test_that(".ridge_linear_solve is numerically stable with correlated predictors", {
  n <- 100
  
  # Create highly correlated predictors
  base <- rnorm(n)
  X <- cbind(
    base + rnorm(n, sd = 0.01),
    base + rnorm(n, sd = 0.01),
    base + rnorm(n, sd = 0.01),
    rnorm(n)  # One independent predictor
  )
  
  # True coefficients
  beta_true <- c(1, 1, 1, 2)
  Y <- X %*% beta_true + rnorm(n, sd = 0.1)
  
  # Test with varying lambda
  lambdas <- c(0, 0.001, 0.01, 0.1, 1, 10)
  
  for (lambda in lambdas) {
    result <- fmriparametric:::.ridge_linear_solve(
      X, matrix(Y, ncol = 1), lambda
    )
    
    expect_true(all(is.finite(result)))
    
    # As lambda increases, coefficients should shrink
    if (lambda > 0) {
      expect_true(sum(result^2) < sum(beta_true^2) * 2)
    }
  }
})

test_that(".ridge_linear_solve handles single column edge cases", {
  n <- 30
  
  # Case 1: Single predictor
  X_single <- matrix(rnorm(n), nrow = n, ncol = 1)
  Y <- matrix(2 * X_single + rnorm(n, sd = 0.1), nrow = n, ncol = 1)
  
  result_single <- fmriparametric:::.ridge_linear_solve(
    X_single, Y, lambda = 0.01
  )
  
  expect_equal(length(result_single), 1)
  expect_true(abs(result_single - 2) < 0.5)  # Should be close to 2
  
  # Case 2: Multiple Y columns with single X
  Y_multi <- cbind(
    2 * X_single + rnorm(n, sd = 0.1),
    -3 * X_single + rnorm(n, sd = 0.1),
    0.5 * X_single + rnorm(n, sd = 0.1)
  )
  
  result_multi_y <- fmriparametric:::.ridge_linear_solve(
    X_single, Y_multi, lambda = 0.01
  )
  
  expect_equal(dim(result_multi_y), c(1, 3))
  expect_true(abs(result_multi_y[1, 1] - 2) < 0.5)
  expect_true(abs(result_multi_y[1, 2] - (-3)) < 0.5)
  expect_true(abs(result_multi_y[1, 3] - 0.5) < 0.5)
})