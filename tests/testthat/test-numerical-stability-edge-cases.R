library(fmriparametric)

test_that(".fast_batch_convolution handles edge cases", {
  # Test 1: Empty signal
  kernels <- matrix(rnorm(30 * 3), nrow = 30, ncol = 3)
  expect_error(
    fmriparametric:::.fast_batch_convolution(numeric(0), kernels, 50),
    NA  # Should not error, but return sensible result
  )
  
  # Test 2: Single time point signal
  signal_single <- 1.0
  result <- fmriparametric:::.fast_batch_convolution(signal_single, kernels, 50)
  expect_equal(dim(result), c(50, 3))
  expect_true(all(is.finite(result)))
  
  # Test 3: Zero signal
  signal_zero <- rep(0, 100)
  result <- fmriparametric:::.fast_batch_convolution(signal_zero, kernels, 100)
  expect_equal(dim(result), c(100, 3))
  expect_true(all(result == 0))  # Convolution with zero should be zero
  
  # Test 4: Single kernel column
  kernel_single <- matrix(rnorm(30), ncol = 1)
  signal <- rnorm(100)
  result <- fmriparametric:::.fast_batch_convolution(signal, kernel_single, 100)
  expect_equal(dim(result), c(100, 1))
  expect_true(all(is.finite(result)))
  
  # Test 5: Very long signal (triggers FFT path)
  signal_long <- rnorm(500)
  kernels_long <- matrix(rnorm(50 * 2), nrow = 50, ncol = 2)
  result <- fmriparametric:::.fast_batch_convolution(signal_long, kernels_long, 500)
  expect_equal(dim(result), c(500, 2))
  expect_true(all(is.finite(result)))
  
  # Test 6: Output length shorter than signal
  signal <- rnorm(100)
  kernels <- matrix(rnorm(20 * 2), nrow = 20, ncol = 2)
  result <- fmriparametric:::.fast_batch_convolution(signal, kernels, 50)
  expect_equal(dim(result), c(50, 2))
  
  # Test 7: Extreme values
  signal_extreme <- c(rep(1e10, 50), rep(-1e10, 50))
  kernels_norm <- matrix(rnorm(20 * 2), nrow = 20, ncol = 2)
  result <- fmriparametric:::.fast_batch_convolution(signal_extreme, kernels_norm, 100)
  expect_true(all(is.finite(result)))
})

test_that(".ridge_linear_solve handles edge cases", {
  # Test 1: Square singular matrix
  X_singular <- matrix(c(1, 1, 1, 1), nrow = 2)
  Y <- matrix(c(1, 2), nrow = 2)
  
  # Without ridge penalty, might fail; with ridge, should work
  result <- fmriparametric:::.ridge_linear_solve(X_singular, Y, lambda_ridge = 0.1)
  expect_equal(dim(result), c(2, 1))
  expect_true(all(is.finite(result)))
  
  # Test 2: Tall matrix (overdetermined)
  X_tall <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)
  Y_tall <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  result <- fmriparametric:::.ridge_linear_solve(X_tall, Y_tall, lambda_ridge = 0.01)
  expect_equal(dim(result), c(3, 5))
  expect_true(all(is.finite(result)))
  
  # Test 3: Wide matrix (underdetermined)
  X_wide <- matrix(rnorm(10 * 30), nrow = 10, ncol = 30)
  Y_wide <- matrix(rnorm(10 * 2), nrow = 10, ncol = 2)
  # Need strong regularization for underdetermined system
  result <- fmriparametric:::.ridge_linear_solve(X_wide, Y_wide, lambda_ridge = 1.0)
  expect_equal(dim(result), c(30, 2))
  expect_true(all(is.finite(result)))
  
  # Test 4: Single column X and Y
  X_single <- matrix(rnorm(50), ncol = 1)
  Y_single <- matrix(rnorm(50), ncol = 1)
  result <- fmriparametric:::.ridge_linear_solve(X_single, Y_single, lambda_ridge = 0)
  expect_equal(dim(result), c(1, 1))
  expect_true(is.finite(result[1, 1]))
  
  # Test 5: Near-zero X values
  X_tiny <- matrix(rnorm(50 * 3, sd = 1e-10), nrow = 50, ncol = 3)
  Y_normal <- matrix(rnorm(50 * 2), nrow = 50, ncol = 2)
  result <- fmriparametric:::.ridge_linear_solve(X_tiny, Y_normal, lambda_ridge = 0.1)
  expect_true(all(is.finite(result)))
  
  # Test 6: Identity matrix (should be easy to solve)
  X_identity <- diag(10)
  Y_test <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)
  result <- fmriparametric:::.ridge_linear_solve(X_identity, Y_test, lambda_ridge = 0)
  # Solution should be Y itself
  expect_equal(result, Y_test, tolerance = 1e-10)
})

test_that(".batch_convolution handles degenerate cases", {
  # Test 1: Empty kernels matrix
  signals <- matrix(rnorm(50), ncol = 1)
  expect_error(
    fmriparametric:::.batch_convolution(signals, matrix(nrow = 0, ncol = 0), 50),
    "kernels matrix is empty"
  )
  
  # Test 2: Mismatched dimensions
  signals <- matrix(rnorm(50 * 2), nrow = 50, ncol = 2)
  kernels <- matrix(rnorm(30 * 3), nrow = 30, ncol = 3)
  # This should work - it sums signals first
  result <- fmriparametric:::.batch_convolution(signals, kernels, 60)
  expect_equal(dim(result), c(60, 3))
  
  # Test 3: All-zero signals
  signals_zero <- matrix(0, nrow = 50, ncol = 3)
  kernels <- matrix(rnorm(20 * 2), nrow = 20, ncol = 2)
  result <- fmriparametric:::.batch_convolution(signals_zero, kernels, 50)
  expect_true(all(result == 0))
  
  # Test 4: Single signal column with multiple kernels
  signal_single <- matrix(c(rep(0, 20), 1, rep(0, 29)), ncol = 1)
  kernels_multi <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  result <- fmriparametric:::.batch_convolution(signal_single, kernels_multi, 50)
  expect_equal(dim(result), c(50, 5))
  expect_true(all(is.finite(result)))
})

test_that("HRF interface functions handle boundary conditions", {
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  
  # Test 1: Parameters exactly at lower bounds
  params_lower <- bounds$lower
  t_eval <- seq(0, 30, by = 0.5)
  
  # This should work but with adjusted sigma
  hrf_vals <- fmriparametric:::.lwu_hrf_function(t_eval, params_lower, bounds)
  expect_true(all(is.finite(hrf_vals)))
  expect_true(all(hrf_vals >= 0))  # HRF should be non-negative
  
  # Test 2: Parameters exactly at upper bounds
  params_upper <- bounds$upper
  hrf_vals <- fmriparametric:::.lwu_hrf_function(t_eval, params_upper, bounds)
  expect_true(all(is.finite(hrf_vals)))
  
  # Test 3: Taylor basis at boundaries (with epsilon adjustment)
  params_near_lower <- bounds$lower + c(0.02, 0.02, 0.02)
  basis <- fmriparametric:::.lwu_hrf_taylor_basis_function(params_near_lower, t_eval, bounds)
  expect_equal(dim(basis), c(length(t_eval), 4))
  expect_true(all(is.finite(basis)))
  
  # Test 4: Very short time vector
  t_short <- c(0, 1)
  hrf_short <- fmriparametric:::.lwu_hrf_function(t_short, c(6, 2.5, 0.35), bounds)
  expect_length(hrf_short, 2)
  expect_true(all(is.finite(hrf_short)))
})

test_that("Fit metrics handle edge cases correctly", {
  # Test 1: Perfect fit (R² = 1)
  y_true <- matrix(1:10, ncol = 1)
  y_pred <- y_true  # Perfect prediction
  metrics <- fmriparametric:::.calculate_fit_metrics(
    y_true, y_pred, n_predictors = 2, has_intercept = TRUE
  )
  expect_equal(metrics$r_squared, 1.0)
  expect_equal(metrics$mse, 0.0)
  
  # Test 2: Constant y_true with intercept
  y_const <- matrix(5, nrow = 10, ncol = 2)
  y_pred_good <- y_const
  metrics <- fmriparametric:::.calculate_fit_metrics(
    y_const, y_pred_good, n_predictors = 1, has_intercept = TRUE
  )
  expect_equal(metrics$r_squared, c(1.0, 1.0))
  
  # Test 3: Constant y_true without intercept
  y_const <- matrix(5, nrow = 10, ncol = 1)
  y_pred_zero <- matrix(0, nrow = 10, ncol = 1)
  metrics <- fmriparametric:::.calculate_fit_metrics(
    y_const, y_pred_zero, n_predictors = 1, has_intercept = FALSE
  )
  # Without intercept, can't fit constant well
  expect_true(metrics$r_squared < 0.5)
  
  # Test 4: More predictors than observations
  y_small <- matrix(rnorm(3), ncol = 1)
  y_pred_small <- matrix(rnorm(3), ncol = 1)
  metrics <- fmriparametric:::.calculate_fit_metrics(
    y_small, y_pred_small, n_predictors = 5, has_intercept = TRUE
  )
  expect_true(is.na(metrics$r_squared_adj))  # Adjusted R² undefined
  
  # Test 5: Pre-computed TSS
  y_residual <- matrix(rnorm(20), ncol = 2)
  y_pred <- matrix(rnorm(20), ncol = 2)
  tss_original <- c(100, 200)  # From original data before residualization
  metrics <- fmriparametric:::.calculate_fit_metrics(
    y_residual, y_pred, n_predictors = 3, 
    has_intercept = TRUE, precomputed_tss = tss_original
  )
  expect_equal(metrics$tss, tss_original)
  
  # Test 6: Negative R² clamping
  y_true <- matrix(rnorm(20), ncol = 1)
  y_pred <- matrix(rnorm(20, mean = 10), ncol = 1)  # Very bad prediction
  metrics <- fmriparametric:::.calculate_fit_metrics(
    y_true, y_pred, n_predictors = 1, has_intercept = TRUE
  )
  expect_true(metrics$r_squared >= 0)  # Should be clamped to [0, 1]
  expect_true(metrics$r_squared <= 1)
})