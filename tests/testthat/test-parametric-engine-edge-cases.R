# Edge case tests for parametric engine
library(fmriparametric)
library(testthat)

test_that(".parametric_engine handles zero stimulus gracefully", {
  n_time <- 20
  n_vox <- 3
  
  Y_proj <- matrix(rnorm(n_time * n_vox), ncol = n_vox)
  S_zero <- matrix(0, nrow = n_time, ncol = 1)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  
  result <- fmriparametric:::.parametric_engine(
    Y_proj = Y_proj,
    S_target_proj = S_zero,
    hrf_eval_times = seq(0, 30, length.out = 31),
    hrf_interface = hrf_interface,
    theta_seed = hrf_interface$default_seed(),
    theta_bounds = hrf_interface$default_bounds(),
    lambda_ridge = 0.01,
    epsilon_beta = 1e-6,
    baseline_model = "intercept",
    verbose = FALSE
  )
  
  # Should return valid results with zero R-squared
  expect_equal(nrow(result$theta_hat), n_vox)
  expect_equal(result$r_squared, rep(0, n_vox))
  expect_true(all(is.finite(result$theta_hat)))
  expect_true(all(result$beta0 == 0))
})

test_that(".parametric_engine avoids divide-by-zero with tiny amplitudes", {
  n_time <- 30
  n_vox <- 4
  
  # Create data that will produce very small amplitudes
  Y_tiny <- matrix(c(
    rep(1e-10, n_time),      # Positive tiny
    rep(-1e-10, n_time),     # Negative tiny
    rep(0, n_time),          # Exact zero
    rnorm(n_time, sd = 0.1)  # Normal scale
  ), ncol = n_vox)
  
  S <- matrix(c(rep(0, 10), rep(1, 5), rep(0, 15)), ncol = 1)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  
  result <- fmriparametric:::.parametric_engine(
    Y_proj = Y_tiny,
    S_target_proj = S,
    hrf_eval_times = seq(0, 30, length.out = 31),
    hrf_interface = hrf_interface,
    theta_seed = hrf_interface$default_seed(),
    theta_bounds = hrf_interface$default_bounds(),
    lambda_ridge = 0.01,
    epsilon_beta = 1e-6,
    baseline_model = "intercept",
    verbose = FALSE
  )
  
  # All results should be finite
  expect_true(all(is.finite(result$theta_hat)))
  expect_true(all(is.finite(result$beta0)))
  expect_true(all(is.finite(result$r_squared)))
  
  # Parameters should stay within bounds
  bounds <- hrf_interface$default_bounds()
  for (i in 1:n_vox) {
    expect_true(all(result$theta_hat[i, ] >= bounds$lower))
    expect_true(all(result$theta_hat[i, ] <= bounds$upper))
  }
})

test_that(".parametric_engine handles constant Y (zero variance)", {
  n_time <- 25
  n_vox <- 2
  
  Y_const <- matrix(5, nrow = n_time, ncol = n_vox)
  S <- matrix(runif(n_time), ncol = 1)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  
  result <- fmriparametric:::.parametric_engine(
    Y_proj = Y_const,
    S_target_proj = S,
    hrf_eval_times = seq(0, 30, length.out = 31),
    hrf_interface = hrf_interface,
    theta_seed = hrf_interface$default_seed(),
    theta_bounds = hrf_interface$default_bounds(),
    lambda_ridge = 0.01,
    epsilon_beta = 1e-6,
    baseline_model = "intercept",
    verbose = FALSE
  )
  
  # R-squared should be 0 for constant Y
  expect_equal(result$r_squared, rep(0, n_vox))
  expect_true(all(is.finite(result$theta_hat)))
})

test_that(".parametric_engine handles extreme parameter seeds", {
  n_time <- 40
  n_vox <- 1
  
  Y <- matrix(rnorm(n_time), ncol = n_vox)
  S <- matrix(c(rep(0, 10), rep(1, 5), rep(0, 25)), ncol = 1)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  bounds <- hrf_interface$default_bounds()
  
  # Test with seeds at bounds
  seeds_to_test <- list(
    lower = bounds$lower,
    upper = bounds$upper,
    mixed = c(bounds$lower[1], bounds$upper[2], bounds$lower[3])
  )
  
  for (seed_name in names(seeds_to_test)) {
    seed <- seeds_to_test[[seed_name]]
    
    result <- fmriparametric:::.parametric_engine(
      Y_proj = Y,
      S_target_proj = S,
      hrf_eval_times = seq(0, 30, length.out = 31),
      hrf_interface = hrf_interface,
      theta_seed = seed,
      theta_bounds = bounds,
      lambda_ridge = 0.01,
      epsilon_beta = 1e-6,
      baseline_model = "intercept",
      verbose = FALSE
    )
    
    expect_true(all(is.finite(result$theta_hat)))
    expect_true(all(result$theta_hat >= bounds$lower))
    expect_true(all(result$theta_hat <= bounds$upper))
  }
})

test_that(".parametric_engine works with different baseline models", {
  n_time <- 30
  n_vox <- 2
  
  Y <- matrix(rnorm(n_time * n_vox), ncol = n_vox)
  S <- matrix(rnorm(n_time), ncol = 1)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  
  # Test with intercept
  result_intercept <- fmriparametric:::.parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    hrf_eval_times = seq(0, 30, length.out = 31),
    hrf_interface = hrf_interface,
    theta_seed = hrf_interface$default_seed(),
    theta_bounds = hrf_interface$default_bounds(),
    lambda_ridge = 0.01,
    epsilon_beta = 1e-6,
    baseline_model = "intercept",
    verbose = FALSE
  )
  
  # Test without intercept
  result_no_intercept <- fmriparametric:::.parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    hrf_eval_times = seq(0, 30, length.out = 31),
    hrf_interface = hrf_interface,
    theta_seed = hrf_interface$default_seed(),
    theta_bounds = hrf_interface$default_bounds(),
    lambda_ridge = 0.01,
    epsilon_beta = 1e-6,
    baseline_model = NULL,
    verbose = FALSE
  )
  
  # Both should produce valid results
  expect_true(all(is.finite(result_intercept$theta_hat)))
  expect_true(all(is.finite(result_no_intercept$theta_hat)))
  
  # Results should differ (intercept affects fit)
  expect_false(all(result_intercept$theta_hat == result_no_intercept$theta_hat))
})

test_that(".parametric_engine handles sparse event designs", {
  n_time <- 100
  n_vox <- 2
  
  # Very sparse stimulus (single brief event)
  S_sparse <- matrix(0, nrow = n_time, ncol = 1)
  S_sparse[50] <- 1
  
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  bounds <- hrf_interface$default_bounds()
  
  # Generate response
  hrf <- hrf_interface$hrf_function(
    seq(0, 30, length.out = 31),
    c(6, 3, 0.5)
  )
  conv_result <- fmriparametric:::.fast_batch_convolution(
    matrix(hrf, ncol = 1),
    S_sparse,
    n_time
  )
  
  Y <- matrix(0, nrow = n_time, ncol = n_vox)
  Y[, 1] <- conv_result[, 1] + rnorm(n_time, sd = 0.1)
  Y[, 2] <- rnorm(n_time)  # Pure noise
  
  result <- fmriparametric:::.parametric_engine(
    Y_proj = Y,
    S_target_proj = S_sparse,
    hrf_eval_times = seq(0, 30, length.out = 31),
    hrf_interface = hrf_interface,
    theta_seed = c(5, 4, 0.3),
    theta_bounds = bounds,
    lambda_ridge = 0.01,
    epsilon_beta = 1e-6,
    baseline_model = "intercept",
    verbose = FALSE
  )
  
  expect_true(all(is.finite(result$theta_hat)))
  expect_true(result$r_squared[1] > result$r_squared[2])  # Signal should fit better than noise
})