# Edge case tests for Gauss-Newton refinement
library(fmriparametric)
library(testthat)

test_that(".calculate_objective_gn handles singular systems", {
  # Case 1: Zero stimulus matrix
  n_time <- 10
  y <- rnorm(n_time)
  S <- matrix(0, nrow = n_time, ncol = 1)
  t_hrf <- seq(0, 30, length.out = 31)
  theta <- c(lag = 6, width = 3, undershoot = 0.5)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  
  obj <- fmriparametric:::.calculate_objective_gn(
    theta = theta,
    y = y,
    S = S,
    t_hrf = t_hrf,
    hrf_interface = hrf_interface,
    n_time = n_time
  )
  
  expect_equal(obj, Inf)
  
  # Case 2: Near-zero predictor magnitude
  S_small <- matrix(1e-15, nrow = n_time, ncol = 1)
  
  obj_small <- fmriparametric:::.calculate_objective_gn(
    theta = theta,
    y = y,
    S = S_small,
    t_hrf = t_hrf,
    hrf_interface = hrf_interface,
    n_time = n_time
  )
  
  expect_equal(obj_small, Inf)
})

test_that(".get_jacobian_and_residuals returns NULL for singular cases", {
  n_time <- 10
  y <- rnorm(n_time)
  S <- matrix(0, nrow = n_time, ncol = 1)
  t_hrf <- seq(0, 30, length.out = 31)
  theta <- c(lag = 6, width = 3, undershoot = 0.5)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  
  result <- fmriparametric:::.get_jacobian_and_residuals(
    theta = theta,
    y = y,
    S = S,
    t_hrf = t_hrf,
    hrf_interface = hrf_interface,
    n_time = n_time
  )
  
  expect_null(result)
})

test_that(".gauss_newton_refinement handles degenerate inputs gracefully", {
  n_time <- 20
  n_vox <- 2
  
  # Case 1: Zero stimulus
  Y_proj <- matrix(rnorm(n_time * n_vox), ncol = n_vox)
  S_zero <- matrix(0, nrow = n_time, ncol = 1)
  theta_start <- matrix(c(6, 3, 0.5), nrow = 1)
  theta_start <- theta_start[rep(1, n_vox), , drop = FALSE]
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  
  result_zero <- fmriparametric:::.gauss_newton_refinement(
    theta_hat_voxel = theta_start,
    r2_voxel = rep(0.1, n_vox),
    Y_proj = Y_proj,
    S_target_proj = S_zero,
    scan_times = seq_len(n_time),
    hrf_eval_times = seq(0, 30, length.out = 31),
    hrf_interface = hrf_interface,
    theta_bounds = hrf_interface$default_bounds(),
    queue_labels = rep("hard_GN", n_vox),
    max_iter_gn = 5,
    verbose = FALSE
  )
  
  expect_true(all(result_zero$convergence_status == "singular_system"))
  expect_equal(result_zero$iteration_counts, rep(0, n_vox))
  
  # Case 2: Constant Y (zero variance)
  Y_const <- matrix(5, nrow = n_time, ncol = n_vox)
  S_normal <- matrix(rnorm(n_time), nrow = n_time, ncol = 1)
  
  result_const <- fmriparametric:::.gauss_newton_refinement(
    theta_hat_voxel = theta_start,
    r2_voxel = rep(0, n_vox),
    Y_proj = Y_const,
    S_target_proj = S_normal,
    scan_times = seq_len(n_time),
    hrf_eval_times = seq(0, 30, length.out = 31),
    hrf_interface = hrf_interface,
    theta_bounds = hrf_interface$default_bounds(),
    queue_labels = rep("moderate", n_vox),
    max_iter_gn = 5,
    verbose = FALSE
  )
  
  expect_equal(result_const$r2, rep(0, n_vox))
  expect_true(all(result_const$convergence_status %in% c("converged", "max_iter")))
})

test_that(".gauss_newton_refinement respects parameter bounds", {
  n_time <- 50
  n_vox <- 3
  
  # Generate data with true parameters near bounds
  theta_true <- matrix(c(
    2, 0.1, 0.01,    # Near lower bounds
    15, 10, 1.49,    # Near upper bounds  
    8, 5, 0.8        # Middle values
  ), ncol = 3, byrow = TRUE)
  
  Y_proj <- matrix(0, nrow = n_time, ncol = n_vox)
  S <- matrix(c(rep(0, 10), rep(1, 5), rep(0, 35)), ncol = 1)
  t_hrf <- seq(0, 30, length.out = 31)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  bounds <- hrf_interface$default_bounds()
  
  # Generate Y from true parameters
  for (i in 1:n_vox) {
    hrf <- hrf_interface$hrf_function(
      t_hrf,
      theta_true[i, ]
    )
    conv_result <- fmriparametric:::.fast_batch_convolution(
      matrix(hrf, ncol = 1),
      S,
      n_time
    )
    Y_proj[, i] <- conv_result[, 1] + rnorm(n_time, sd = 0.1)
  }
  
  # Start far from true values
  theta_start <- matrix(c(8, 4, 0.7), nrow = 1)
  theta_start <- theta_start[rep(1, n_vox), , drop = FALSE]
  
  result <- fmriparametric:::.gauss_newton_refinement(
    theta_hat_voxel = theta_start,
    r2_voxel = rep(0.5, n_vox),
    Y_proj = Y_proj,
    S_target_proj = S,
    scan_times = seq_len(n_time),
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_interface,
    theta_bounds = bounds,
    queue_labels = rep("hard_GN", n_vox),
    max_iter_gn = 20,
    verbose = FALSE
  )
  
  # Check all parameters are within bounds
  for (i in 1:n_vox) {
    expect_true(all(result$theta_hat[i, ] >= bounds$lower))
    expect_true(all(result$theta_hat[i, ] <= bounds$upper))
  }
  
  # Should improve fit
  expect_true(all(result$r2 >= 0.5))
})

test_that(".gauss_newton_refinement handles numerical precision issues", {
  n_time <- 30
  n_vox <- 2
  
  # Case 1: Very small scale data
  Y_small <- matrix(rnorm(n_time * n_vox, sd = 1e-8), ncol = n_vox)
  S <- matrix(c(rep(0, 5), rep(1, 3), rep(0, 22)), ncol = 1)
  theta_start <- matrix(c(6, 3, 0.5), nrow = 1)
  theta_start <- theta_start[rep(1, n_vox), , drop = FALSE]
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  
  result_small <- fmriparametric:::.gauss_newton_refinement(
    theta_hat_voxel = theta_start,
    r2_voxel = rep(0.1, n_vox),
    Y_proj = Y_small,
    S_target_proj = S,
    scan_times = seq_len(n_time),
    hrf_eval_times = seq(0, 30, length.out = 31),
    hrf_interface = hrf_interface,
    theta_bounds = hrf_interface$default_bounds(),
    queue_labels = rep("moderate", n_vox),
    max_iter_gn = 5,
    lambda_ridge = 0.1,  # Higher ridge for stability
    verbose = FALSE
  )
  
  expect_true(all(is.finite(result_small$theta_hat)))
  expect_true(all(is.finite(result_small$r2)))
  
  # Case 2: Very large scale data
  Y_large <- matrix(rnorm(n_time * n_vox, sd = 1e6), ncol = n_vox)
  
  result_large <- fmriparametric:::.gauss_newton_refinement(
    theta_hat_voxel = theta_start,
    r2_voxel = rep(0.1, n_vox),
    Y_proj = Y_large,
    S_target_proj = S,
    scan_times = seq_len(n_time),
    hrf_eval_times = seq(0, 30, length.out = 31),
    hrf_interface = hrf_interface,
    theta_bounds = hrf_interface$default_bounds(),
    queue_labels = rep("hard_GN", n_vox),
    max_iter_gn = 5,
    verbose = FALSE
  )
  
  expect_true(all(is.finite(result_large$theta_hat)))
  expect_true(all(is.finite(result_large$r2)))
})