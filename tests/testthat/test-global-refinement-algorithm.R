library(testthat)
library(fmriparametric)

# Test that global iterative refinement improves fit for a simple dataset

test_that("global refinement improves parameter estimates", {
  set.seed(42)
  n_time <- 60
  onsets <- rep(0, n_time)
  onsets[c(10, 30, 50)] <- 1

  true_theta <- c(tau = 6, sigma = 2.5, rho = 0.4)
  t_hrf <- seq(0, 30, by = 0.5)
  
  # Get proper bounds for HRF function
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  true_hrf <- fmriparametric:::.lwu_hrf_function(t_hrf, true_theta, bounds)

  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  signal <- conv_full[seq_len(n_time)]
  Y <- matrix(signal + rnorm(n_time, sd = 0.05), ncol = 1)

  event_mat <- matrix(onsets, ncol = 1)

  bad_seed <- c(9, 4, 0.1)

  fit_no_refine <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = event_mat,
    parametric_hrf = "lwu",
    theta_seed = bad_seed,
    global_refinement = FALSE,
    verbose = FALSE
  )

  fit_refine <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = event_mat,
    parametric_hrf = "lwu",
    theta_seed = bad_seed,
    global_refinement = TRUE,
    global_passes = 2,
    verbose = FALSE
  )

  # Refined version should have better R-squared
  expect_gt(mean(fit_refine$r_squared), mean(fit_no_refine$r_squared))

  # Refined version should be closer to true parameters
  dist_no_ref <- sum((coef(fit_no_refine)[1, ] - true_theta)^2)
  dist_refine <- sum((coef(fit_refine)[1, ] - true_theta)^2)
  expect_lt(dist_refine, dist_no_ref)
})

test_that("global refinement works with multiple voxels", {
  set.seed(123)
  n_time <- 80
  n_vox <- 5
  
  # Create event design
  onsets <- rep(0, n_time)
  onsets[c(15, 35, 55, 75)] <- 1
  event_mat <- matrix(onsets, ncol = 1)
  
  # Generate data with different true parameters per voxel
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  t_hrf <- seq(0, 30, by = 0.5)
  Y <- matrix(NA, n_time, n_vox)
  
  true_params <- matrix(c(
    5.5, 2.0, 0.3,
    6.5, 2.5, 0.4,
    7.0, 3.0, 0.2,
    5.0, 2.2, 0.5,
    6.0, 2.8, 0.35
  ), nrow = n_vox, byrow = TRUE)
  
  for (v in 1:n_vox) {
    true_hrf <- fmriparametric:::.lwu_hrf_function(t_hrf, true_params[v, ], bounds)
    conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
    signal <- conv_full[seq_len(n_time)]
    Y[, v] <- signal + rnorm(n_time, sd = 0.08)
  }
  
  # Bad starting seed for all voxels
  bad_seed <- c(8, 4, 0.1)
  
  fit_no_refine <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = event_mat,
    parametric_hrf = "lwu",
    theta_seed = bad_seed,
    global_refinement = FALSE,
    verbose = FALSE
  )
  
  fit_refine <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = event_mat,
    parametric_hrf = "lwu",
    theta_seed = bad_seed,
    global_refinement = TRUE,
    global_passes = 3,
    verbose = FALSE
  )
  
  # Expect improvement in average R-squared
  expect_gt(mean(fit_refine$r_squared), mean(fit_no_refine$r_squared))
  
  # Expect improvement in parameter accuracy for most voxels
  errors_no_refine <- rowSums((coef(fit_no_refine) - true_params)^2)
  errors_refine <- rowSums((coef(fit_refine) - true_params)^2)
  
  # Most voxels should improve
  improved_count <- sum(errors_refine < errors_no_refine)
  expect_gt(improved_count, n_vox / 2)
})

test_that("global refinement convergence works properly", {
  set.seed(456)
  n_time <- 100
  
  # Create clear signal
  onsets <- rep(0, n_time)
  onsets[seq(20, 80, by = 20)] <- 1
  event_mat <- matrix(onsets, ncol = 1)
  
  # Generate high-quality synthetic data
  true_theta <- c(6, 2.5, 0.35)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  t_hrf <- seq(0, 30, by = 0.5)
  true_hrf <- fmriparametric:::.lwu_hrf_function(t_hrf, true_theta, bounds)
  
  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  signal <- conv_full[seq_len(n_time)]
  Y <- matrix(signal + rnorm(n_time, sd = 0.02), ncol = 1)  # Low noise
  
  # Test different numbers of global passes
  results <- list()
  for (passes in c(1, 2, 5)) {
    results[[paste0("passes_", passes)]] <- estimate_parametric_hrf(
      fmri_data = Y,
      event_model = event_mat,
      parametric_hrf = "lwu",
      theta_seed = c(8, 4, 0.1),  # Bad seed
      global_refinement = TRUE,
      global_passes = passes,
      verbose = FALSE
    )
  }
  
  # More passes should generally lead to better or equal performance
  r2_1 <- results$passes_1$r_squared[1]
  r2_2 <- results$passes_2$r_squared[1]
  r2_5 <- results$passes_5$r_squared[1]
  
  expect_gte(r2_2, r2_1 - 0.01)  # Allow small tolerance for numerical differences
  expect_gte(r2_5, r2_2 - 0.01)
  
  # All should achieve decent fit with good data
  expect_gt(r2_1, 0.5)
  expect_gt(r2_2, 0.5)
  expect_gt(r2_5, 0.5)
})

test_that("global refinement handles difficult optimization cases", {
  set.seed(789)
  n_time <- 60
  n_vox <- 3
  
  # Create challenging scenario with weak signal
  onsets <- rep(0, n_time)
  onsets[c(15, 35, 55)] <- 1
  event_mat <- matrix(onsets, ncol = 1)
  
  # Generate weak signal with more noise
  true_theta <- c(5.5, 2.0, 0.25)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  t_hrf <- seq(0, 30, by = 0.5)
  true_hrf <- fmriparametric:::.lwu_hrf_function(t_hrf, true_theta, bounds)
  
  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  signal <- conv_full[seq_len(n_time)]
  
  Y <- matrix(NA, n_time, n_vox)
  for (v in 1:n_vox) {
    # Add substantial noise
    Y[, v] <- 0.5 * signal + rnorm(n_time, sd = 0.3)
  }
  
  # Very bad starting seed
  bad_seed <- c(12, 5, 0.05)
  
  # Should not crash and should produce finite results
  expect_no_error({
    fit <- estimate_parametric_hrf(
      fmri_data = Y,
      event_model = event_mat,
      parametric_hrf = "lwu",
      theta_seed = bad_seed,
      global_refinement = TRUE,
      global_passes = 2,
      verbose = FALSE
    )
  })
  
  # Results should be finite and within bounds
  params <- coef(fit)
  expect_true(all(is.finite(params)))
  expect_true(all(params[, 1] >= bounds$lower[1] & params[, 1] <= bounds$upper[1]))
  expect_true(all(params[, 2] >= bounds$lower[2] & params[, 2] <= bounds$upper[2]))
  expect_true(all(params[, 3] >= bounds$lower[3] & params[, 3] <= bounds$upper[3]))
})

test_that("global refinement shows improvement over baseline", {
  set.seed(999)
  n_time <- 60
  n_vox <- 8
  
  # Create realistic test scenario
  onsets <- rep(0, n_time)
  onsets[c(12, 25, 38, 51)] <- 1  # Well-spaced events
  event_mat <- matrix(onsets, ncol = 1)
  
  # Generate data with known good HRF
  true_theta <- c(6.2, 2.3, 0.38)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  t_hrf <- seq(0, 30, by = 0.5)
  true_hrf <- fmriparametric:::.lwu_hrf_function(t_hrf, true_theta, bounds)
  
  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  signal <- conv_full[seq_len(n_time)]
  
  Y <- matrix(NA, n_time, n_vox)
  for (v in 1:n_vox) {
    # Add realistic noise level
    amplitude <- 1 + v * 0.1  # Varying amplitudes
    Y[, v] <- amplitude * signal + rnorm(n_time, sd = 0.15)
  }
  
  # Test with moderately bad seed
  moderately_bad_seed <- c(8.5, 3.5, 0.15)
  
  fit_global <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = event_mat,
    parametric_hrf = "lwu",
    theta_seed = moderately_bad_seed,
    global_refinement = TRUE,
    global_passes = 3,
    verbose = FALSE
  )
  
  fit_no_global <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = event_mat,
    parametric_hrf = "lwu",
    theta_seed = moderately_bad_seed,
    global_refinement = FALSE,
    verbose = FALSE
  )
  
  # Global refinement should generally improve fit quality
  expect_gte(mean(fit_global$r_squared), mean(fit_no_global$r_squared) - 0.02)
  
  # Both should produce valid results
  expect_true(all(is.finite(coef(fit_global))))
  expect_true(all(is.finite(coef(fit_no_global))))
  expect_true(all(fit_global$r_squared >= 0))
  expect_true(all(fit_no_global$r_squared >= 0))
})

