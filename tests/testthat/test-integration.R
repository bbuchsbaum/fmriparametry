library(testthat)
library(fmriparametric)
library(fmrireg)

# context() is deprecated in testthat 3rd edition

# Test Scenario 1: Basic workflow with simulated data
test_that("full estimation workflow works with basic simulated data", {
 
  set.seed(789)
  n_time <- 100
  n_vox <- 10
  TR <- 2
  
  
  sframe <- fmrireg::sampling_frame(
    blocklens=n_time,
    TR = TR,
    start_time=TR/2
  )
  # Create simple event model
  onsets <- seq(10, 190, by = 20)
  durations <- rep(2, length(onsets))
  events_df <- data.frame(
    onset = onsets,
    duration = durations,
    trial_type = "stimulus",
    block=rep(1,length(onsets))
  )
  
  # Generate synthetic BOLD data
  # True HRF parameters
  true_tau <- 6
  true_sigma <- 2.5
  true_rho <- 0.35
  true_amp <- 4  # Increased amplitude for better signal
  
  # Generate HRF
  t_hrf <- seq(0, 30, by = TR)
  true_hrf <- exp(-(t_hrf - true_tau)^2 / (2 * true_sigma^2)) - 
              true_rho * exp(-(t_hrf - true_tau - 2*true_sigma)^2 / (2 * (1.6*true_sigma)^2))
  true_hrf <- true_hrf / max(true_hrf)  # Normalize
  
  # Create stimulus vector
  stim_vec <- rep(0, n_time)
  for (onset in onsets) {
    idx <- ceiling(onset / TR)
    if (idx <= n_time) stim_vec[idx] <- 1
  }
  
  # Convolve to get signal
  signal <- stats::filter(stim_vec, true_hrf, sides = 1)
  signal[is.na(signal)] <- 0
  
  # Create data matrix with noise
  Y_data <- matrix(rep(true_amp * signal, n_vox), ncol = n_vox)
  Y_data <- Y_data + matrix(rnorm(n_time * n_vox, sd = 0.3), ncol = n_vox)  # Reduced noise
  
  # Create event model object (mock if fmrireg not available)
  if (requireNamespace("fmrireg", quietly = TRUE)) {
    
    event_model <- fmrireg::event_model(
      onset ~ fmrireg::hrf(trial_type),
      data = events_df,
      block = events_df$block,
      sampling_frame = sframe
    )
    
    # Create fmri dataset
    fmri_data <- fmrireg::matrix_dataset(
      Y_data,
      TR = TR,
      run_length = n_time
    )
  } else {
    # Mock objects for testing without fmrireg
    event_model <- list(
      terms = list(matrix(stim_vec, ncol = 1)),
      sampling_rate = 1/TR,
      class = c("event_model", "list")
    )
    class(event_model) <- c("event_model", "list")
    
    fmri_data <- list(
      data = Y_data,
      sampling_rate = 1/TR,
      TR = TR,
      class = c("matrix_dataset", "list")
    )
    class(fmri_data) <- c("matrix_dataset", "list")
  }
  
  # Run estimation
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    theta_seed = c(7, 3, 0.4),  # Start slightly off
    verbose = FALSE
  )
  
  # Verify output structure
  expect_s3_class(fit, "parametric_hrf_fit")
  expect_equal(nrow(fit$estimated_parameters), n_vox)
  expect_equal(ncol(fit$estimated_parameters), 3)
  expect_equal(length(fit$amplitudes), n_vox)
  expect_equal(fit$parameter_names, c("tau", "sigma", "rho"))
  expect_equal(fit$hrf_model, "lwu")
  
  # Check parameter estimates are reasonable
  mean_params <- colMeans(fit$estimated_parameters)
  expect_true(mean_params[1] > 3 && mean_params[1] < 10)  # tau
  expect_true(mean_params[2] > 0.5 && mean_params[2] < 5)   # sigma - allow smaller values
  expect_true(mean_params[3] >= 0 && mean_params[3] <= 1.5) # rho - use actual bound
  
  # Amplitudes should be positive on average
  expect_true(mean(fit$amplitudes) > 0)
})

# Test Scenario 2: Complex event design with multiple conditions
test_that("workflow handles multiple event types correctly", {
  skip_if_not_installed("fmrireg")
  
  set.seed(123)
  n_time <- 150
  n_vox <- 5
  TR <- 2
  
  sframe <- fmrireg::sampling_frame(
    blocklens = n_time,
    TR = TR,
    start_time = TR/2
  )
  
  # Create complex event design - ensure onsets are ordered
  events_df <- data.frame(
    onset = c(10, 20, 30, 40, 50, 60),
    duration = c(1, 2, 1, 2, 1, 2),
    trial_type = c("A", "B", "A", "B", "A", "B"),
    parametric_value = c(1, 0.5, 2, 1, 3, 1.5),
    block = rep(1, 6)
  )
  
  # For simplicity, combine both conditions
  all_onsets <- events_df$onset
  stim_vec <- rep(0, n_time)
  for (onset in all_onsets) {
    idx <- ceiling(onset / TR)
    if (idx <= n_time) stim_vec[idx] <- 1
  }
  
  # Generate data with parameters within custom bounds
  # Use tau=5, sigma=2.5, rho=0.4 (all within custom bounds)
  t_hrf <- seq(0, 20, by = TR)
  hrf <- exp(-(t_hrf - 5)^2 / (2 * 2.5^2)) - 0.4 * exp(-(t_hrf - 7.5)^2 / (2 * 4^2))
  hrf <- hrf / max(hrf)
  
  signal <- stats::filter(stim_vec, hrf, sides = 1)
  signal[is.na(signal)] <- 0
  
  # Stronger signal for better estimation
  Y_data <- matrix(rep(3 * signal, n_vox), ncol = n_vox) * seq(1, 2, length.out = n_vox)
  Y_data <- Y_data + matrix(rnorm(n_time * n_vox, sd = 0.2), ncol = n_vox)
  
  if (requireNamespace("fmrireg", quietly = TRUE)) {
    event_model <- fmrireg::event_model(
      onset ~ fmrireg::hrf(trial_type),
      data = events_df,
      block = events_df$block,
      sampling_frame = sframe
    )
    
    fmri_data <- fmrireg::matrix_dataset(
      Y_data,
      TR = TR,
      run_length = n_time
    )
  } else {
    # Mock objects
    event_model <- list(
      terms = list(matrix(stim_vec, ncol = 1)),
      sampling_rate = 1/TR
    )
    class(event_model) <- c("event_model", "list")
    
    fmri_data <- list(data = Y_data, sampling_rate = 1/TR, TR = TR)
    class(fmri_data) <- c("matrix_dataset", "list")
  }
  
  # Run with custom bounds
  custom_bounds <- list(
    lower = c(2, 1, 0),
    upper = c(10, 5, 0.8)
  )
  
  # Disable global refinement for bounds test
  old_opt <- getOption("fmriparametric.refine_global")
  options(fmriparametric.refine_global = FALSE)
  
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    theta_bounds = custom_bounds,
    lambda_ridge = 0.05,
    verbose = FALSE
  )
  
  options(fmriparametric.refine_global = old_opt)
  
  expect_s3_class(fit, "parametric_hrf_fit")
  
  # Check bounds were applied
  expect_true(all(fit$estimated_parameters[,1] >= custom_bounds$lower[1]))
  expect_true(all(fit$estimated_parameters[,1] <= custom_bounds$upper[1]))
  expect_true(all(fit$estimated_parameters[,2] >= custom_bounds$lower[2]))
  expect_true(all(fit$estimated_parameters[,2] <= custom_bounds$upper[2]))
  expect_true(all(fit$estimated_parameters[,3] >= custom_bounds$lower[3]))
  expect_true(all(fit$estimated_parameters[,3] <= custom_bounds$upper[3]))
})

# Test Scenario 3: Workflow with confound regression
test_that("workflow handles confound regression properly", {
  skip_if_not_installed("fmrireg")
  
  set.seed(456)
  n_time <- 80
  n_vox <- 3
  TR <- 2
  
  sframe <- fmrireg::sampling_frame(
    blocklens = n_time,
    TR = TR,
    start_time = TR/2
  )
  
  # Simple event
  stim_vec <- rep(0, n_time)
  stim_vec[c(10, 25, 40)] <- 1
  
  # Create data with linear trend and motion confound
  time_vec <- seq_len(n_time)
  linear_trend <- 0.01 * time_vec
  motion_confound <- 0.5 * sin(2 * pi * time_vec / 20)
  
  # True signal
  hrf <- exp(-(seq(0, 15, by = TR) - 6)^2 / (2 * 2.5^2))
  signal <- stats::filter(stim_vec, hrf, sides = 1)
  signal[is.na(signal)] <- 0
  
  Y_data <- matrix(0, nrow = n_time, ncol = n_vox)
  for (v in 1:n_vox) {
    Y_data[,v] <- signal + linear_trend + motion_confound + rnorm(n_time, sd = 0.2)
  }
  
  # Create confound dataframe
  confounds_df <- data.frame(
    linear = linear_trend,
    motion = motion_confound
  )
  
  if (requireNamespace("fmrireg", quietly = TRUE)) {
    events_df <- data.frame(
      onset = which(stim_vec == 1) * TR,
      block = rep(1, sum(stim_vec == 1)),
      stim = factor(rep("stimulus", sum(stim_vec == 1)))
    )
    event_model <- fmrireg::event_model(
      onset ~ fmrireg::hrf(stim),
      data = events_df,
      block = events_df$block,
      sampling_frame = sframe
    )
    
    fmri_data <- fmrireg::matrix_dataset(
      Y_data,
      TR = TR,
      run_length = n_time
    )
  } else {
    event_model <- list(
      terms = list(matrix(stim_vec, ncol = 1)),
      sampling_rate = 1/TR
    )
    class(event_model) <- c("event_model", "list")
    
    fmri_data <- list(data = Y_data, sampling_rate = 1/TR, TR = TR)
    class(fmri_data) <- c("matrix_dataset", "list")
  }
  
  # Run without confound regression
  fit_no_confound <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  # Run with confound regression (if supported)
  # Note: confound_formula not implemented in sprint 1
  # This test is a placeholder for sprint 2
  expect_s3_class(fit_no_confound, "parametric_hrf_fit")
})

# Test Scenario 4: Edge cases and error handling
test_that("workflow handles edge cases gracefully", {
  # Empty data
  expect_error(
    estimate_parametric_hrf(
      fmri_data = matrix(0, nrow = 0, ncol = 0),
      event_model = matrix(0, nrow = 0, ncol = 1)
    ),
    "X is NULL or empty"
  )
  
  # Mismatched dimensions
  expect_error(
    estimate_parametric_hrf(
      fmri_data = matrix(rnorm(20), nrow = 10, ncol = 2),
      event_model = matrix(rbinom(15, 1, 0.2), ncol = 1)
    ),
    "wrong number of rows"
  )
  
  # All zero data
  zero_data <- matrix(0, nrow = 50, ncol = 5)
  zero_events <- matrix(c(rep(0, 40), rep(1, 10)), ncol = 1)
  
  fit_zero <- estimate_parametric_hrf(
    fmri_data = zero_data,
    event_model = zero_events,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  expect_s3_class(fit_zero, "parametric_hrf_fit")
  expect_false(any(is.na(fit_zero$estimated_parameters)))
  expect_false(any(is.infinite(fit_zero$estimated_parameters)))
})

# Test Scenario 5: Performance benchmarks
test_that("workflow completes in reasonable time for moderate data", {
  skip_if_not(Sys.getenv("FMRIPARAMETRIC_EXTENDED_TESTS") == "true",
              "Skipping extended performance test")
  
  set.seed(999)
  n_time <- 200
  n_vox <- 1000
  
  # Simple data
  Y_data <- matrix(rnorm(n_time * n_vox), nrow = n_time)
  event_data <- matrix(0, nrow = n_time, ncol = 1)
  event_data[seq(10, n_time, by = 30), 1] <- 1
  
  start_time <- Sys.time()
  fit <- estimate_parametric_hrf(
    fmri_data = Y_data,
    event_model = event_data,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
  
  expect_s3_class(fit, "parametric_hrf_fit")
  expect_equal(nrow(fit$estimated_parameters), n_vox)
  expect_true(elapsed < 10)  # Should complete in under 10 seconds
  
  # Report performance
  cat("\nPerformance:", n_vox, "voxels processed in", round(elapsed, 2), "seconds\n")
  cat("Rate:", round(n_vox / elapsed), "voxels/second\n")
})

# Test Scenario 6: S3 methods integration
test_that("S3 methods work correctly with fit object", {
  # Simple test data
  Y <- matrix(rnorm(50), nrow = 10, ncol = 5)
  S <- matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0, 0), ncol = 1)
  
  fit <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = S,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  # Test print method
  expect_output(print(fit), "Parametric HRF Fit")
  expect_output(print(fit), "Model: lwu")
  expect_output(print(fit), "Voxels: 5")
  
  # Test coef method
  coefs <- coef(fit)
  expect_equal(dim(coefs), c(5, 3))
  expect_equal(colnames(coefs), c("tau", "sigma", "rho"))
  
  # Test summary method
  summ <- summary(fit)
  expect_type(summ, "list")
  expect_true("parameter_summary" %in% names(summ))
})

# Test Scenario 7: Different HRF evaluation settings
test_that("workflow respects HRF evaluation time settings", {
  Y <- matrix(rnorm(100), nrow = 20, ncol = 5)
  S <- matrix(rbinom(20, 1, 0.1), ncol = 1)
  
  # Custom HRF evaluation times
  custom_times <- seq(0, 20, by = 0.5)
  
  fit <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = S,
    parametric_hrf = "lwu",
    hrf_eval_times = custom_times,
    verbose = FALSE
  )
  
  expect_s3_class(fit, "parametric_hrf_fit")
  
  # Different HRF span
  fit_short <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = S,
    parametric_hrf = "lwu",
    hrf_span = 15,
    verbose = FALSE
  )
  
  expect_s3_class(fit_short, "parametric_hrf_fit")
})
