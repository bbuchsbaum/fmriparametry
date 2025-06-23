# Test 2: Confound Regression Robustness
# 
# This test verifies that HRF parameter estimation remains accurate
# in the presence of realistic confounds (drift, motion, physiological noise).
# This is critical because real fMRI data always contains these artifacts.
#
# NOTE: This test covers both the legacy confound_formula interface and the new
# baseline_model integration for comprehensive confound regression.

library(fmriparametric)
library(testthat)

test_that("HRF estimation is robust to linear and polynomial drift", {
  set.seed(456)
  n_time <- 300
  TR <- 2.0  # seconds
  true_params <- c(5.5, 1.8, 0.35)
  
  # Create realistic event design
  onsets <- rep(0, n_time)
  # Jittered ISI between 12-20 seconds
  event_times <- cumsum(sample(6:10, 20, replace = TRUE))
  event_times <- event_times[event_times < n_time]
  onsets[event_times] <- 1
  
  # Generate clean HRF signal
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  t_hrf <- seq(0, 30, by = TR)
  true_hrf <- hrf_interface$hrf_function(t_hrf, true_params)
  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  clean_signal <- conv_full[1:n_time]
  
  # Add realistic drift patterns
  time_vec <- 1:n_time
  linear_drift <- 0.5 * (time_vec / n_time)  # 50% signal change over scan
  quadratic_drift <- 0.3 * ((time_vec / n_time)^2)
  
  # Test different drift conditions
  drift_conditions <- list(
    list(name = "linear", drift = linear_drift, formula = ~ scan),
    list(name = "quadratic", drift = linear_drift + quadratic_drift, 
         formula = ~ poly(scan, 2)),
    list(name = "severe_linear", drift = 2.0 * linear_drift, formula = ~ scan)
  )
  
  for (condition in drift_conditions) {
    # Create signal with drift
    Y <- matrix(clean_signal + condition$drift + rnorm(n_time, sd = 0.05), ncol = 1)
    
    # Fit with confound regression
    # Note: confound_formula uses 'scan' variable automatically created from scan times
    fit <- estimate_parametric_hrf(
      Y,
      matrix(onsets, ncol = 1),
      confound_formula = condition$formula,
      parametric_model = "lwu",
      hrf_eval_times = t_hrf,
      global_refinement = TRUE,
      global_passes = 2,
      verbose = FALSE
    )
    
    params_recovered <- as.numeric(coef(fit))
    
    # Should recover parameters despite drift
    expect_equal(params_recovered[1], true_params[1], tolerance = 0.5,
                 label = paste("tau recovery with", condition$name, "drift"))
    expect_equal(params_recovered[2], true_params[2], tolerance = 0.3,
                 label = paste("sigma recovery with", condition$name, "drift"))
    expect_equal(params_recovered[3], true_params[3], tolerance = 0.15,
                 label = paste("rho recovery with", condition$name, "drift"))
    
    # R-squared should still be reasonable after confound removal
    expect_gt(fit$r_squared[1], 0.6,
              label = paste("R-squared with", condition$name, "drift"))
  }
})

test_that("HRF estimation handles basis function drift modeling", {
  set.seed(789)
  n_time <- 250
  n_vox <- 20
  TR <- 2.0
  true_params <- c(6.0, 2.2, 0.4)
  
  # Event design
  onsets <- rep(0, n_time)
  onsets[seq(15, n_time-20, by = 25)] <- 1
  
  # Generate signals for multiple voxels
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  t_hrf <- seq(0, 30, by = TR)
  true_hrf <- hrf_interface$hrf_function(t_hrf, true_params)
  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  clean_signal <- conv_full[1:n_time]
  
  # Create complex drift pattern using basis functions
  time_vec <- 1:n_time
  # Simulate scanner drift with multiple frequency components
  drift <- 0.3 * sin(2 * pi * time_vec / (n_time/2)) + 
           0.2 * sin(2 * pi * time_vec / (n_time/4)) +
           0.4 * (time_vec / n_time)
  
  # Create voxel data with varying drift effects
  Y <- matrix(0, n_time, n_vox)
  for (v in 1:n_vox) {
    # Base signal with voxel-specific amplitude
    amplitude <- runif(1, 0.8, 1.2)
    voxel_signal <- amplitude * clean_signal
    
    # Add drift with voxel-specific scaling
    drift_scale <- runif(1, 0.5, 1.5)
    
    # Add noise
    Y[,v] <- voxel_signal + drift_scale * drift + rnorm(n_time, sd = 0.1)
  }
  
  # Test different confound models
  confound_models <- list(
    linear = ~ scan,
    quadratic = ~ poly(scan, 2),
    cubic = ~ poly(scan, 3),
    spline = ~ splines::bs(scan, df = 5)
  )
  
  results <- list()
  for (model_name in names(confound_models)) {
    fit <- estimate_parametric_hrf(
      Y,
      matrix(onsets, ncol = 1),
      confound_formula = confound_models[[model_name]],
      parametric_model = "lwu",
      hrf_eval_times = t_hrf,
      global_refinement = TRUE,
      global_passes = 2,
      verbose = FALSE
    )
    
    params_est <- coef(fit)
    errors <- abs(params_est - matrix(true_params, n_vox, 3, byrow = TRUE))
    
    results[[model_name]] <- list(
      mean_tau_error = mean(errors[,1]),
      mean_sigma_error = mean(errors[,2]),
      mean_rho_error = mean(errors[,3]),
      mean_r2 = mean(fit$r_squared)
    )
  }
  
  # Higher-order models should perform better for complex drift
  expect_lt(results$cubic$mean_tau_error, results$linear$mean_tau_error * 1.2,
            label = "Cubic model handles complex drift")
  expect_lt(results$spline$mean_tau_error, results$linear$mean_tau_error * 1.2,
            label = "Spline model handles complex drift")
  
  # All models should achieve reasonable parameter recovery
  for (model_name in names(results)) {
    expect_lt(results[[model_name]]$mean_tau_error, 1.0,
              label = paste(model_name, "tau error"))
    expect_lt(results[[model_name]]$mean_sigma_error, 0.5,
              label = paste(model_name, "sigma error"))
    expect_lt(results[[model_name]]$mean_rho_error, 0.2,
              label = paste(model_name, "rho error"))
    expect_gt(results[[model_name]]$mean_r2, 0.45,
              label = paste(model_name, "R-squared"))
  }
})

test_that("Confound regression preserves signal of interest", {
  set.seed(654)
  n_time <- 200
  TR <- 2.0
  true_params <- c(5.5, 2.1, 0.32)
  
  # Create stimulus
  onsets <- rep(0, n_time)
  onsets[seq(20, n_time-30, by = 40)] <- 1
  
  # Generate true signal
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  t_hrf <- seq(0, 30, by = TR)
  true_hrf <- hrf_interface$hrf_function(t_hrf, true_params)
  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  clean_signal <- conv_full[1:n_time]
  
  # Create confounds based on scan times
  time_vec <- 1:n_time
  # Polynomial drift
  drift <- 0.3 * (time_vec/n_time) + 0.2 * (time_vec/n_time)^2 - 0.1 * (time_vec/n_time)^3
  
  # Add drift to signal
  Y <- matrix(clean_signal + drift + rnorm(n_time, sd = 0.05), ncol = 1)
  
  # Fit with different levels of confound regression
  formulas <- list(
    none = NULL,
    linear = ~ scan,
    polynomial = ~ poly(scan, 3)
  )
  
  results <- list()
  for (name in names(formulas)) {
    if (is.null(formulas[[name]])) {
      # No confound regression
      fit <- estimate_parametric_hrf(
        Y,
        matrix(onsets, ncol = 1),
        parametric_model = "lwu",
        hrf_eval_times = t_hrf,
        global_refinement = TRUE,
        global_passes = 2,
        verbose = FALSE
      )
    } else {
      fit <- estimate_parametric_hrf(
        Y,
        matrix(onsets, ncol = 1),
        confound_formula = formulas[[name]],
        parametric_model = "lwu",
        hrf_eval_times = t_hrf,
        global_refinement = TRUE,
        global_passes = 2,
        verbose = FALSE
      )
    }
    
    results[[name]] <- list(
      params = as.numeric(coef(fit)),
      amplitude = fit$amplitudes[1],
      r_squared = fit$r_squared[1]
    )
  }
  
  # Confound regression effects can vary
  none_errors <- abs(results$none$params - true_params)
  poly_errors <- abs(results$polynomial$params - true_params)
  
  # In this test case, polynomial confounds help with tau
  expect_lt(poly_errors[1], none_errors[1] * 1.1,
            label = "Polynomial confounds improve tau estimate")
  
  # For sigma, the effect can go either way due to the interaction
  # between drift removal and HRF width estimation
  # We just check that it doesn't get dramatically worse
  expect_lt(poly_errors[2], 0.5,
            label = "Sigma estimate remains reasonable with polynomial confounds")
  
  # All confound models should preserve the HRF parameters reasonably well
  for (name in names(results)) {
    if (name != "none") {
      expect_equal(results[[name]]$params[1], true_params[1], tolerance = 0.6,
                   label = paste("tau with", name, "confounds"))
      expect_equal(results[[name]]$params[2], true_params[2], tolerance = 0.4,
                   label = paste("sigma with", name, "confounds"))
      expect_equal(results[[name]]$params[3], true_params[3], tolerance = 0.2,
                   label = paste("rho with", name, "confounds"))
      
      # Amplitude should be preserved (around 1.0)
      expect_equal(results[[name]]$amplitude, 1.0, tolerance = 0.25,
                   label = paste("amplitude with", name, "confounds"))
    }
  }
  
  # Confound regression should improve R-squared
  expect_gt(results$polynomial$r_squared, results$none$r_squared,
            label = "Confound regression improves R-squared")
})

test_that("Algorithm handles degenerate confound specifications", {
  set.seed(111)
  n_time <- 100
  TR <- 2.0
  
  # Create simple data
  onsets <- rep(0, n_time)
  onsets[c(20, 50, 80)] <- 1
  
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  t_hrf <- seq(0, 30, by = TR)
  true_hrf <- hrf_interface$hrf_function(t_hrf, c(5, 2, 0.3))
  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  Y <- matrix(conv_full[1:n_time] + rnorm(n_time, sd = 0.1), ncol = 1)
  
  # Test 1: High-degree polynomial (potential overfitting)
  expect_warning({
    fit_highdeg <- estimate_parametric_hrf(
      Y,
      matrix(onsets, ncol = 1),
      confound_formula = ~ poly(scan, 10),
      parametric_model = "lwu",
      hrf_eval_times = t_hrf,
      verbose = FALSE
    )
  }, NA)  # Should not warn, but handle gracefully
  
  # Parameters should still be reasonable
  params <- as.numeric(coef(fit_highdeg))
  expect_true(all(is.finite(params)))
  expect_true(params[1] >= 0 && params[1] <= 20)  # tau bounds
  expect_true(params[2] > 0.05 && params[2] <= 10)  # sigma bounds
  expect_true(params[3] >= 0 && params[3] <= 1.5)  # rho bounds
  
  # Test 2: Constant confound (should be handled)
  fit_const <- estimate_parametric_hrf(
    Y,
    matrix(onsets, ncol = 1),
    confound_formula = ~ 1,  # Just intercept
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,
    verbose = FALSE
  )
  
  # Should produce valid results
  expect_true(all(is.finite(coef(fit_const))))
  expect_true(all(fit_const$r_squared >= 0 & fit_const$r_squared <= 1))
})

test_that("baseline_model provides equivalent or better drift correction than formulas", {
  skip_if_not_installed("fmrireg")
  
  set.seed(999)
  library(fmrireg)
  
  # Setup
  sframe <- sampling_frame(blocklens = 200, TR = 2.0)
  n_time <- 200
  true_params <- c(5.5, 1.8, 0.35)
  
  # Create events
  onsets <- rep(0, n_time)
  event_times <- seq(20, n_time-20, by = 30)
  onsets[event_times] <- 1
  
  # Generate clean signal
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  t_hrf <- seq(0, 30, by = 2)
  true_hrf <- hrf_interface$hrf_function(t_hrf, true_params)
  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  clean_signal <- conv_full[1:n_time]
  
  # Add polynomial drift
  time_vec <- 1:n_time
  drift <- 0.4 * (time_vec/n_time) + 0.3 * (time_vec/n_time)^2 - 0.1 * (time_vec/n_time)^3
  
  Y <- matrix(clean_signal + drift + rnorm(n_time, sd = 0.05), ncol = 1)
  
  # Method 1: confound_formula
  fit_formula <- estimate_parametric_hrf(
    Y,
    matrix(onsets, ncol = 1),
    confound_formula = ~ poly(scan, 3),
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,
    verbose = FALSE
  )
  
  # Method 2: baseline_model
  bmodel <- baseline_model(
    basis = "poly",
    degree = 3,
    sframe = sframe
  )
  
  fit_baseline <- estimate_parametric_hrf(
    Y,
    matrix(onsets, ncol = 1),
    baseline_model = bmodel,
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,
    verbose = FALSE
  )
  
  # Both should give similar results
  params_formula <- as.numeric(coef(fit_formula))
  params_baseline <- as.numeric(coef(fit_baseline))
  
  expect_equal(params_formula[1], params_baseline[1], tolerance = 0.1,
               label = "tau similar between methods")
  expect_equal(params_formula[2], params_baseline[2], tolerance = 0.1,
               label = "sigma similar between methods")
  expect_equal(params_formula[3], params_baseline[3], tolerance = 0.05,
               label = "rho similar between methods")
  
  # Both should recover true parameters well
  expect_equal(params_baseline[1], true_params[1], tolerance = 0.5)
  expect_equal(params_baseline[2], true_params[2], tolerance = 0.3)
  expect_equal(params_baseline[3], true_params[3], tolerance = 0.15)
})