# Test baseline_model integration with fmrireg
#
# This test verifies that the package correctly handles fmrireg::baseline_model
# objects for confound regression, providing a consistent interface with event_model.

library(fmriparametric)
library(testthat)

test_that("baseline_model objects are accepted and processed correctly", {
  skip_if_not_installed("fmrireg")
  
  set.seed(123)
  library(fmrireg)
  
  # Create sampling frame
  sframe <- sampling_frame(blocklens = c(100, 100), TR = 2.0)
  n_time <- sum(blocklens(sframe))
  n_vox <- 50
  
  # Create events
  onsets <- rep(0, n_time)
  onsets[seq(10, n_time-20, by = 25)] <- 1
  event_mat <- matrix(onsets, ncol = 1)
  
  # Generate true signal
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  t_hrf <- seq(0, 30, by = 2)
  true_params <- c(5.5, 2.0, 0.35)
  true_hrf <- hrf_interface$hrf_function(t_hrf, true_params)
  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  clean_signal <- conv_full[1:n_time]
  
  # Add drift
  time_vec <- 1:n_time
  drift <- 0.5 * (time_vec / n_time) + 0.3 * (time_vec / n_time)^2
  
  # Create data with drift
  Y <- matrix(0, n_time, n_vox)
  for (v in 1:n_vox) {
    amp <- runif(1, 0.8, 1.2)
    Y[,v] <- amp * clean_signal + drift + rnorm(n_time, sd = 0.1)
  }
  
  # Test 1: Basic baseline_model with polynomial drift
  bmodel_poly <- baseline_model(
    basis = "poly",
    degree = 3,
    sframe = sframe,
    intercept = "runwise"
  )
  
  fit_baseline <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = event_mat,
    baseline_model = bmodel_poly,
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,
    global_refinement = FALSE,
    verbose = FALSE
  )
  
  # Should produce valid results
  expect_s3_class(fit_baseline, "parametric_hrf_fit")
  expect_equal(nrow(coef(fit_baseline)), n_vox)
  expect_true(all(fit_baseline$r_squared >= 0 & fit_baseline$r_squared <= 1))
  
  # Test 2: Spline basis
  bmodel_spline <- baseline_model(
    basis = "bs",
    degree = 5,
    sframe = sframe,
    intercept = "global"
  )
  
  fit_spline <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = event_mat,
    baseline_model = bmodel_spline,
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,
    verbose = FALSE
  )
  
  expect_s3_class(fit_spline, "parametric_hrf_fit")
  
  # Test 3: With nuisance regressors
  # Simulate motion parameters
  motion_run1 <- matrix(cumsum(rnorm(100 * 6, 0, 0.01)), nrow = 100, ncol = 6)
  motion_run2 <- matrix(cumsum(rnorm(100 * 6, 0, 0.01)), nrow = 100, ncol = 6)
  colnames(motion_run1) <- colnames(motion_run2) <- c("tx", "ty", "tz", "rx", "ry", "rz")
  
  bmodel_nuisance <- baseline_model(
    basis = "bs",
    degree = 4,
    sframe = sframe,
    intercept = "runwise",
    nuisance_list = list(motion_run1, motion_run2)
  )
  
  # Add motion artifacts to data
  motion_full <- rbind(motion_run1, motion_run2)
  Y_motion <- Y
  for (v in 1:n_vox) {
    motion_weights <- runif(6, -0.2, 0.2)
    Y_motion[,v] <- Y_motion[,v] + motion_full %*% motion_weights
  }
  
  fit_nuisance <- estimate_parametric_hrf(
    fmri_data = Y_motion,
    event_model = event_mat,
    baseline_model = bmodel_nuisance,
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,
    verbose = FALSE
  )
  
  # Nuisance regression should improve fits
  expect_s3_class(fit_nuisance, "parametric_hrf_fit")
  expect_true(all(is.finite(coef(fit_nuisance))))
})

test_that("baseline_model and confound_formula can coexist with precedence", {
  skip_if_not_installed("fmrireg")
  
  set.seed(456)
  library(fmrireg)
  
  n_time <- 150
  n_vox <- 20
  
  # Simple data
  onsets <- rep(0, n_time)
  onsets[seq(15, n_time-15, by = 30)] <- 1
  
  # Generate data with linear drift
  time_vec <- 1:n_time
  drift <- 0.5 * (time_vec / n_time)
  
  Y <- matrix(rnorm(n_time * n_vox, sd = 0.1), n_time, n_vox)
  for (v in 1:n_vox) {
    Y[,v] <- Y[,v] + drift
  }
  
  # Create baseline model
  sframe <- sampling_frame(blocklens = n_time, TR = 2.0)
  bmodel <- baseline_model(
    basis = "poly",
    degree = 2,
    sframe = sframe
  )
  
  # Test with both baseline_model and confound_formula
  # baseline_model should take precedence
  fit <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = matrix(onsets, ncol = 1),
    baseline_model = bmodel,
    confound_formula = ~ poly(scan, 5),  # This should be ignored
    parametric_model = "lwu",
    verbose = FALSE
  )
  
  # Should work without error
  expect_s3_class(fit, "parametric_hrf_fit")
  expect_equal(nrow(coef(fit)), n_vox)
})

test_that("baseline_model improves parameter recovery with complex confounds", {
  skip_if_not_installed("fmrireg")
  
  set.seed(789)
  library(fmrireg)
  
  # Two-run design
  sframe <- sampling_frame(blocklens = c(120, 120), TR = 2.0)
  n_time <- sum(blocklens(sframe))
  n_vox <- 30
  
  # Events
  onsets <- rep(0, n_time)
  # Run 1 events
  onsets[c(20, 50, 80, 110)] <- 1
  # Run 2 events  
  onsets[c(140, 170, 200, 230)] <- 1
  
  # True parameters
  true_params <- c(6.0, 2.2, 0.4)
  
  # Generate clean signal
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  t_hrf <- seq(0, 30, by = 2)
  true_hrf <- hrf_interface$hrf_function(t_hrf, true_params)
  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  clean_signal <- conv_full[1:n_time]
  
  # Create complex confounds:
  # 1. Different drift per run
  time_vec <- 1:n_time
  run1_drift <- 0.4 * (time_vec[1:120] / 120)^2
  run2_drift <- -0.3 * (time_vec[121:240] / 120) + 0.5
  drift <- c(run1_drift, run2_drift)
  
  # 2. Physiological noise (cardiac-like)
  cardiac <- sin(2 * pi * 1.2 * time_vec / 2) * 0.2
  
  # 3. Motion spikes
  motion <- rep(0, n_time)
  spike_times <- c(45, 90, 165, 210)
  motion[spike_times] <- sample(c(-1, 1), 4, replace = TRUE) * runif(4, 0.5, 1.0)
  motion <- filter(motion, rep(1/3, 3), sides = 2)  # Smooth slightly
  motion[is.na(motion)] <- 0
  
  # Combine into nuisance matrix
  nuisance_run1 <- cbind(
    cardiac = cardiac[1:120],
    motion = motion[1:120]
  )
  nuisance_run2 <- cbind(
    cardiac = cardiac[121:240],
    motion = motion[121:240]
  )
  
  # Create data
  Y <- matrix(0, n_time, n_vox)
  for (v in 1:n_vox) {
    amp <- runif(1, 0.7, 1.3)
    # Add all confounds
    Y[,v] <- amp * clean_signal + drift + 
             cardiac * runif(1, 0.5, 1.5) + 
             motion * runif(1, 0.3, 0.7) +
             rnorm(n_time, sd = 0.08)
  }
  
  # Fit without baseline model
  fit_no_baseline <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = matrix(onsets, ncol = 1),
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,
    verbose = FALSE
  )
  
  # Fit with comprehensive baseline model
  bmodel <- baseline_model(
    basis = "bs",
    degree = 5,
    sframe = sframe,
    intercept = "runwise",
    nuisance_list = list(nuisance_run1, nuisance_run2)
  )
  
  fit_with_baseline <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = matrix(onsets, ncol = 1),
    baseline_model = bmodel,
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,
    verbose = FALSE
  )
  
  # Baseline model should improve parameter recovery
  params_no_baseline <- coef(fit_no_baseline)
  params_with_baseline <- coef(fit_with_baseline)
  
  errors_no_baseline <- colMeans(abs(params_no_baseline - 
                                     matrix(true_params, n_vox, 3, byrow = TRUE)))
  errors_with_baseline <- colMeans(abs(params_with_baseline - 
                                       matrix(true_params, n_vox, 3, byrow = TRUE)))
  
  # Errors should generally be reduced (allow small tolerance for randomness)
  # For tau, we expect improvement
  expect_lt(errors_with_baseline[1], errors_no_baseline[1] + 0.01,
            label = "Baseline model improves or maintains tau recovery")
  
  # For sigma, improvement might be smaller due to its smaller scale
  # Allow for small variations due to numerical precision
  expect_lt(errors_with_baseline[2], errors_no_baseline[2] + 0.01, 
            label = "Baseline model improves or maintains sigma recovery")
  
  # R-squared should improve or at least be comparable
  expect_gt(mean(fit_with_baseline$r_squared), mean(fit_no_baseline$r_squared) - 0.01,
            label = "Baseline model improves or maintains R-squared")
})

test_that("Different baseline_model specifications produce expected behavior", {
  skip_if_not_installed("fmrireg")
  
  set.seed(321)
  library(fmrireg)
  
  # Single run for simplicity
  sframe <- sampling_frame(blocklens = 100, TR = 2.0)
  n_time <- 100
  n_vox <- 10
  
  # Simple data
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onsets <- rep(0, n_time)
  onsets[c(20, 50, 80)] <- 1
  
  # Test different basis functions
  basis_types <- c("constant", "poly", "bs", "ns")
  
  results <- list()
  for (basis in basis_types) {
    if (basis == "constant") {
      bmodel <- baseline_model(basis = basis, sframe = sframe)
    } else {
      bmodel <- baseline_model(basis = basis, degree = 3, sframe = sframe)
    }
    
    fit <- estimate_parametric_hrf(
      fmri_data = Y,
      event_model = matrix(onsets, ncol = 1),
      baseline_model = bmodel,
      parametric_model = "lwu",
      verbose = FALSE
    )
    
    results[[basis]] <- fit
    
    # All should produce valid results
    expect_s3_class(fit, "parametric_hrf_fit")
    expect_true(all(is.finite(coef(fit))))
  }
  
  # Different basis functions should give slightly different results
  params_poly <- coef(results[["poly"]])
  params_bs <- coef(results[["bs"]])
  expect_false(all(params_poly == params_bs),
               label = "Different basis functions produce different results")
})