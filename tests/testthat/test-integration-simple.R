# Integration tests for core functionality

library(testthat)
library(fmriparametric)

test_that("estimate_parametric_hrf integration with minimal data", {
  set.seed(42)
  
  # Very simple test data
  n_time <- 50
  n_vox <- 5
  
  # Simple synthetic fMRI data
  fmri_data <- matrix(rnorm(n_time * n_vox, mean = 100, sd = 10), 
                      nrow = n_time, ncol = n_vox)
  
  # Simple event model - just a few events
  event_model <- matrix(0, nrow = n_time, ncol = 1)
  event_model[c(10, 25, 40), 1] <- 1
  
  # Should run without errors
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  # Basic structure checks
  expect_s3_class(fit, "parametric_hrf_fit")
  expect_true(is.matrix(coef(fit)))
  expect_equal(nrow(coef(fit)), n_vox)
  expect_equal(ncol(coef(fit)), 3)
  expect_equal(colnames(coef(fit)), c("tau", "sigma", "rho"))
  
  # Parameters should be numeric and finite
  params <- coef(fit)
  expect_true(all(is.finite(params)))
  
  # R-squared should be computed and valid
  expect_true(is.numeric(fit$r_squared))
  expect_equal(length(fit$r_squared), n_vox)
  expect_true(all(fit$r_squared >= 0 & fit$r_squared <= 1))
  
  # Amplitudes should be computed
  expect_true(is.numeric(fit$amplitudes))
  expect_equal(length(fit$amplitudes), n_vox)
  expect_true(all(is.finite(fit$amplitudes)))
})

test_that("estimate_parametric_hrf with realistic synthetic data", {
  # Generate more realistic test data with known HRF
  set.seed(123)
  n_time <- 60
  n_vox <- 3
  
  # Create event times
  event_model <- matrix(0, nrow = n_time, ncol = 1)
  event_model[c(10, 30, 50), 1] <- 1
  
  # Generate synthetic data with true LWU HRF
  true_params <- c(tau = 6, sigma = 2.5, rho = 0.35)
  
  # Simple LWU function for test data generation
  lwu_hrf <- function(t, tau, sigma, rho) {
    main <- exp(-(t - tau)^2 / (2 * sigma^2))
    undershoot <- rho * exp(-(t - tau - 2*sigma)^2 / (2 * (1.6*sigma)^2))
    hrf <- main - undershoot
    hrf[t < 0] <- 0
    hrf
  }
  
  # Generate HRF and convolve with events
  t_hrf <- seq(0, 30, length.out = 61)
  hrf_true <- lwu_hrf(t_hrf, true_params[1], true_params[2], true_params[3])
  conv_signal <- convolve(event_model[, 1], rev(hrf_true), type = "open")[1:n_time]
  
  # Create realistic fMRI data
  fmri_data <- matrix(0, n_time, n_vox)
  for (v in 1:n_vox) {
    amplitude <- 1 + v * 0.5  # Different amplitudes per voxel
    baseline <- 100
    noise_sd <- sd(conv_signal) * amplitude / 4  # SNR = 4
    fmri_data[, v] <- baseline + amplitude * conv_signal + rnorm(n_time, sd = noise_sd)
  }
  
  # Estimate parameters
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  # Basic structure checks
  params <- coef(fit)
  expect_s3_class(fit, "parametric_hrf_fit")
  expect_equal(dim(params), c(n_vox, 3))
  expect_true(all(is.finite(params)))
  
  # R-squared should be reasonably high for good synthetic data
  expect_true(all(fit$r_squared >= 0))
  expect_true(mean(fit$r_squared) > 0.1)  # Should have some explanatory power
  
  # Amplitudes should be positive (since we added positive signal)
  expect_true(all(fit$amplitudes > 0))
})

test_that("S3 methods work correctly", {
  set.seed(123)
  
  # Generate test data
  fmri_data <- matrix(rnorm(30 * 3), nrow = 30, ncol = 3)
  event_model <- matrix(rbinom(30, 1, 0.2), ncol = 1)
  
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  # Test print method
  expect_output(print(fit), "Parametric HRF Fit")
  
  # Test coef method  
  params <- coef(fit)
  expect_true(is.matrix(params))
  expect_equal(dim(params), c(3, 3))
  
  # Test summary method
  summary_fit <- summary(fit)
  expect_s3_class(summary_fit, "summary_parametric_hrf_fit")
})

test_that("estimation handles edge cases gracefully", {
  # Test with very few events
  set.seed(789)
  n_time <- 30
  fmri_data <- matrix(rnorm(30 * 2, mean = 100), nrow = 30, ncol = 2)
  event_model <- matrix(0, nrow = 30, ncol = 1)
  event_model[15, 1] <- 1  # Only one event
  
  expect_no_error({
    fit1 <- estimate_parametric_hrf(
      fmri_data = fmri_data,
      event_model = event_model,
      parametric_hrf = "lwu",
      verbose = FALSE
    )
  })
  
  expect_s3_class(fit1, "parametric_hrf_fit")
  expect_equal(nrow(coef(fit1)), 2)
  expect_true(all(is.finite(coef(fit1))))
})

test_that("basic parameter bounds are reasonable", {
  # Test with a controlled example to verify basic functionality
  set.seed(456)
  
  n_time <- 40
  n_vox <- 2
  
  # Create data with some signal
  event_model <- matrix(0, nrow = n_time, ncol = 1)
  event_model[c(10, 25), 1] <- 1
  
  fmri_data <- matrix(100, n_time, n_vox)
  # Add a simple response pattern
  for (i in 1:n_vox) {
    fmri_data[11:15, i] <- fmri_data[11:15, i] + 2  # Simple response after event
    fmri_data[26:30, i] <- fmri_data[26:30, i] + 1.5  # Simple response after event
    fmri_data[, i] <- fmri_data[, i] + rnorm(n_time, sd = 0.5)  # Add noise
  }
  
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  params <- coef(fit)
  
  # Basic checks - parameters should be finite and have reasonable magnitudes
  expect_true(all(is.finite(params)))
  expect_true(all(abs(params) < 1000))  # No extreme values
  
  # R-squared should be reasonable for data with signal
  expect_true(all(fit$r_squared >= 0))
  expect_true(any(fit$r_squared > 0.1))  # At least some explanatory power
})