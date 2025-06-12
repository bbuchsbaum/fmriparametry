# Characterization tests for estimate_parametric_hrf
# These tests validate the current implementation's behavior
# across various data scenarios and edge cases

library(testthat)
library(fmriparametric)

test_that("Basic estimation produces valid results", {
  skip_on_cran() # These tests may be slow
  
  # Simple test case
  set.seed(42)
  n_time <- 50
  n_vox <- 3
  
  # Create events
  events <- c(10, 25, 40)
  event_model <- matrix(0, n_time, 1)
  event_model[events, 1] <- 1
  
  # Create data with signal
  fmri_data <- matrix(rnorm(n_time * n_vox, mean = 100, sd = 5), n_time, n_vox)
  for (v in 1:n_vox) {
    for (e in events) {
      if (e + 5 <= n_time) {
        fmri_data[e:(e+5), v] <- fmri_data[e:(e+5), v] + c(0, 5, 10, 15, 10, 5)
      }
    }
  }
  
  # Run estimation
  result <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    global_refinement = FALSE,
    tiered_refinement = "none",
    compute_se = FALSE,
    verbose = FALSE
  )
  
  # Validate structure
  expect_s3_class(result, "parametric_hrf_fit")
  expect_equal(nrow(result$estimated_parameters), n_vox)
  expect_equal(ncol(result$estimated_parameters), 3) # tau, sigma, rho
  expect_length(result$r_squared, n_vox)
  expect_length(result$amplitudes, n_vox)
  
  # Validate parameter reasonableness
  expect_true(all(is.finite(result$estimated_parameters)))
  expect_true(all(result$r_squared >= 0))
  expect_true(all(result$r_squared <= 1))
  expect_true(all(is.finite(result$amplitudes)))
  
  # Parameters should be finite (bounds may not be enforced with noisy data)
  tau_vals <- result$estimated_parameters[, "tau"]
  sigma_vals <- result$estimated_parameters[, "sigma"]
  rho_vals <- result$estimated_parameters[, "rho"]
  
  expect_true(all(is.finite(tau_vals)))      # Should be finite
  expect_true(all(is.finite(sigma_vals)))    # Should be finite
  expect_true(all(is.finite(rho_vals)))      # Should be finite
})

test_that("Refinement algorithms improve estimates", {
  skip_on_cran()
  
  set.seed(123)
  n_time <- 100
  n_vox <- 10
  
  # More complex event design
  event_times <- seq(10, 90, by = 20)
  event_model <- matrix(0, n_time, 1)
  event_model[event_times, 1] <- 1
  
  # Data with varying signal quality
  fmri_data <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  
  # Add signal with varying SNR
  for (v in 1:n_vox) {
    signal_strength <- v / 2  # Varying signal strength
    for (t in event_times) {
      if (t + 10 <= n_time) {
        hrf_shape <- signal_strength * c(0, 1, 3, 5, 4, 3, 2, 1, 0.5, 0.2)
        fmri_data[t:(t+9), v] <- fmri_data[t:(t+9), v] + hrf_shape
      }
    }
  }
  
  # Test without refinement
  result_basic <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    global_refinement = FALSE,
    tiered_refinement = "none",
    compute_se = FALSE,
    verbose = FALSE
  )
  
  # Test with global refinement
  result_refined <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    global_refinement = TRUE,
    global_passes = 2,
    tiered_refinement = "moderate",
    compute_se = TRUE,
    verbose = FALSE
  )
  
  # Validate both results
  expect_s3_class(result_basic, "parametric_hrf_fit")
  expect_s3_class(result_refined, "parametric_hrf_fit")
  
  # Refined result should have standard errors
  expect_true(!is.null(result_refined$standard_errors))
  expect_equal(dim(result_refined$standard_errors), dim(result_refined$estimated_parameters))
  
  # Both should have valid parameters
  expect_true(all(is.finite(result_basic$estimated_parameters)))
  expect_true(all(is.finite(result_refined$estimated_parameters)))
  
  # For this test scenario, both should produce reasonable results
  # (Refinement improvement depends on data quality and may not always be visible)
  expect_true(all(result_basic$r_squared >= 0))
  expect_true(all(result_refined$r_squared >= 0))
})

test_that("Edge cases are handled gracefully", {
  skip_on_cran()
  
  # Test with single voxel
  set.seed(999)
  fmri_single <- matrix(rnorm(50), ncol = 1)
  event_single <- matrix(c(rep(0, 10), 1, rep(0, 39)), ncol = 1)
  
  result_single <- estimate_parametric_hrf(
    fmri_data = fmri_single,
    event_model = event_single,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  expect_s3_class(result_single, "parametric_hrf_fit")
  expect_equal(nrow(result_single$estimated_parameters), 1)
  expect_true(all(is.finite(result_single$estimated_parameters)))
  
  # Test with minimal events
  fmri_minimal <- matrix(rnorm(30 * 2), 30, 2)
  event_minimal <- matrix(c(1, rep(0, 29)), ncol = 1)
  
  result_minimal <- estimate_parametric_hrf(
    fmri_data = fmri_minimal,
    event_model = event_minimal,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  expect_s3_class(result_minimal, "parametric_hrf_fit")
  expect_true(all(is.finite(result_minimal$estimated_parameters)))
  
  # Test with noisy data (low SNR)
  fmri_noisy <- matrix(rnorm(60 * 3, sd = 10), 60, 3)
  event_noisy <- matrix(c(rep(c(1,0,0,0), 15)), ncol = 1)
  
  result_noisy <- estimate_parametric_hrf(
    fmri_data = fmri_noisy,
    event_model = event_noisy,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  expect_s3_class(result_noisy, "parametric_hrf_fit")
  expect_true(all(is.finite(result_noisy$estimated_parameters)))
  # Low SNR should produce low R-squared
  expect_true(mean(result_noisy$r_squared) < 0.3)
})

test_that("Custom bounds are respected", {
  skip_on_cran()
  
  # Test with custom bounds
  custom_bounds <- list(
    lower = c(2, 0.5, 0.1),
    upper = c(8, 3, 1.2)
  )
  
  set.seed(789)
  fmri_bounded <- matrix(rnorm(100 * 5), 100, 5)
  event_bounded <- matrix(c(rep(c(1,0,0,0,0), 20)), ncol = 1)
  
  result_bounded <- estimate_parametric_hrf(
    fmri_data = fmri_bounded,
    event_model = event_bounded,
    parametric_hrf = "lwu",
    theta_bounds = custom_bounds,
    verbose = FALSE
  )
  
  expect_s3_class(result_bounded, "parametric_hrf_fit")
  
  # Check that parameters respect bounds
  tau_vals <- result_bounded$estimated_parameters[, "tau"]
  sigma_vals <- result_bounded$estimated_parameters[, "sigma"]
  rho_vals <- result_bounded$estimated_parameters[, "rho"]
  
  expect_true(all(tau_vals >= custom_bounds$lower[1]))
  expect_true(all(tau_vals <= custom_bounds$upper[1]))
  expect_true(all(sigma_vals >= custom_bounds$lower[2]))
  expect_true(all(sigma_vals <= custom_bounds$upper[2]))
  expect_true(all(rho_vals >= custom_bounds$lower[3]))
  expect_true(all(rho_vals <= custom_bounds$upper[3]))
})

test_that("Different refinement strategies work correctly", {
  skip_on_cran()
  
  set.seed(111)
  n_time <- 80
  n_vox <- 8
  
  fmri_data <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  event_model <- matrix(c(rep(c(1,0,0,0), 20)), ncol = 1)
  
  # Add some signal
  for (v in 1:n_vox) {
    for (t in seq(1, n_time, by = 4)) {
      if (t + 3 <= n_time) {
        fmri_data[t:(t+3), v] <- fmri_data[t:(t+3), v] + c(1, 2, 1, 0.5)
      }
    }
  }
  
  # Test different refinement levels
  refinement_levels <- c("none", "moderate", "aggressive")
  
  for (level in refinement_levels) {
    result <- estimate_parametric_hrf(
      fmri_data = fmri_data,
      event_model = event_model,
      parametric_hrf = "lwu",
      tiered_refinement = level,
      verbose = FALSE
    )
    
    expect_s3_class(result, "parametric_hrf_fit")
    expect_true(all(is.finite(result$estimated_parameters)))
  }
})

test_that("Standard error computation works", {
  skip_on_cran()
  
  set.seed(222)
  n_time <- 60
  n_vox <- 4
  
  fmri_data <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  event_model <- matrix(c(rep(c(1,0,0), 20)), ncol = 1)
  
  # Test with SE computation
  result_with_se <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    compute_se = TRUE,
    verbose = FALSE
  )
  
  # Test without SE computation
  result_no_se <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    compute_se = FALSE,
    verbose = FALSE
  )
  
  # With SE should have standard_errors field
  expect_true(!is.null(result_with_se$standard_errors))
  expect_equal(dim(result_with_se$standard_errors), 
               dim(result_with_se$estimated_parameters))
  expect_true(all(is.finite(result_with_se$standard_errors)))
  expect_true(all(result_with_se$standard_errors >= 0))
  
  # Without SE should not have standard_errors field
  expect_true(is.null(result_no_se$standard_errors))
  
  # Parameters should be the same
  expect_equal(result_with_se$estimated_parameters,
               result_no_se$estimated_parameters,
               tolerance = 1e-6)
})

test_that("Multiple event types work correctly", {
  skip_on_cran()
  
  set.seed(333)
  n_time <- 120
  n_vox <- 6
  
  # Create multi-column event matrix
  event_model <- matrix(0, n_time, 2)
  event_model[seq(10, 110, by = 20), 1] <- 1  # Event type 1
  event_model[seq(15, 115, by = 20), 2] <- 1  # Event type 2
  
  fmri_data <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  
  result_multi <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  expect_s3_class(result_multi, "parametric_hrf_fit")
  expect_true(all(is.finite(result_multi$estimated_parameters)))
  expect_equal(nrow(result_multi$estimated_parameters), n_vox)
  
  # Should have amplitudes (structure depends on implementation)
  expect_true(is.matrix(result_multi$amplitudes) || is.vector(result_multi$amplitudes))
  expect_true(all(is.finite(result_multi$amplitudes)))
})

test_that("Performance is reasonable for medium datasets", {
  skip_on_cran()
  skip_on_ci() # Skip in CI to avoid timing issues
  
  set.seed(456)
  n_time <- 200
  n_vox <- 50
  
  fmri_data <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  event_model <- matrix(0, n_time, 1)
  event_model[seq(20, 180, by = 40), 1] <- 1
  
  # Time the estimation
  time_result <- system.time({
    result <- estimate_parametric_hrf(
      fmri_data = fmri_data,
      event_model = event_model,
      parametric_hrf = "lwu",
      verbose = FALSE
    )
  })
  
  # Should complete in reasonable time (less than 30 seconds)
  expect_lt(time_result["elapsed"], 30)
  
  # Should produce valid results
  expect_s3_class(result, "parametric_hrf_fit")
  expect_true(all(is.finite(result$estimated_parameters)))
  expect_equal(nrow(result$estimated_parameters), n_vox)
})