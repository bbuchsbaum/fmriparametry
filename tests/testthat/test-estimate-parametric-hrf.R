test_that("estimate_parametric_hrf works with basic inputs", {
  # Create simple test data
  set.seed(123)
  n_time <- 100
  n_voxels <- 5
  
  # Simulate fMRI data
  fmri_data <- matrix(rnorm(n_time * n_voxels), nrow = n_time, ncol = n_voxels)
  
  # Create simple event design
  event_times <- seq(10, 90, by = 20)
  event_design <- matrix(0, nrow = n_time, ncol = 1)
  event_design[event_times, 1] <- 1
  
  # Test basic functionality
  result <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    parametric_hrf = "lwu",
    verbose = FALSE,

  )
  
  expect_s3_class(result, "parametric_hrf_fit")
  expect_equal(nrow(result$estimated_parameters), n_voxels)
  expect_equal(ncol(result$estimated_parameters), 3) # LWU has 3 parameters
  expect_true(all(is.finite(result$estimated_parameters)))
})

test_that("estimate_parametric_hrf validates inputs correctly", {
  fmri_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
  event_design <- matrix(c(1,0,0,1,0,0,1,0,0,1), ncol = 1)
  
  # Test invalid HRF model
  expect_error(
    estimate_parametric_hrf(fmri_data, event_design, parametric_hrf = "invalid"),
    "HRF model 'invalid' is not registered"
  )
  
  # Test invalid theta_seed length
  expect_error(
    estimate_parametric_hrf(fmri_data, event_design, theta_seed = c(1, 2)),
    "theta_seed.*must have length 3"
  )
  
  # Test invalid theta_bounds structure
  expect_error(
    estimate_parametric_hrf(fmri_data, event_design, 
                           theta_bounds = list(lower = c(1, 2, 3))),
    "theta_bounds must have both"
  )
  
  # Test non-logical verbose
  expect_error(
    estimate_parametric_hrf(fmri_data, event_design, verbose = "yes"),
    "verbose.*must be logical"
  )
})

test_that("estimate_parametric_hrf returns expected structure", {
  set.seed(456)
  fmri_data <- matrix(rnorm(200), nrow = 20, ncol = 10)
  event_design <- matrix(rbinom(20, 1, 0.2), ncol = 1)
  
  result <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    verbose = FALSE
  )
  
  # Check structure - current fit object contains a richer set of fields
  expect_named(
    result,
    c(
      "estimated_parameters", "amplitudes", "parameter_names", "hrf_model",
      "r_squared", "residuals", "parameter_ses", "convergence_info",
      "metadata", "parameters", "convergence", "standard_errors",
      "se_amplitudes", "fit_quality", "refinement_info"
    )
  )
  
  # Check dimensions
  expect_equal(nrow(result$estimated_parameters), ncol(fmri_data))
  expect_equal(length(result$amplitudes), ncol(fmri_data))
  
  # Check metadata
  expect_equal(result$metadata$parametric_hrf, "lwu")
  expect_equal(result$metadata$n_voxels, 10)
  expect_equal(result$metadata$n_timepoints, 20)
})

test_that("estimate_parametric_hrf handles edge cases", {
  # Test with single voxel
  fmri_single <- matrix(rnorm(50), ncol = 1)
  event_design <- matrix(rbinom(50, 1, 0.1), ncol = 1)
  
  result_single <- estimate_parametric_hrf(
    fmri_data = fmri_single,
    event_model = event_design,
    verbose = FALSE
  )
  
  expect_s3_class(result_single, "parametric_hrf_fit")
  expect_equal(nrow(result_single$estimated_parameters), 1)
  
  # Test with no events (should handle gracefully)
  no_events <- matrix(0, nrow = 50, ncol = 1)
  
  expect_warning(
    result_no_events <- estimate_parametric_hrf(
      fmri_data = fmri_single,
      event_model = no_events,
      verbose = FALSE
    ),
    "No events detected|Design matrix may be singular"
  )
})

test_that("estimate_parametric_hrf respects parameter bounds", {
  set.seed(789)
  fmri_data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  event_design <- matrix(c(rep(c(1,0), 10)), ncol = 1)
  
  # Custom bounds
  custom_bounds <- list(
    lower = c(2, 0.5, 0.1),
    upper = c(10, 5, 0.8)
  )
  
  result <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    theta_bounds = custom_bounds,
    verbose = FALSE
  )
  
  # Check all parameters are within bounds
  params <- result$estimated_parameters
  expect_true(all(params[,1] >= custom_bounds$lower[1] & 
                 params[,1] <= custom_bounds$upper[1]))
  expect_true(all(params[,2] >= custom_bounds$lower[2] & 
                 params[,2] <= custom_bounds$upper[2]))
  expect_true(all(params[,3] >= custom_bounds$lower[3] & 
                 params[,3] <= custom_bounds$upper[3]))
})