# Test bounds enforcement throughout the fitting process

test_that("parameter bounds are enforced during initial estimation", {
  set.seed(123)
  
  # Create test data
  n_time <- 100
  n_voxels <- 5
  fmri_data <- matrix(rnorm(n_time * n_voxels), nrow = n_time, ncol = n_voxels)
  
  # Create event design with some events
  event_design <- matrix(0, nrow = n_time, ncol = 1)
  event_design[seq(10, 90, by = 20), 1] <- 1
  
  # Test with restrictive bounds
  restrictive_bounds <- list(
    lower = c(2, 1, 0.1),    # tau >= 2, sigma >= 1, rho >= 0.1
    upper = c(8, 3, 0.5)     # tau <= 8, sigma <= 3, rho <= 0.5
  )
  
  # Fit with restrictive bounds
  fit_restricted <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    theta_bounds = restrictive_bounds,
    verbose = FALSE
  )
  
  # Check all parameters are within bounds
  params <- fit_restricted$estimated_parameters
  
  # tau bounds
  expect_true(all(params[, 1] >= restrictive_bounds$lower[1]))
  expect_true(all(params[, 1] <= restrictive_bounds$upper[1]))
  
  # sigma bounds
  expect_true(all(params[, 2] >= restrictive_bounds$lower[2]))
  expect_true(all(params[, 2] <= restrictive_bounds$upper[2]))
  
  # rho bounds
  expect_true(all(params[, 3] >= restrictive_bounds$lower[3]))
  expect_true(all(params[, 3] <= restrictive_bounds$upper[3]))
})

test_that("bounds enforcement works with extreme initial seeds", {
  set.seed(456)
  
  # Create test data
  n_time <- 50
  fmri_data <- matrix(rnorm(n_time * 2), nrow = n_time, ncol = 2)
  event_design <- matrix(rbinom(n_time, 1, 0.1), ncol = 1)
  
  # Standard bounds
  standard_bounds <- list(
    lower = c(0, 0.1, 0),
    upper = c(20, 10, 1.5)
  )
  
  # Test with seed outside bounds (too low)
  low_seed <- c(-5, -2, -1)  # All parameters below lower bounds
  
  fit_low <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    theta_seed = low_seed,
    theta_bounds = standard_bounds,
    verbose = FALSE
  )
  
  # Test with seed outside bounds (too high)
  high_seed <- c(30, 15, 3)  # All parameters above upper bounds
  
  fit_high <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    theta_seed = high_seed,
    theta_bounds = standard_bounds,
    verbose = FALSE
  )
  
  # Both should produce valid parameters within bounds
  for (fit in list(fit_low, fit_high)) {
    params <- fit$estimated_parameters
    expect_true(all(params[, 1] >= standard_bounds$lower[1]))
    expect_true(all(params[, 1] <= standard_bounds$upper[1]))
    expect_true(all(params[, 2] >= standard_bounds$lower[2]))
    expect_true(all(params[, 2] <= standard_bounds$upper[2]))
    expect_true(all(params[, 3] >= standard_bounds$lower[3]))
    expect_true(all(params[, 3] <= standard_bounds$upper[3]))
  }
})

test_that("bounds are enforced during refinement stages", {
  skip_if_not(requireNamespace("fmrireg", quietly = TRUE),
              "fmrireg not available")
  
  set.seed(789)
  
  # Create challenging data that might push parameters to boundaries
  n_time <- 100
  n_voxels <- 10
  
  # Generate data with varying SNR to trigger refinement
  fmri_data <- matrix(0, nrow = n_time, ncol = n_voxels)
  for (v in 1:n_voxels) {
    # Add different noise levels
    noise_level <- 0.1 + (v - 1) * 0.2
    fmri_data[, v] <- rnorm(n_time, sd = noise_level)
  }
  
  # Add signal
  event_times <- seq(10, 90, by = 20)
  for (t in event_times) {
    if (t + 10 <= n_time) {
      # Add a simple response
      fmri_data[t:(t+10), ] <- fmri_data[t:(t+10), ] + 
        matrix(rep(seq(0, 1, length.out = 11), n_voxels), ncol = n_voxels)
    }
  }
  
  event_design <- matrix(0, nrow = n_time, ncol = 1)
  event_design[event_times, 1] <- 1
  
  # Use tight bounds that refinement might want to violate
  tight_bounds <- list(
    lower = c(3, 1.5, 0.2),
    upper = c(7, 2.5, 0.4)
  )
  
  # Fit with aggressive refinement
  fit_refined <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    theta_bounds = tight_bounds,
    tiered_refinement = "aggressive",
    global_refinement = TRUE,
    verbose = FALSE
  )
  
  # Check bounds are respected even after refinement
  params <- fit_refined$estimated_parameters
  
  expect_true(all(params[, 1] >= tight_bounds$lower[1] - 1e-10))  # Allow tiny numerical errors
  expect_true(all(params[, 1] <= tight_bounds$upper[1] + 1e-10))
  expect_true(all(params[, 2] >= tight_bounds$lower[2] - 1e-10))
  expect_true(all(params[, 2] <= tight_bounds$upper[2] + 1e-10))
  expect_true(all(params[, 3] >= tight_bounds$lower[3] - 1e-10))
  expect_true(all(params[, 3] <= tight_bounds$upper[3] + 1e-10))
})

test_that("sigma lower bound of 0.051 is enforced for numerical stability", {
  set.seed(321)
  
  # Create test data
  n_time <- 50
  fmri_data <- matrix(rnorm(n_time * 3), nrow = n_time, ncol = 3)
  event_design <- matrix(rbinom(n_time, 1, 0.15), ncol = 1)
  
  # Try to set sigma bound below 0.051
  invalid_bounds <- list(
    lower = c(0, 0.01, 0),    # sigma lower bound is too small
    upper = c(20, 10, 1.5)
  )
  
  # The validation should warn about this
  expect_warning(
    fit <- estimate_parametric_hrf(
      fmri_data = fmri_data,
      event_model = event_design,
      theta_bounds = invalid_bounds,
      verbose = FALSE
    ),
    regexp = "sigma.*0\\.05|physiologically implausible",
    ignore.case = TRUE
  )
  
  # Even with the invalid bound request, sigma should never go below 0.051
  if (exists("fit")) {
    params <- fit$estimated_parameters
    expect_true(all(params[, 2] >= 0.051))
  }
})

test_that("bounds enforcement handles edge case at boundaries", {
  set.seed(654)
  
  # Create synthetic data designed to push parameters to boundaries
  n_time <- 80
  n_voxels <- 4
  
  # Create event design
  event_design <- matrix(0, nrow = n_time, ncol = 1)
  event_design[c(10, 30, 50, 70), 1] <- 1
  
  # Generate HRF with parameters at boundaries
  true_tau <- 2      # At lower bound
  true_sigma <- 10   # At upper bound
  true_rho <- 0      # At lower bound
  
  # Generate true HRF
  t_hrf <- seq(0, 30, length.out = 61)
  # Need to provide bounds for the function
  test_bounds <- list(lower = c(0, 0.051, 0), upper = c(20, 10, 1.5))
  true_hrf <- fmriparametric:::.lwu_hrf_function(
    t_hrf, 
    c(true_tau, true_sigma, true_rho),
    bounds = test_bounds
  )
  
  # Convolve with events to create signal
  fmri_data <- matrix(0, nrow = n_time, ncol = n_voxels)
  for (v in 1:n_voxels) {
    # Use fast_batch_convolution to create signal
    signal <- fmriparametric:::.fast_batch_convolution(
      matrix(event_design[, 1], ncol = 1),
      matrix(true_hrf, ncol = 1),
      n_time
    )
    
    # Add noise
    fmri_data[, v] <- signal[, 1] + rnorm(n_time, sd = 0.1)
  }
  
  # Set bounds that match the true parameters at boundaries
  boundary_bounds <- list(
    lower = c(2, 0.1, 0),
    upper = c(20, 10, 1.5)
  )
  
  # Fit the model
  fit_boundary <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    theta_bounds = boundary_bounds,
    theta_seed = c(5, 5, 0.5),  # Start away from boundaries
    verbose = FALSE
  )
  
  # Check that the algorithm can find parameters at boundaries
  params <- fit_boundary$estimated_parameters
  
  # At least some parameters should be close to boundaries
  # (within 10% of the boundary range)
  tau_range <- boundary_bounds$upper[1] - boundary_bounds$lower[1]
  sigma_range <- boundary_bounds$upper[2] - boundary_bounds$lower[2]
  rho_range <- boundary_bounds$upper[3] - boundary_bounds$lower[3]
  
  at_lower_tau <- any(params[, 1] < boundary_bounds$lower[1] + 0.1 * tau_range)
  at_upper_sigma <- any(params[, 2] > boundary_bounds$upper[2] - 0.1 * sigma_range)
  at_lower_rho <- any(params[, 3] < boundary_bounds$lower[3] + 0.1 * rho_range)
  
  # At least one parameter should be near a boundary
  expect_true(at_lower_tau || at_upper_sigma || at_lower_rho,
              info = "Algorithm should be able to find solutions near boundaries")
})

# Test to verify bounds are passed correctly through the pipeline
test_that("user-specified bounds propagate through all stages", {
  # This test addresses Gemini's critical finding about inconsistent bounds
  
  set.seed(999)
  n_time <- 60
  fmri_data <- matrix(rnorm(n_time * 2), nrow = n_time, ncol = 2)
  event_design <- matrix(rbinom(n_time, 1, 0.1), ncol = 1)
  
  # Create custom bounds that differ from defaults
  custom_bounds <- list(
    lower = c(1, 0.5, 0.05),   # Different from default c(0, 0.051, 0)
    upper = c(15, 8, 1.2)      # Different from default c(20, 10, 1.5)
  )
  
  # Trace that bounds are used by adding a custom warning in validation
  original_validate <- fmriparametric:::.validate_theta_bounds
  
  # Run fit
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    theta_bounds = custom_bounds,
    verbose = FALSE
  )
  
  # Verify parameters respect custom bounds, not defaults
  params <- fit$estimated_parameters
  
  # Should respect custom bounds
  expect_true(all(params[, 1] >= custom_bounds$lower[1]))
  expect_true(all(params[, 1] <= custom_bounds$upper[1]))
  expect_true(all(params[, 2] >= custom_bounds$lower[2]))
  expect_true(all(params[, 2] <= custom_bounds$upper[2]))
  expect_true(all(params[, 3] >= custom_bounds$lower[3]))
  expect_true(all(params[, 3] <= custom_bounds$upper[3]))
  
  # If using default bounds, tau could be < 1, which would indicate a bug
  expect_false(any(params[, 1] < custom_bounds$lower[1]),
               info = "Parameters should not violate user-specified lower bounds")
})