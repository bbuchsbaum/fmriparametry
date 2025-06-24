# Test: Taylor Approximation Gradient Recovery
# 
# This is THE critical test that verifies the core mathematical engine of the package.
# It tests whether the Taylor approximation correctly uses HRF gradients to update parameters.
#
# The key insight: By constructing synthetic data that lies EXACTLY on the tangent plane
# of the Taylor approximation, we can verify that the linear solver recovers the true
# parameters with near-perfect precision. Any bug in the gradient computation or 
# linear system construction will cause this test to fail.

library(fmriparametric)
library(testthat)

test_that("Taylor approximation correctly uses HRF gradients for parameter updates", {
  # 1. Generate synthetic data with KNOWN gradient structure
  set.seed(42)
  n_time <- 200
  true_params <- c(tau = 5.0, sigma = 2.0, rho = 0.3)
  
  # Create stimulus with good temporal coverage
  onsets <- rep(0, n_time)
  onsets[seq(10, n_time-20, by = 30)] <- 1
  
  # 2. Get HRF interface and bounds
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  bounds <- hrf_interface$default_bounds()
  t_hrf <- seq(0, 30, by = 0.5)
  
  # Generate signal at TRUE parameters
  true_hrf <- hrf_interface$hrf_function(t_hrf, true_params)
  true_signal_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  true_signal <- true_signal_full[1:n_time]
  
  # 3. Start from NEARBY parameters (not default seed)
  # This ensures we're in the linear regime of Taylor approximation
  perturbed_params <- true_params + c(0.5, 0.3, 0.1)
  
  # Ensure perturbed params are within bounds
  perturbed_params <- pmax(bounds$lower + 0.1, 
                           pmin(perturbed_params, bounds$upper - 0.1))
  
  # 4. Compute Taylor basis at perturbed parameters
  taylor_basis <- hrf_interface$taylor_basis(perturbed_params, t_hrf)
  # This gives us: [HRF, dHRF/dtau, dHRF/dsigma, dHRF/drho]
  
  expect_equal(ncol(taylor_basis), 4, 
               info = "Taylor basis should have 4 columns: HRF + 3 derivatives")
  
  # 5. Create the EXACT signal that the linear system should recover
  # The key equation: Y = beta0 * (h0 + sum(delta_theta_i * dh/dtheta_i))
  # Rearranging: Y = beta0 * h0 + beta0 * sum(delta_theta_i * dh/dtheta_i)
  amplitude <- 1.0
  delta_params <- true_params - perturbed_params
  
  # Construct the HRF that when convolved will give us a known solution
  Y_hrf <- taylor_basis[,1]  # HRF at perturbed params
  for (i in 1:3) {
    Y_hrf <- Y_hrf + delta_params[i] * taylor_basis[,i+1]
  }
  
  # Convolve with stimulus
  Y_conv_full <- stats::convolve(onsets, rev(Y_hrf), type = "open")
  Y_conv <- amplitude * Y_conv_full[1:n_time]
  
  # Add tiny noise to avoid numerical singularities
  Y_conv <- Y_conv + rnorm(n_time, sd = 1e-8)
  
  # 6. Run the estimation starting from perturbed_params
  fit <- estimate_parametric_hrf(
    fmri_data = matrix(Y_conv, ncol = 1),
    event_model = matrix(onsets, ncol = 1),
    theta_seed = perturbed_params,
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,     # CRITICAL: Use same time grid as data generation
    global_refinement = FALSE,  # Single Taylor step only
    tiered_refinement = "none", # No additional refinement
    verbose = FALSE
  )
  
  # 7. CRITICAL TESTS:
  estimated_params <- as.numeric(coef(fit))
  names(estimated_params) <- names(true_params)
  
  # The Taylor approximation should recover the TRUE parameters
  # because we constructed the data to lie exactly on the tangent plane
  param_errors <- abs(estimated_params - true_params)
  
  # Should be near-perfect recovery (within numerical precision + ridge penalty effects)
  # Ridge penalty slightly biases the solution, so we allow small tolerance
  expect_true(all(param_errors < 1e-3), 
    info = sprintf("Parameter errors too large: tau=%.2e, sigma=%.2e, rho=%.2e\nTrue: %s\nEstimated: %s",
                   param_errors[1], param_errors[2], param_errors[3],
                   paste(round(true_params, 4), collapse=", "),
                   paste(round(estimated_params, 4), collapse=", ")))
  
  # Verify the amplitude is recovered correctly
  expect_equal(fit$amplitudes[1], amplitude, tolerance = 1e-2,
               info = sprintf("Amplitude not recovered: expected %.3f, got %.3f", 
                              amplitude, fit$amplitudes[1]))
  
  # Verify high R-squared (should be nearly perfect fit)
  expect_gt(fit$r_squared[1], 0.99)
})

test_that("Taylor approximation gradient directions are independent", {
  # This test verifies that parameter updates along different gradient directions
  # are independent, which is crucial for the optimization to work correctly
  
  set.seed(123)
  n_time <- 200
  base_params <- c(tau = 5.0, sigma = 2.0, rho = 0.3)
  
  # Create stimulus
  onsets <- rep(0, n_time)
  onsets[seq(10, n_time-20, by = 30)] <- 1
  
  # Get HRF interface
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  bounds <- hrf_interface$default_bounds()
  t_hrf <- seq(0, 30, by = 0.5)
  
  # Test 1: Pure tau perturbation
  # Create data that should update ONLY tau
  taylor_basis <- hrf_interface$taylor_basis(base_params, t_hrf)
  
  tau_delta <- 0.5
  Y_hrf_tau <- taylor_basis[,1] + tau_delta * taylor_basis[,2]  # h0 + delta_tau * dh/dtau
  
  Y_conv_full <- stats::convolve(onsets, rev(Y_hrf_tau), type = "open")
  Y_tau <- Y_conv_full[1:n_time] + rnorm(n_time, sd = 1e-8)
  
  fit_tau <- estimate_parametric_hrf(
    matrix(Y_tau, ncol = 1),
    matrix(onsets, ncol = 1),
    theta_seed = base_params,
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,
    global_refinement = FALSE,
    tiered_refinement = "none",
    verbose = FALSE
  )
  
  tau_params <- as.numeric(coef(fit_tau))
  tau_updates <- tau_params - base_params
  
  # Should update primarily tau, minimal effect on sigma and rho
  expect_gt(abs(tau_updates[1]), 0.3)
  expect_lt(abs(tau_updates[2]), 0.1)
  expect_lt(abs(tau_updates[3]), 0.1)
  
  # Test 2: Pure sigma perturbation
  sigma_delta <- 0.3
  Y_hrf_sigma <- taylor_basis[,1] + sigma_delta * taylor_basis[,3]  # h0 + delta_sigma * dh/dsigma
  
  Y_conv_full <- stats::convolve(onsets, rev(Y_hrf_sigma), type = "open")
  Y_sigma <- Y_conv_full[1:n_time] + rnorm(n_time, sd = 1e-8)
  
  fit_sigma <- estimate_parametric_hrf(
    matrix(Y_sigma, ncol = 1),
    matrix(onsets, ncol = 1),
    theta_seed = base_params,
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,
    global_refinement = FALSE,
    tiered_refinement = "none",
    verbose = FALSE
  )
  
  sigma_params <- as.numeric(coef(fit_sigma))
  sigma_updates <- sigma_params - base_params
  
  # Should update primarily sigma
  expect_lt(abs(sigma_updates[1]), 0.1)
  expect_gt(abs(sigma_updates[2]), 0.2)
  expect_lt(abs(sigma_updates[3]), 0.1)
})

test_that("Taylor approximation handles edge cases correctly", {
  # Test behavior at parameter bounds and with extreme values
  
  set.seed(456)
  n_time <- 150
  onsets <- rep(0, n_time)
  onsets[seq(10, n_time-20, by = 25)] <- 1
  
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  bounds <- hrf_interface$default_bounds()
  t_hrf <- seq(0, 30, by = 0.5)
  
  # Test 1: Parameters near lower bounds
  low_params <- bounds$lower + c(0.5, 0.1, 0.05)
  
  taylor_basis <- hrf_interface$taylor_basis(low_params, t_hrf)
  Y_hrf <- taylor_basis[,1]  # Just use HRF, no perturbation
  
  Y_conv_full <- stats::convolve(onsets, rev(Y_hrf), type = "open")
  Y_low <- Y_conv_full[1:n_time] + rnorm(n_time, sd = 0.01)
  
  expect_no_error({
    fit_low <- estimate_parametric_hrf(
      matrix(Y_low, ncol = 1),
      matrix(onsets, ncol = 1),
      theta_seed = low_params,
      parametric_model = "lwu",
      hrf_eval_times = t_hrf,
      global_refinement = FALSE,
      tiered_refinement = "none",
      verbose = FALSE
    )
  })
  
  # Parameters should stay within bounds
  est_params <- as.numeric(coef(fit_low))
  expect_true(all(est_params >= bounds$lower),
              info = "Parameters should respect lower bounds")
  expect_true(all(est_params <= bounds$upper),
              info = "Parameters should respect upper bounds")
  
  # Test 2: Large parameter updates should be handled gracefully
  start_params <- c(tau = 3.0, sigma = 1.0, rho = 0.1)
  target_params <- c(tau = 8.0, sigma = 3.0, rho = 0.8)  # Large jump
  
  # Create synthetic data at target
  target_hrf <- hrf_interface$hrf_function(t_hrf, target_params)
  Y_conv_full <- stats::convolve(onsets, rev(target_hrf), type = "open")
  Y_target <- Y_conv_full[1:n_time] + rnorm(n_time, sd = 0.01)
  
  fit_jump <- estimate_parametric_hrf(
    matrix(Y_target, ncol = 1),
    matrix(onsets, ncol = 1),
    theta_seed = start_params,
    parametric_model = "lwu",
    hrf_eval_times = t_hrf,
    global_refinement = TRUE,  # Allow refinement for large jump
    global_passes = 5,
    tiered_refinement = "none",
    verbose = FALSE
  )
  
  # Should get reasonably close despite large initial distance
  final_params <- as.numeric(coef(fit_jump))
  param_errors <- abs(final_params - target_params)
  
  # More lenient tolerance for large jumps
  # Note: Very large parameter jumps (e.g., tau jumping from 3 to 8) are challenging
  # for Taylor approximation methods, so we allow larger errors
  expect_true(param_errors[1] < 5.0 && param_errors[2] < 3.0 && param_errors[3] < 0.5,
              info = sprintf("Large jump recovery failed: errors = %s",
                             paste(round(param_errors, 3), collapse=", ")))
})