library(testthat)

context("Global re-centering functionality")

# Test scenario 1: Basic re-centering behavior
test_that("global re-centering improves estimates", {
  set.seed(123)
  n_time <- 100
  n_vox <- 20
  
  # True parameters
  true_params <- c(tau = 6, sigma = 2.5, rho = 0.35)
  
  # Create simple HRF interface
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      # Simple Gaussian-like HRF with derivatives
      hrf_val <- exp(-(t_hrf - theta0[1])^2 / (2 * theta0[2]^2))
      d_tau <- (t_hrf - theta0[1]) * hrf_val / theta0[2]^2
      d_sigma = (t_hrf - theta0[1])^2 * hrf_val / theta0[2]^3
      d_rho <- -0.3 * exp(-(t_hrf - theta0[1] - 2*theta0[2])^2 / (2 * (1.6*theta0[2])^2))
      cbind(hrf_val, d_tau, d_sigma, d_rho)
    }
  )
  
  # Generate stimulus
  S <- matrix(0, nrow = n_time, ncol = 1)
  S[seq(10, n_time, by = 20), 1] <- 1
  
  # Generate data with true HRF
  t_hrf <- 0:15
  true_hrf <- exp(-(t_hrf - true_params[1])^2 / (2 * true_params[2]^2))
  Y <- matrix(0, nrow = n_time, ncol = n_vox)
  
  for (v in 1:n_vox) {
    # Add some variability across voxels
    vox_params <- true_params * runif(3, 0.9, 1.1)
    vox_hrf <- exp(-(t_hrf - vox_params[1])^2 / (2 * vox_params[2]^2))
    conv_signal <- stats::filter(S[,1], vox_hrf, sides = 1)
    conv_signal[is.na(conv_signal)] <- 0
    Y[, v] <- conv_signal + rnorm(n_time, sd = 0.1)
  }
  
  # Start with poor seed
  poor_seed <- c(tau = 8, sigma = 4, rho = 0.6)
  
  # Single pass result
  res_single <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:n_time,
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_interface,
    theta_seed = poor_seed,
    theta_bounds = list(lower = c(1, 0.5, 0), upper = c(15, 10, 1.5)),
    recenter_global_passes = 1,
    compute_residuals = FALSE,
    verbose = FALSE
  )
  
  # Multi-pass result
  res_multi <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:n_time,
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_interface,
    theta_seed = poor_seed,
    theta_bounds = list(lower = c(1, 0.5, 0), upper = c(15, 10, 1.5)),
    recenter_global_passes = 3,
    compute_residuals = FALSE,
    verbose = FALSE
  )
  
  # Check that R² improved
  expect_true(mean(res_multi$r_squared) > mean(res_single$r_squared))
  
  # Check convergence info
  expect_true(length(res_multi$convergence_info$trajectory) > 1)
  expect_true(is.numeric(res_multi$convergence_info$final_global_theta))
  
  # Check that parameters moved closer to truth
  error_single <- mean(abs(colMeans(res_single$theta_hat) - true_params))
  error_multi <- mean(abs(colMeans(res_multi$theta_hat) - true_params))
  expect_true(error_multi < error_single)
})

# Test scenario 2: Convergence detection
test_that("re-centering detects convergence", {
  set.seed(456)
  
  # Create interface that converges quickly
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      # Interface that pushes toward fixed point
      target <- c(6, 2.5, 0.35)
      hrf_val <- rep(1, length(t_hrf))
      # Derivatives point toward target
      d_tau <- 0.5 * (target[1] - theta0[1]) * hrf_val
      d_sigma <- 0.5 * (target[2] - theta0[2]) * hrf_val
      d_rho <- 0.5 * (target[3] - theta0[3]) * hrf_val
      cbind(hrf_val, d_tau, d_sigma, d_rho)
    }
  )
  
  Y <- matrix(rnorm(50), nrow = 10, ncol = 5)
  S <- matrix(c(1, rep(0, 9)), ncol = 1)
  
  res <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:10,
    hrf_eval_times = 0:5,
    hrf_interface = hrf_interface,
    theta_seed = c(5, 2, 0.3),
    theta_bounds = list(lower = c(1, 0.5, 0), upper = c(10, 5, 1)),
    recenter_global_passes = 10,  # Allow many iterations
    recenter_epsilon = 0.01,
    compute_residuals = FALSE,
    verbose = FALSE
  )
  
  # Should converge before max iterations
  expect_true(res$convergence_info$converged)
  expect_true(res$convergence_info$n_iterations < 10)
})

# Test scenario 3: Data-driven initialization
test_that("data_median initialization works correctly", {
  skip_if_not_installed("fmrireg", minimum_version = "0.2.0")
  
  set.seed(789)
  n_time <- 50
  n_vox <- 30
  
  # Create data with clear structure
  Y <- matrix(0, nrow = n_time, ncol = n_vox)
  S <- matrix(0, nrow = n_time, ncol = 1)
  S[c(5, 15, 25, 35), 1] <- 1
  
  # Good voxels with consistent parameters
  good_params <- c(tau = 5, sigma = 2, rho = 0.3)
  for (v in 1:20) {
    hrf <- exp(-(0:10 - good_params[1])^2 / (2 * good_params[2]^2))
    signal <- stats::filter(S[,1], hrf, sides = 1)
    signal[is.na(signal)] <- 0
    Y[, v] <- signal + rnorm(n_time, sd = 0.05)
  }
  
  # Bad voxels with noise
  Y[, 21:30] <- matrix(rnorm(n_time * 10, sd = 1), ncol = 10)
  
  # Mock event model
  event_model <- list(
    terms = list(S),
    sampling_rate = 1
  )
  class(event_model) <- c("event_model", "list")
  
  # Mock fmri data
  fmri_data <- list(
    data = Y,
    sampling_rate = 1,
    TR = 1
  )
  class(fmri_data) <- c("matrix_dataset", "list")
  
  # Test with consolidated estimator
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    theta_seed = "data_median",
    global_refinement = TRUE,
    global_passes = 2,
    verbose = FALSE
  )
  
  expect_s3_class(fit, "parametric_hrf_fit")
  expect_equal(fit$metadata$theta_seed, "data_median")
  
  # Check that good voxels have high R²
  expect_true(mean(fit$r_squared[1:20]) > mean(fit$r_squared[21:30]))
})

# Test scenario 4: R² threshold behavior
test_that("R² threshold correctly filters voxels", {
  set.seed(321)
  
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      cbind(rep(1, length(t_hrf)), 
            rep(0.1, length(t_hrf)),
            rep(0.1, length(t_hrf)), 
            rep(0.1, length(t_hrf)))
    }
  )
  
  # Create data with varying quality
  n_vox <- 50
  Y <- matrix(0, nrow = 20, ncol = n_vox)
  S <- matrix(rbinom(20, 1, 0.2), ncol = 1)
  
  # First half: good signal
  Y[, 1:25] <- 2 * S %*% t(rep(1, 25)) + matrix(rnorm(20*25, sd = 0.1), ncol = 25)
  # Second half: pure noise
  Y[, 26:50] <- matrix(rnorm(20*25, sd = 1), ncol = 25)
  
  res <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:20,
    hrf_eval_times = 0:5,
    hrf_interface = hrf_interface,
    theta_seed = c(5, 2, 0.3),
    theta_bounds = list(lower = c(1, 0.5, 0), upper = c(10, 5, 1)),
    recenter_global_passes = 2,
    r2_threshold = 0.2,
    compute_residuals = FALSE,
    verbose = FALSE
  )
  
  # Check R² distribution
  expect_true(mean(res$r_squared[1:25]) > 0.5)
  expect_true(mean(res$r_squared[26:50]) < 0.2)
  
  # Global theta should be influenced mainly by good voxels
  expect_true(length(res$convergence_info$trajectory) >= 2)
})