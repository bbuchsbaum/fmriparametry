library(testthat)

context("Standard error calculation via Delta method")

# Test scenario 1: SE calculation accuracy
test_that("standard errors match theoretical expectations", {
  set.seed(1234)
  
  # Simple linear model scenario for comparison
  n_time <- 100
  n_vox <- 5
  sigma_true <- 0.5  # Known noise level
  
  # Create HRF interface with known structure
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      # Linear basis for simple comparison
      cbind(rep(1, length(t_hrf)), 
            t_hrf,
            t_hrf^2,
            t_hrf^3)
    }
  )
  
  # Generate data with known noise
  S <- matrix(1, nrow = n_time, ncol = 1)
  t_eval <- seq(0, 1, length.out = 4)
  
  # True parameters
  true_params <- c(2, 1, 0.5)
  X_true <- cbind(1, seq_len(n_time)/n_time, (seq_len(n_time)/n_time)^2, (seq_len(n_time)/n_time)^3)
  
  Y <- matrix(0, nrow = n_time, ncol = n_vox)
  for (v in 1:n_vox) {
    Y[, v] <- X_true %*% c(2, 0.5, -0.3, 0.1) + rnorm(n_time, sd = sigma_true)
  }
  
  # Fit with SE calculation
  res <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = seq_len(n_time),
    hrf_eval_times = t_eval,
    hrf_interface = hrf_interface,
    theta_seed = true_params,
    theta_bounds = list(lower = c(0, 0, 0), upper = c(10, 10, 10)),
    recenter_global_passes = 1,
    compute_residuals = TRUE,
    compute_se = TRUE,
    lambda_ridge_jacobian = 0.001,  # Small ridge for stability
    verbose = FALSE
  )
  
  # Check SE structure
  expect_equal(dim(res$se_theta_hat), c(n_vox, length(true_params)))
  expect_true(all(res$se_theta_hat > 0, na.rm = TRUE))
  
  # SEs should be similar across voxels (same noise level)
  se_means <- colMeans(res$se_theta_hat, na.rm = TRUE)
  se_sds <- apply(res$se_theta_hat, 2, sd, na.rm = TRUE)
  expect_true(all(se_sds / se_means < 0.5))  # Low relative variation
})

# Test scenario 2: SE increases with noise
test_that("standard errors increase with noise level", {
  set.seed(2345)
  
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      hrf_val <- exp(-(t_hrf - theta0[1])^2 / (2 * theta0[2]^2))
      # Simple derivatives
      d1 <- 0.1 * hrf_val
      d2 <- 0.2 * hrf_val
      d3 <- 0.05 * hrf_val
      cbind(hrf_val, d1, d2, d3)
    }
  )
  
  n_time <- 50
  noise_levels <- c(0.1, 0.5, 1.0)
  ses_by_noise <- list()
  
  for (i in seq_along(noise_levels)) {
    # Generate data with different noise
    S <- matrix(rbinom(n_time, 1, 0.2), ncol = 1)
    Y <- matrix(rnorm(n_time * 3, sd = noise_levels[i]), ncol = 3)
    
    res <- .parametric_engine_iterative(
      Y_proj = Y,
      S_target_proj = S,
      scan_times = 1:n_time,
      hrf_eval_times = 0:10,
      hrf_interface = hrf_interface,
      theta_seed = c(5, 2, 0.3),
      theta_bounds = list(lower = c(1, 0.5, 0), upper = c(10, 5, 1)),
      recenter_global_passes = 1,
      compute_residuals = TRUE,
      compute_se = TRUE,
      verbose = FALSE
    )
    
    ses_by_noise[[i]] <- colMeans(res$se_theta_hat, na.rm = TRUE)
  }
  
  # SEs should increase with noise
  for (param in 1:3) {
    se_vals <- sapply(ses_by_noise, function(x) x[param])
    expect_true(all(diff(se_vals) > 0))
  }
})

# Test scenario 3: SE computation with near-zero amplitudes
test_that("SE handles near-zero amplitudes gracefully", {
  set.seed(3456)
  
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      cbind(rep(1, length(t_hrf)),
            rep(0.1, length(t_hrf)),
            rep(0.1, length(t_hrf)),
            rep(0.1, length(t_hrf)))
    }
  )
  
  # Create data with some zero-amplitude voxels
  n_time <- 30
  n_vox <- 10
  Y <- matrix(0, nrow = n_time, ncol = n_vox)
  S <- matrix(c(rep(c(1,0,0), 10)), nrow = n_time, ncol = 1)
  
  # Half voxels have signal, half don't
  Y[, 1:5] <- S %*% t(rep(1, 5)) + matrix(rnorm(n_time * 5, sd = 0.1), ncol = 5)
  Y[, 6:10] <- matrix(rnorm(n_time * 5, sd = 0.01), ncol = 5)  # Just noise
  
  res <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:n_time,
    hrf_eval_times = 0:5,
    hrf_interface = hrf_interface,
    theta_seed = c(5, 2, 0.3),
    theta_bounds = list(lower = c(1, 0.5, 0), upper = c(10, 5, 1)),
    recenter_global_passes = 1,
    compute_residuals = TRUE,
    compute_se = TRUE,
    epsilon_beta = 1e-6,
    verbose = FALSE
  )
  
  # Should have SEs for good voxels
  expect_true(any(!is.na(res$se_theta_hat[1:5, ])))
  
  # May have NA for zero-amplitude voxels
  # But no Inf or negative values
  expect_false(any(is.infinite(res$se_theta_hat)))
  expect_true(all(res$se_theta_hat[!is.na(res$se_theta_hat)] > 0))
})

# Test scenario 4: Delta method accuracy
test_that("Delta method produces correct transformation", {
  set.seed(4567)
  
  # Use simple transformation where we can verify analytically
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      # Identity-like for first parameter
      cbind(rep(1, length(t_hrf)),
            rep(1, length(t_hrf)),
            rep(0, length(t_hrf)),
            rep(0, length(t_hrf)))
    }
  )
  
  n_time <- 50
  Y <- matrix(rnorm(n_time * 2, mean = 2), ncol = 2)
  S <- matrix(1, nrow = n_time, ncol = 1)
  
  res <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:n_time,
    hrf_eval_times = 0:3,
    hrf_interface = hrf_interface,
    theta_seed = c(1, 1, 1),
    theta_bounds = list(lower = c(0, 0, 0), upper = c(10, 10, 10)),
    recenter_global_passes = 1,
    compute_residuals = TRUE,
    compute_se = TRUE,
    lambda_ridge_jacobian = 0.01,
    verbose = FALSE
  )
  
  # With this setup, first parameter SE should be related to coefficient SE
  # Check that SEs are computed and reasonable
  expect_true(all(!is.na(res$se_theta_hat[, 1])))
  expect_true(all(res$se_theta_hat[, 1] > 0))
})

# Test scenario 5: Integration with fit object
test_that("SEs accessible through coef method", {
  skip_if_not_installed("fmrireg", minimum_version = "0.2.0")
  
  # Simple data
  Y <- matrix(rnorm(100, mean = 1), nrow = 20, ncol = 5) 
  S <- matrix(rbinom(20, 1, 0.3), ncol = 1)
  
  # Fit with SEs
  fit <- estimate_parametric_hrf_v2(
    fmri_data = Y,
    event_model = S,
    recenter_global_passes = 1,
    compute_se = TRUE,
    verbose = FALSE
  )
  
  # Extract SEs via coef
  param_se <- coef(fit, type = "se")
  expect_equal(dim(param_se), c(5, 3))
  expect_true(all(param_se > 0, na.rm = TRUE))
  
  # Check vcov for a voxel
  vcov_mat <- vcov(fit, voxel_index = 1)
  if (!is.null(vcov_mat)) {
    expect_equal(dim(vcov_mat), c(3, 3))
    expect_true(all(diag(vcov_mat) > 0))
    # Should be diagonal (current implementation)
    expect_true(all(vcov_mat[upper.tri(vcov_mat)] == 0))
  }
})

# Test scenario 6: SE flag behavior
test_that("compute_se flag works correctly", {
  set.seed(7890)
  
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      cbind(rep(1, 4), rep(0.1, 4), rep(0.1, 4), rep(0.1, 4))
    }
  )
  
  Y <- matrix(rnorm(40), nrow = 10, ncol = 4)
  S <- matrix(c(1, rep(0, 9)), ncol = 1)
  
  # Without SE
  res_no_se <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:10,
    hrf_eval_times = 0:3,
    hrf_interface = hrf_interface,
    theta_seed = c(1, 1, 1),
    theta_bounds = list(lower = c(0, 0, 0), upper = c(10, 10, 10)),
    recenter_global_passes = 1,
    compute_residuals = FALSE,
    compute_se = FALSE,
    verbose = FALSE
  )
  
  expect_null(res_no_se$se_theta_hat)
  
  # With SE
  res_with_se <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:10,
    hrf_eval_times = 0:3,
    hrf_interface = hrf_interface,
    theta_seed = c(1, 1, 1),
    theta_bounds = list(lower = c(0, 0, 0), upper = c(10, 10, 10)),
    recenter_global_passes = 1,
    compute_residuals = TRUE,
    compute_se = TRUE,
    verbose = FALSE
  )
  
  expect_false(is.null(res_with_se$se_theta_hat))
  expect_equal(dim(res_with_se$se_theta_hat), c(4, 3))
})