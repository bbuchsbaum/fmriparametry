library(testthat)

context("R-squared and residuals computation")

# Test scenario 1: R² calculation accuracy
test_that("R-squared calculation is accurate", {
  set.seed(111)
  
  # Create perfect fit scenario
  n_time <- 50
  n_vox <- 10
  
  # Design matrix
  X <- cbind(1, rnorm(n_time), rnorm(n_time))
  
  # True coefficients
  true_beta <- matrix(c(2, 0.5, -0.3), nrow = 3, ncol = n_vox)
  
  # Generate Y with varying noise levels
  Y <- X %*% true_beta
  noise_levels <- seq(0, 0.5, length.out = n_vox)
  for (v in 1:n_vox) {
    Y[, v] <- Y[, v] + rnorm(n_time, sd = noise_levels[v])
  }
  
  # Simple HRF interface that returns X as Taylor basis
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      # Return identity-like basis
      cbind(rep(1, length(t_hrf)), 
            seq_along(t_hrf) / length(t_hrf),
            (seq_along(t_hrf) / length(t_hrf))^2)
    }
  )
  
  # Dummy stimulus
  S <- matrix(1, nrow = n_time, ncol = 1)
  
  res <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:n_time,
    hrf_eval_times = 1:3,
    hrf_interface = hrf_interface,
    theta_seed = c(1, 1, 1),
    theta_bounds = list(lower = c(0, 0, 0), upper = c(10, 10, 10)),
    recenter_global_passes = 1,
    compute_residuals = TRUE,
    verbose = FALSE
  )
  
  # Manual R² calculation
  Y_mean <- colMeans(Y)
  SS_tot <- colSums((Y - matrix(Y_mean, nrow = n_time, ncol = n_vox, byrow = TRUE))^2)
  SS_res <- colSums(res$residuals^2)
  r2_manual <- 1 - SS_res / SS_tot
  
  # Compare
  expect_equal(res$r_squared, r2_manual, tolerance = 1e-6)
  
  # R² should decrease with noise
  expect_true(all(diff(res$r_squared) < 0.1))  # Mostly decreasing
  expect_true(res$r_squared[1] > 0.9)  # Low noise = high R²
  expect_true(res$r_squared[n_vox] < 0.7)  # High noise = lower R²
})

# Test scenario 2: Residuals properties
test_that("residuals have correct properties", {
  set.seed(222)
  
  # Generate simple linear data
  n_time <- 100
  n_vox <- 5
  t <- seq(0, 1, length.out = n_time)
  
  Y <- matrix(0, nrow = n_time, ncol = n_vox)
  for (v in 1:n_vox) {
    Y[, v] <- 2 + 3*t + rnorm(n_time, sd = 0.1)
  }
  
  # Linear HRF interface
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      cbind(rep(1, length(t_hrf)), t_hrf)
    }
  )
  
  S <- matrix(1, nrow = n_time, ncol = 1)
  
  res <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = t,
    hrf_eval_times = t[1:2],
    hrf_interface = hrf_interface,
    theta_seed = c(1, 1),
    theta_bounds = list(lower = c(0, 0), upper = c(10, 10)),
    recenter_global_passes = 1,
    compute_residuals = TRUE,
    verbose = FALSE
  )
  
  # Check residual properties
  expect_equal(dim(res$residuals), c(n_time, n_vox))
  
  # Mean should be near zero
  expect_true(all(abs(colMeans(res$residuals)) < 0.05))
  
  # Should be uncorrelated with fitted values
  fitted <- Y - res$residuals
  cors <- sapply(1:n_vox, function(v) cor(fitted[,v], res$residuals[,v]))
  expect_true(all(abs(cors) < 0.1))
  
  # Variance should match noise level
  expect_true(all(apply(res$residuals, 2, sd) < 0.2))
})

# Test scenario 3: Residuals with poor fit
test_that("residuals capture misspecification", {
  set.seed(333)
  
  # Generate nonlinear data
  n_time <- 50
  t <- seq(0, 2*pi, length.out = n_time)
  Y <- matrix(sin(t), nrow = n_time, ncol = 3)
  Y <- Y + matrix(rnorm(n_time * 3, sd = 0.05), ncol = 3)
  
  # Linear HRF interface (misspecified)
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      cbind(rep(1, length(t_hrf)), t_hrf, t_hrf^2)
    }
  )
  
  S <- matrix(1, nrow = n_time, ncol = 1)
  
  res <- .parametric_engine_iterative(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:n_time,
    hrf_eval_times = 1:3,
    hrf_interface = hrf_interface,
    theta_seed = c(1, 1, 1),
    theta_bounds = list(lower = c(0, 0, 0), upper = c(10, 10, 10)),
    recenter_global_passes = 1,
    compute_residuals = TRUE,
    verbose = FALSE
  )
  
  # R² should be poor
  expect_true(all(res$r_squared < 0.5))
  
  # Residuals should show systematic pattern
  # Check autocorrelation
  acf_vals <- acf(res$residuals[,1], plot = FALSE)$acf[2:5]
  expect_true(any(abs(acf_vals) > 0.3))  # Significant autocorrelation
})

# Test scenario 4: Integration with S3 methods
test_that("residuals method works with fit object", {
  skip_if_not_installed("fmrireg", minimum_version = "0.2.0")
  
  # Simple test data
  Y <- matrix(rnorm(100), nrow = 20, ncol = 5)
  S <- matrix(rbinom(20, 1, 0.3), ncol = 1)
  
  # Create fit object
  fit <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = S,
    global_refinement = TRUE,
    global_passes = 1,
    compute_se = FALSE,
    verbose = FALSE
  )
  
  # Extract residuals
  resid <- residuals(fit)
  expect_equal(dim(resid), dim(Y))
  
  # Check fitted values
  fitted_vals <- fitted(fit, Y_proj = Y)
  expect_equal(dim(fitted_vals), dim(Y))
  
  # Fitted + residuals should equal original
  expect_equal(fitted_vals + resid, Y, tolerance = 1e-10)
  
  # Check R² is stored
  expect_true(!is.null(fit$r_squared))
  expect_length(fit$r_squared, ncol(Y))
  expect_true(all(fit$r_squared >= -1 & fit$r_squared <= 1))
})

# Test scenario 5: Boundary cases
test_that("R² and residuals handle edge cases", {
  set.seed(555)
  
  hrf_interface <- list(
    taylor_basis = function(theta0, t_hrf) {
      cbind(rep(1, length(t_hrf)), rep(0.1, length(t_hrf)))
    }
  )
  
  # Case 1: All zero data
  Y_zero <- matrix(0, nrow = 10, ncol = 3)
  S <- matrix(c(1, rep(0, 9)), ncol = 1)
  
  res_zero <- .parametric_engine_iterative(
    Y_proj = Y_zero,
    S_target_proj = S,
    scan_times = 1:10,
    hrf_eval_times = 0:2,
    hrf_interface = hrf_interface,
    theta_seed = c(1, 1),
    theta_bounds = list(lower = c(0, 0), upper = c(10, 10)),
    recenter_global_passes = 1,
    compute_residuals = TRUE,
    verbose = FALSE
  )
  
  # R² undefined when SS_tot = 0, often set to 0 or -Inf
  expect_true(all(res_zero$r_squared <= 0))
  expect_true(all(abs(res_zero$residuals) < 1e-10))
  
  # Case 2: Constant data
  Y_const <- matrix(5, nrow = 10, ncol = 3)
  
  res_const <- .parametric_engine_iterative(
    Y_proj = Y_const,
    S_target_proj = S,
    scan_times = 1:10,
    hrf_eval_times = 0:2,
    hrf_interface = hrf_interface,
    theta_seed = c(1, 1),
    theta_bounds = list(lower = c(0, 0), upper = c(10, 10)),
    recenter_global_passes = 1,
    compute_residuals = TRUE,
    verbose = FALSE
  )
  
  # Should handle constant data gracefully
  expect_false(any(is.na(res_const$r_squared)))
  expect_false(any(is.na(res_const$residuals)))
})