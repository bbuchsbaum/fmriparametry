library(testthat)

context("parametric engine")

# Test Scenario 1: Basic functionality with dummy interface
test_that("engine returns correct structure with basic inputs", {
  # Dummy interface that returns a simple Taylor basis
  hrf_iface <- list(
    taylor_basis = function(theta0, t_hrf) {
      matrix(rep(c(1, 0.1, 0.2, 0.05), each = length(t_hrf)), nrow = length(t_hrf), byrow = FALSE)
    }
  )
  
  Y <- matrix(rnorm(20), nrow = 10, ncol = 2)
  S <- matrix(rbinom(10, 1, 0.2), ncol = 1)
  scan_t <- seq_len(10)
  t_hrf <- c(0, 1)
  
  res <- .parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = scan_t,
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_iface,
    theta_seed = c(1, 1, 1),
    theta_bounds = list(lower = c(0, 0, 0), upper = c(2, 2, 2))
  )
  
  expect_type(res, "list")
  expect_true(all(c("theta_hat", "beta0") %in% names(res)))
  expect_equal(nrow(res$theta_hat), ncol(Y))
  expect_equal(ncol(res$theta_hat), 3)
  expect_length(res$beta0, ncol(Y))
})

# Test Scenario 2: Perfect recovery with known parameters
test_that("engine recovers known parameters when seeded at truth", {
  set.seed(123)
  n_time <- 100
  n_vox <- 5
  
  # True parameters for LWU model
  true_params <- c(tau = 6, sigma = 2.5, rho = 0.35)
  true_amp <- 2.0
  
  # Create HRF interface that returns identity at true params
  hrf_iface <- list(
    taylor_basis = function(theta0, t_hrf) {
      # At true params, derivatives should be zero
      if (all(abs(theta0 - true_params) < 1e-10)) {
        # HRF value = 1, all derivatives = 0
        matrix(c(rep(1, length(t_hrf)), rep(0, length(t_hrf) * 3)), 
               nrow = length(t_hrf), ncol = 4)
      } else {
        # Simple linear approximation
        hrf_val <- rep(1, length(t_hrf))
        deriv1 <- rep(0.1, length(t_hrf))
        deriv2 <- rep(0.2, length(t_hrf))
        deriv3 <- rep(0.05, length(t_hrf))
        cbind(hrf_val, deriv1, deriv2, deriv3)
      }
    }
  )
  
  # Generate stimulus
  S <- matrix(0, nrow = n_time, ncol = 1)
  S[seq(10, n_time, by = 20), 1] <- 1
  
  # Generate data with known convolution
  X_true <- matrix(0, nrow = n_time, ncol = 1)
  hrf_true <- rep(1, 10)  # Simple HRF shape
  for (i in which(S[,1] == 1)) {
    idx <- i:(min(i + length(hrf_true) - 1, n_time))
    X_true[idx, 1] <- X_true[idx, 1] + hrf_true[1:length(idx)]
  }
  
  Y <- matrix(rep(true_amp * X_true[,1], n_vox), ncol = n_vox)
  Y <- Y + matrix(rnorm(n_time * n_vox, sd = 0.01), nrow = n_time)  # Small noise
  
  res <- .parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = seq_len(n_time),
    hrf_eval_times = 0:9,
    hrf_interface = hrf_iface,
    theta_seed = true_params,
    theta_bounds = list(lower = c(0, 0.1, 0), upper = c(15, 10, 1.5))
  )
  
  # Should recover parameters very close to seed (since we're at truth)
  expect_true(all(abs(res$theta_hat - matrix(true_params, nrow = n_vox, ncol = 3, byrow = TRUE)) < 0.1))
  expect_true(all(abs(res$beta0 - true_amp) < 0.5))
})

# Test Scenario 3: Convergence from offset seed
test_that("engine improves parameters when starting from offset seed", {
  set.seed(456)
  n_time <- 50
  
  # Create interface with non-zero derivatives
  hrf_iface <- list(
    taylor_basis = function(theta0, t_hrf) {
      # Simulate realistic Taylor expansion
      hrf_val <- exp(-(t_hrf - theta0[1])^2 / (2 * theta0[2]^2))
      # Approximate derivatives
      d_tau <- 0.1 * (t_hrf - theta0[1]) * hrf_val / theta0[2]^2
      d_sigma <- 0.1 * (t_hrf - theta0[1])^2 * hrf_val / theta0[2]^3
      d_rho <- -0.05 * exp(-(t_hrf - theta0[1] - 2*theta0[2])^2 / (2 * (1.6*theta0[2])^2))
      cbind(hrf_val, d_tau, d_sigma, d_rho)
    }
  )
  
  # True parameters
  true_params <- c(tau = 5, sigma = 2, rho = 0.3)
  offset_seed <- c(tau = 7, sigma = 3, rho = 0.5)  # Start away from truth
  
  # Generate stimulus and data
  S <- matrix(0, nrow = n_time, ncol = 1)
  S[c(5, 20, 35), 1] <- 1
  
  # Use true HRF to generate data
  t_hrf <- 0:15
  true_hrf <- exp(-(t_hrf - true_params[1])^2 / (2 * true_params[2]^2)) - 
              true_params[3] * exp(-(t_hrf - true_params[1] - 2*true_params[2])^2 / (2 * (1.6*true_params[2])^2))
  
  Y <- stats::filter(S[,1], true_hrf, sides = 1)
  Y[is.na(Y)] <- 0
  Y <- matrix(Y + rnorm(n_time, sd = 0.1), ncol = 1)
  
  res <- .parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = seq_len(n_time),
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_iface,
    theta_seed = offset_seed,
    theta_bounds = list(lower = c(0, 0.5, 0), upper = c(15, 10, 1.5))
  )
  
  # Parameters should move toward truth
  initial_error <- sum((offset_seed - true_params)^2)
  final_error <- sum((res$theta_hat[1,] - true_params)^2)
  expect_true(final_error < initial_error)
})

# Test Scenario 4: Parameter bounds enforcement
test_that("engine respects parameter bounds", {
  hrf_iface <- list(
    taylor_basis = function(theta0, t_hrf) {
      # Interface that would push parameters outside bounds
      cbind(rep(1, length(t_hrf)), 
            rep(-10, length(t_hrf)),  # Strong negative gradient for tau
            rep(10, length(t_hrf)),   # Strong positive gradient for sigma  
            rep(5, length(t_hrf)))    # Strong positive gradient for rho
    }
  )
  
  Y <- matrix(rnorm(50, mean = 10), nrow = 10, ncol = 5)
  S <- matrix(c(rep(0, 4), 1, rep(0, 5)), ncol = 1)
  
  bounds <- list(lower = c(2, 1, 0.1), upper = c(8, 4, 0.8))
  
  res <- .parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:10,
    hrf_eval_times = 0:5,
    hrf_interface = hrf_iface,
    theta_seed = c(5, 2.5, 0.4),
    theta_bounds = bounds
  )
  
  # All parameters should be within bounds
  expect_true(all(res$theta_hat[,1] >= bounds$lower[1]))
  expect_true(all(res$theta_hat[,1] <= bounds$upper[1]))
  expect_true(all(res$theta_hat[,2] >= bounds$lower[2]))
  expect_true(all(res$theta_hat[,2] <= bounds$upper[2]))
  expect_true(all(res$theta_hat[,3] >= bounds$lower[3]))
  expect_true(all(res$theta_hat[,3] <= bounds$upper[3]))
})

# Test Scenario 5: Numerical stability
test_that("engine handles poorly conditioned data gracefully", {
  hrf_iface <- list(
    taylor_basis = function(theta0, t_hrf) {
      # Nearly collinear basis
      base <- rep(1, length(t_hrf))
      cbind(base, base * 1.001, base * 0.999, base * 1.002)
    }
  )
  
  # Very small signal
  Y <- matrix(rnorm(30, sd = 1e-6), nrow = 10, ncol = 3)
  S <- matrix(rbinom(10, 1, 0.3), ncol = 1)
  
  expect_no_error({
    res <- .parametric_engine(
      Y_proj = Y,
      S_target_proj = S,
      scan_times = 1:10,
      hrf_eval_times = 0:3,
      hrf_interface = hrf_iface,
      theta_seed = c(5, 2, 0.3),
      theta_bounds = list(lower = c(0, 0.1, 0), upper = c(10, 5, 1)),
      lambda_ridge = 0.1  # Higher ridge for stability
    )
  })
  
  # No NaN or Inf values
  expect_false(any(is.na(res$theta_hat)))
  expect_false(any(is.infinite(res$theta_hat)))
  expect_false(any(is.na(res$beta0)))
  expect_false(any(is.infinite(res$beta0)))
})

# Test Scenario 6: Near-zero amplitude handling
test_that("engine handles near-zero amplitudes without NaNs", {
  hrf_iface <- list(
    taylor_basis = function(theta0, t_hrf) {
      matrix(c(rep(1, length(t_hrf)), 
               rep(0.1, length(t_hrf)), 
               rep(0.2, length(t_hrf)), 
               rep(0.05, length(t_hrf))), 
             nrow = length(t_hrf), ncol = 4)
    }
  )
  
  # Multiple voxels with varying amplitudes including zero
  Y <- cbind(
    rep(0, 10),                    # Zero signal
    rnorm(10, mean = 0, sd = 1e-8), # Near-zero signal
    rnorm(10, mean = 1, sd = 0.1)   # Normal signal
  )
  S <- matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0, 0), ncol = 1)
  
  res <- .parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:10,
    hrf_eval_times = 0:4,
    hrf_interface = hrf_iface,
    theta_seed = c(5, 2.5, 0.35),
    theta_bounds = list(lower = c(0, 0.5, 0), upper = c(10, 5, 1)),
    epsilon_beta = 1e-6
  )
  
  expect_false(any(is.na(res$theta_hat)))
  expect_false(any(is.infinite(res$theta_hat)))
  expect_false(any(is.na(res$beta0)))
  expect_false(any(is.infinite(res$beta0)))
  
  # Parameters should still be reasonable even for zero-amplitude voxels
  expect_true(all(res$theta_hat >= 0))
})

# Test Scenario 7: Different ridge penalties
test_that("engine behavior changes appropriately with ridge penalty", {
  hrf_iface <- list(
    taylor_basis = function(theta0, t_hrf) {
      cbind(rep(1, length(t_hrf)),
            seq_along(t_hrf) * 0.1,
            seq_along(t_hrf)^2 * 0.01,
            seq_along(t_hrf)^0.5 * 0.2)
    }
  )
  
  Y <- matrix(rnorm(50, mean = 2), nrow = 10, ncol = 5)
  S <- matrix(c(1, 0, 0, 1, 0, 0, 1, 0, 0, 0), ncol = 1)
  
  # Low ridge penalty
  res_low <- .parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:10,
    hrf_eval_times = 0:4,
    hrf_interface = hrf_iface,
    theta_seed = c(5, 2.5, 0.35),
    theta_bounds = list(lower = c(0, 0.5, 0), upper = c(10, 5, 1)),
    lambda_ridge = 0.001
  )
  
  # High ridge penalty
  res_high <- .parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = 1:10,
    hrf_eval_times = 0:4,
    hrf_interface = hrf_iface,
    theta_seed = c(5, 2.5, 0.35),
    theta_bounds = list(lower = c(0, 0.5, 0), upper = c(10, 5, 1)),
    lambda_ridge = 1.0
  )
  
  # Higher ridge should lead to smaller parameter updates
  low_updates <- rowMeans(abs(res_low$theta_hat - matrix(c(5, 2.5, 0.35), nrow = 5, ncol = 3, byrow = TRUE)))
  high_updates <- rowMeans(abs(res_high$theta_hat - matrix(c(5, 2.5, 0.35), nrow = 5, ncol = 3, byrow = TRUE)))
  
  expect_true(mean(high_updates) < mean(low_updates))
})

# Test Scenario 8: Large-scale performance
test_that("engine performs efficiently with many voxels", {
  skip_if_not(Sys.getenv("FMRIPARAMETRIC_EXTENDED_TESTS") == "true",
              "Skipping extended performance test")
  
  hrf_iface <- list(
    taylor_basis = function(theta0, t_hrf) {
      hrf_val <- exp(-(t_hrf - theta0[1])^2 / (2 * theta0[2]^2))
      cbind(hrf_val, 
            0.1 * hrf_val,
            0.2 * hrf_val, 
            0.05 * hrf_val)
    }
  )
  
  n_vox <- 10000
  n_time <- 200
  
  Y <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  S <- matrix(0, nrow = n_time, ncol = 1)
  S[seq(10, n_time, by = 30), 1] <- 1
  
  start_time <- Sys.time()
  res <- .parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = seq_len(n_time),
    hrf_eval_times = 0:20,
    hrf_interface = hrf_iface,
    theta_seed = c(6, 2.5, 0.35),
    theta_bounds = list(lower = c(0, 0.5, 0), upper = c(15, 10, 1.5))
  )
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
  
  expect_equal(nrow(res$theta_hat), n_vox)
  expect_true(elapsed < 10)  # Should complete in under 10 seconds
})
