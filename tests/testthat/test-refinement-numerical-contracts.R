library(testthat)
library(fmriparametric)

test_that("multi-regressor stimulus collapsing preserves linear convolution", {
  set.seed(101)
  n_time <- 80
  signal_mat <- matrix(rnorm(n_time * 3), nrow = n_time, ncol = 3)
  hrf <- exp(-((seq(0, 30, length.out = 61) - 6)^2) / (2 * 2.5^2))

  x_full <- fmriparametric:::.convolve_signal_with_kernels(
    signal = signal_mat,
    kernels = matrix(hrf, ncol = 1),
    output_length = n_time
  )[, 1]
  x_summed <- fmriparametric:::.convolve_signal_with_kernels(
    signal = rowSums(signal_mat),
    kernels = matrix(hrf, ncol = 1),
    output_length = n_time
  )[, 1]

  expect_equal(x_full, x_summed, tolerance = 1e-10)
  expect_equal(
    fmriparametric:::.extract_primary_stimulus(signal_mat),
    rowSums(signal_mat),
    tolerance = 0
  )
})

test_that("Gauss-Newton objective intercept path matches closed-form OLS objective", {
  set.seed(202)
  n_time <- 120
  S <- matrix(0, nrow = n_time, ncol = 2)
  S[seq(8, 96, by = 16), 1] <- 1
  S[seq(12, 100, by = 20), 2] <- 0.8
  theta <- c(6.2, 2.1, 0.35)
  t_hrf <- seq(0, 30, length.out = 61)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")

  h <- hrf_interface$hrf_function(t_hrf, theta)
  x <- fmriparametric:::.convolve_signal_with_kernels(
    signal = S,
    kernels = matrix(h, ncol = 1),
    output_length = n_time
  )[, 1]
  y <- 0.7 + 1.25 * x + rnorm(n_time, sd = 0.03)

  obj <- fmriparametric:::.calculate_objective_gn(
    theta = theta,
    y = y,
    S = S,
    t_hrf = t_hrf,
    hrf_interface = hrf_interface,
    n_time = n_time,
    has_intercept = TRUE
  )

  x_centered <- x - mean(x)
  y_centered <- y - mean(y)
  beta <- sum(x_centered * y_centered) / sum(x_centered^2)
  alpha <- mean(y) - beta * mean(x)
  ref_obj <- sum((y - alpha - beta * x)^2)

  expect_true(is.finite(obj))
  expect_equal(obj, ref_obj, tolerance = 1e-8)

  jac <- fmriparametric:::.get_jacobian_and_residuals(
    theta = theta,
    y = y,
    S = S,
    t_hrf = t_hrf,
    hrf_interface = hrf_interface,
    n_time = n_time,
    has_intercept = TRUE
  )
  expect_false(is.null(jac))
  expect_true(all(is.finite(jac$jacobian)))
  expect_equal(sum(jac$residuals^2), obj, tolerance = 1e-8)
})

test_that("apply_bounds epsilon handling cannot invert narrow bounds", {
  bounds <- list(
    lower = c(1.0, 2.0, 3.0),
    upper = c(1.000001, 2.000001, 3.000001)
  )
  params <- c(0, 10, 2)

  bounded <- fmriparametric:::.apply_bounds(params, bounds, epsilon = 1.0)
  expect_true(all(is.finite(bounded)))
  expect_true(all(bounded >= bounds$lower))
  expect_true(all(bounded <= bounds$upper))
})

test_that("theta grouping supports parameter vectors beyond LWU length", {
  theta_current <- cbind(
    c(1, 1, 2, 2),
    c(0.5, 0.5, 0.7, 0.7),
    c(0.2, 0.2, 0.3, 0.3),
    c(4, 4, 5, 5)
  )

  groups <- fmriparametric:::.compute_theta_grid_groups(theta_current = theta_current)
  expect_true(is.list(groups))
  expect_true(length(groups) >= 1)
  expect_true(all(vapply(groups, function(g) length(g$centroid), integer(1)) == ncol(theta_current)))
})

test_that("parametric engine is invariant to voxel column permutation", {
  set.seed(303)
  n_time <- 90
  n_vox <- 5
  S <- matrix(0, nrow = n_time, ncol = 2)
  S[seq(10, 70, by = 15), 1] <- 1
  S[seq(20, 80, by = 20), 2] <- 1
  t_hrf <- seq(0, 30, length.out = 61)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  theta_seed <- hrf_interface$default_seed()
  bounds <- hrf_interface$default_bounds()
  h <- hrf_interface$hrf_function(t_hrf, theta_seed)
  x <- fmriparametric:::.convolve_signal_with_kernels(S, matrix(h, ncol = 1), n_time)[, 1]

  amps <- c(0.8, 1.1, -0.5, 0.4, 1.7)
  ints <- c(0.2, -0.1, 0.3, 0.0, -0.2)
  Y <- vapply(
    seq_len(n_vox),
    function(v) ints[v] + amps[v] * x + rnorm(n_time, sd = 0.02),
    numeric(n_time)
  )
  Y <- as.matrix(Y)

  fit_a <- fmriparametric:::.parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_interface,
    theta_seed = theta_seed,
    theta_bounds = bounds,
    lambda_ridge = 0.01,
    baseline_model = "intercept"
  )

  perm <- sample.int(n_vox)
  fit_b <- fmriparametric:::.parametric_engine(
    Y_proj = Y[, perm, drop = FALSE],
    S_target_proj = S,
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_interface,
    theta_seed = theta_seed,
    theta_bounds = bounds,
    lambda_ridge = 0.01,
    baseline_model = "intercept"
  )

  inv_perm <- order(perm)
  expect_equal(fit_b$theta_hat[inv_perm, , drop = FALSE], fit_a$theta_hat, tolerance = 1e-8)
  expect_equal(fit_b$beta0[inv_perm], fit_a$beta0, tolerance = 1e-8)
  expect_equal(fit_b$r_squared[inv_perm], fit_a$r_squared, tolerance = 1e-10)
})

test_that("parametric engine is scale-invariant for theta and R^2", {
  set.seed(404)
  n_time <- 120
  S <- matrix(0, nrow = n_time, ncol = 2)
  S[seq(8, 96, by = 16), 1] <- 1
  S[seq(14, 110, by = 20), 2] <- 0.8
  t_hrf <- seq(0, 30, length.out = 61)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  theta_seed <- hrf_interface$default_seed()
  bounds <- hrf_interface$default_bounds()

  h <- hrf_interface$hrf_function(t_hrf, theta_seed)
  x <- fmriparametric:::.convolve_signal_with_kernels(S, matrix(h, ncol = 1), n_time)[, 1]
  Y <- cbind(
    0.3 + 1.1 * x + rnorm(n_time, sd = 0.03),
    -0.2 + 0.7 * x + rnorm(n_time, sd = 0.03),
    0.0 - 0.4 * x + rnorm(n_time, sd = 0.03)
  )

  fit_a <- fmriparametric:::.parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_interface,
    theta_seed = theta_seed,
    theta_bounds = bounds,
    lambda_ridge = 0.01,
    baseline_model = "intercept"
  )

  scale_fac <- 3.75
  fit_b <- fmriparametric:::.parametric_engine(
    Y_proj = scale_fac * Y,
    S_target_proj = S,
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_interface,
    theta_seed = theta_seed,
    theta_bounds = bounds,
    lambda_ridge = 0.01,
    baseline_model = "intercept"
  )

  expect_equal(fit_b$theta_hat, fit_a$theta_hat, tolerance = 1e-8)
  expect_equal(fit_b$r_squared, fit_a$r_squared, tolerance = 1e-10)
  expect_equal(fit_b$r_squared_raw, fit_a$r_squared_raw, tolerance = 1e-10)
  expect_equal(fit_b$beta0, scale_fac * fit_a$beta0, tolerance = 1e-7)
})

test_that("fast convolution remains linear below and above FFT threshold", {
  check_linearity <- function(n_time, kernel_len = 31L) {
    set.seed(500 + n_time)
    s1 <- rnorm(n_time)
    s2 <- rnorm(n_time)
    k1 <- dnorm(seq(-3, 3, length.out = kernel_len))
    k2 <- sin(seq(0, pi, length.out = kernel_len))
    kernels <- cbind(k1, k2)

    lhs <- fmriparametric:::.fast_batch_convolution(
      signal = s1 + s2,
      kernels = kernels,
      output_length = n_time
    )
    rhs <- fmriparametric:::.fast_batch_convolution(
      signal = s1,
      kernels = kernels,
      output_length = n_time
    ) + fmriparametric:::.fast_batch_convolution(
      signal = s2,
      kernels = kernels,
      output_length = n_time
    )
    expect_equal(lhs, rhs, tolerance = 1e-8)
  }

  # Below threshold: direct/C++ path
  check_linearity(n_time = 180L, kernel_len = 21L)
  # Above threshold: FFT path
  check_linearity(n_time = 260L, kernel_len = 31L)
})
