library(testthat)
library(fmriparametric)

test_that("LWU C++ HRF kernel matches explicit formula", {
  t <- seq(0, 30, length.out = 61)
  params <- c(tau = 6.3, sigma = 2.4, rho = 0.41)

  ref <- exp(-((t - params[1])^2) / (2 * params[2]^2)) -
    params[3] * exp(-((t - (params[1] + 2 * params[2]))^2) / (2 * (1.6 * params[2])^2))

  got <- fmriparametric:::.lwu_hrf_formula(
    t = t,
    tau = params[1],
    sigma = params[2],
    rho = params[3],
    normalize = "none"
  )

  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("LWU C++ Taylor basis matches analytic derivatives", {
  t <- seq(0, 30, length.out = 61)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  theta <- c(6.1, 2.2, 0.35)

  theta0 <- fmriparametric:::.apply_bounds(theta, bounds, epsilon = c(1e-6, 1e-6, 1e-6))
  theta0[2] <- max(theta0[2], 0.051)

  ref <- matrix(0, nrow = length(t), ncol = 4)
  tau <- theta0[1]
  sigma <- theta0[2]
  rho <- theta0[3]
  sigma2 <- sigma^2
  sigma3 <- sigma^3
  dt <- t - tau
  u <- t - tau - 2 * sigma
  term1 <- exp(-(dt^2) / (2 * sigma2))
  term2 <- exp(-(u^2) / (5.12 * sigma2))

  ref[, 1] <- term1 - rho * term2
  ref[, 2] <- term1 * (dt / sigma2) - rho * term2 * (u / (2.56 * sigma2))
  ref[, 3] <- term1 * ((dt^2) / sigma3) - rho * term2 * ((u * dt) / (2.56 * sigma3))
  ref[, 4] <- -term2

  got <- fmriparametric:::.lwu_hrf_taylor_basis_function(
    params_vector0 = theta,
    t_hrf_eval = t,
    bounds = bounds
  )

  expect_equal(got, ref, tolerance = 1e-11)

  # Sanity check against central finite differences.
  hrf_fn <- function(p) {
    exp(-((t - p[1])^2) / (2 * p[2]^2)) -
      p[3] * exp(-((t - (p[1] + 2 * p[2]))^2) / (2 * (1.6 * p[2])^2))
  }
  for (k in seq_len(3)) {
    step <- 1e-5
    plus <- theta0
    minus <- theta0
    plus[k] <- plus[k] + step
    minus[k] <- minus[k] - step
    if (k == 2) {
      plus[2] <- max(plus[2], 0.051)
      minus[2] <- max(minus[2], 0.051)
    }
    fd <- (hrf_fn(plus) - hrf_fn(minus)) / (plus[k] - minus[k])
    expect_equal(got[, k + 1], fd, tolerance = 5e-5)
  }
})

test_that("LWU C++ Taylor basis remains finite at bounds", {
  t <- seq(0, 30, length.out = 61)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()

  edge_cases <- list(
    c(bounds$lower[1], bounds$lower[2], bounds$lower[3]),
    c(bounds$upper[1], bounds$upper[2], bounds$upper[3]),
    c(bounds$lower[1], bounds$upper[2], bounds$lower[3])
  )

  for (theta in edge_cases) {
    basis <- fmriparametric:::.lwu_hrf_taylor_basis_function(
      params_vector0 = theta,
      t_hrf_eval = t,
      bounds = bounds
    )
    expect_equal(dim(basis), c(length(t), 4))
    expect_true(all(is.finite(basis)))
  }
})

test_that("LWU GN C++ objective matches R reference", {
  set.seed(1)
  n_time <- 120
  t <- seq(0, 30, length.out = 61)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  theta <- c(6.2, 2.4, 0.32)

  signal <- numeric(n_time)
  signal[seq(10, n_time - 10, by = 20)] <- 1
  h <- fmriparametric:::.lwu_hrf_function(t, theta, bounds)
  x <- fmriparametric:::.convolve_signal_with_kernels(
    signal,
    matrix(h, ncol = 1),
    n_time
  )[, 1]
  y <- 0.9 * x + rnorm(n_time, sd = 0.05)

  denom <- sum(x^2)
  ref <- if (denom < 1e-8) {
    Inf
  } else {
    beta <- as.numeric(crossprod(x, y)) / denom
    sum((y - beta * x)^2)
  }

  got <- fmriparametric:::lwu_gn_objective_cpp(
    theta = theta,
    y = y,
    signal = signal,
    t_hrf_eval = t,
    lower = bounds$lower,
    upper = bounds$upper
  )

  expect_true(is.finite(got))
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("LWU GN C++ jacobian path agrees with R fallback", {
  set.seed(2)
  n_time <- 120
  t <- seq(0, 30, length.out = 61)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  interface <- fmriparametric:::.create_hrf_interface("lwu", user_bounds = bounds)
  theta <- c(6.1, 2.2, 0.38)

  signal <- numeric(n_time)
  signal[seq(8, n_time - 8, by = 16)] <- 1
  h <- interface$hrf_function(t, theta)
  x <- fmriparametric:::.convolve_signal_with_kernels(
    signal,
    matrix(h, ncol = 1),
    n_time
  )[, 1]
  y <- 0.7 * x + rnorm(n_time, sd = 0.04)

  # Explicit R reference implementation for jacobian/residuals.
  basis <- interface$taylor_basis(theta, t)
  X_conv <- fmriparametric:::.convolve_signal_with_kernels(signal, basis, n_time)
  x_hrf <- X_conv[, 1]
  denom <- sum(x_hrf^2)
  beta_ref <- as.numeric(crossprod(x_hrf, y)) / denom
  resid_ref <- y - beta_ref * x_hrf
  jac_ref <- matrix(0, nrow = n_time, ncol = 3)
  for (k in seq_len(3)) {
    dx <- X_conv[, k + 1]
    dbeta <- (as.numeric(crossprod(dx, y)) - 2 * beta_ref * as.numeric(crossprod(x_hrf, dx))) / denom
    jac_ref[, k] <- -beta_ref * dx - dbeta * x_hrf
  }

  # Explicit C++ call
  cpp <- fmriparametric:::lwu_gn_jacobian_cpp(
    theta = theta,
    y = y,
    signal = signal,
    t_hrf_eval = t,
    lower = bounds$lower,
    upper = bounds$upper
  )

  expect_true(isTRUE(cpp$ok))
  expect_equal(cpp$amplitude, beta_ref, tolerance = 1e-10)
  expect_equal(cpp$residuals, resid_ref, tolerance = 1e-10)
  expect_equal(cpp$jacobian, jac_ref, tolerance = 1e-10)
})
