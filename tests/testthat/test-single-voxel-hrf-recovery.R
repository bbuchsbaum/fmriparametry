library(fmriparametric)

# This test verifies that the estimation procedure can recover the
# shape of the HRF in a simple single-voxel scenario.

test_that("single_voxel_sanity_check recovers HRF shape", {
  set.seed(123)
  n_time <- 120
  onsets <- rep(0, n_time)
  onsets[c(20, 50, 80, 110)] <- 1

  true_theta <- c(tau = 6, sigma = 2.5, rho = 0.4)
  t_hrf <- seq(0, 30, by = 0.5)
  true_hrf <- fmriparametric:::.lwu_hrf_function(t_hrf, true_theta)

  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  signal <- conv_full[1:n_time]
  y <- signal + rnorm(n_time, sd = 0.05)

  fit <- single_voxel_sanity_check(y, onsets, hrf_eval_times = t_hrf)

  est_hrf <- fmriparametric:::.lwu_hrf_function(t_hrf, fit$theta_hat)
  cor_hrf <- cor(est_hrf, true_hrf)

  expect_true(cor_hrf > 0.9)
  expect_true(all(abs(fit$theta_hat - true_theta) < c(0.5, 0.5, 0.2)))
  expect_gt(fit$r_squared, 0.8)
})

test_that("estimate_parametric_hrf recovers single voxel HRF", {
  set.seed(456)
  n_time <- 100
  onsets <- rep(0, n_time)
  onsets[c(15, 45, 75)] <- 1

  true_theta <- c(tau = 5.5, sigma = 2.0, rho = 0.3)
  t_hrf <- seq(0, 30, by = 0.5)
  true_hrf <- fmriparametric:::.lwu_hrf_function(t_hrf, true_theta)

  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  signal <- conv_full[1:n_time]
  Y <- matrix(signal + rnorm(n_time, sd = 0.05), ncol = 1)

  fit <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = matrix(onsets, ncol = 1),
    parametric_hrf = "lwu",
    hrf_eval_times = t_hrf,
    verbose = FALSE
  )

  params <- as.numeric(coef(fit))
  est_hrf <- fmriparametric:::.lwu_hrf_function(t_hrf, params)
  cor_hrf <- cor(est_hrf, true_hrf)

  expect_true(cor_hrf > 0.9)
  expect_true(all(abs(params - true_theta) < c(0.5, 0.5, 0.2)))
  expect_gt(mean(fit$r_squared), 0.8)
})
