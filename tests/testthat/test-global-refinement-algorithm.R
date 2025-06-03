library(testthat)
library(fmriparametric)

# Test that global iterative refinement improves fit for a simple dataset

test_that("global refinement improves parameter estimates", {
  set.seed(42)
  n_time <- 60
  onsets <- rep(0, n_time)
  onsets[c(10, 30, 50)] <- 1

  true_theta <- c(tau = 6, sigma = 2.5, rho = 0.4)
  t_hrf <- seq(0, 30, by = 0.5)
  true_hrf <- fmriparametric:::.lwu_hrf_function(t_hrf, true_theta)

  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  signal <- conv_full[seq_len(n_time)]
  Y <- matrix(signal + rnorm(n_time, sd = 0.05), ncol = 1)

  event_mat <- matrix(onsets, ncol = 1)

  bad_seed <- c(9, 4, 0.1)

  fit_no_refine <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = event_mat,
    parametric_hrf = "lwu",
    theta_seed = bad_seed,
    global_refinement = FALSE,
    verbose = FALSE
  )

  fit_refine <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = event_mat,
    parametric_hrf = "lwu",
    theta_seed = bad_seed,
    global_refinement = TRUE,
    global_passes = 2,
    verbose = FALSE
  )

  expect_gt(mean(fit_refine$r_squared), mean(fit_no_refine$r_squared))

  dist_no_ref <- sum((coef(fit_no_refine)[1, ] - true_theta)^2)
  dist_refine <- sum((coef(fit_refine)[1, ] - true_theta)^2)
  expect_lt(dist_refine, dist_no_ref)
})

