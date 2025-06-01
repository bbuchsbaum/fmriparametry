library(testthat)

context("Parallel processing")

# Verify that parallel and sequential fits yield the same results

test_that("Parallel processing matches sequential results", {
  skip_if_not(parallel::detectCores() > 1, "Single core system")
  skip_if_not_installed("future")

  set.seed(101)
  n_time <- 40
  n_vox <- 6
  fmri_data <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  event_model <- matrix(rbinom(n_time, 1, 0.1), ncol = 1)

  fit_seq <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    parallel = FALSE,
    compute_se = FALSE,
    verbose = FALSE
  )

  fit_par <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    parallel = TRUE,
    n_cores = 2,
    compute_se = FALSE,
    verbose = FALSE
  )

  expect_equal(coef(fit_seq), coef(fit_par), tolerance = 1e-6)
  expect_equal(fit_seq$r_squared, fit_par$r_squared, tolerance = 1e-6)
})
