library(testthat)
library(fmriparametric)

context("estimate_parametric_hrf_ultimate wrapper")

test_that("ultimate estimator runs with basic numeric inputs", {
  set.seed(99)
  n_time <- 40
  n_vox <- 4

  fmri_data <- matrix(rnorm(n_time * n_vox, sd = 2), nrow = n_time, ncol = n_vox)
  event_model <- matrix(0, nrow = n_time, ncol = 1)
  event_model[c(8, 20, 32), 1] <- 1

  fit <- fmriparametric:::estimate_parametric_hrf_ultimate(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    iterative_recentering = FALSE,
    verbose = FALSE,
    return_diagnostics = FALSE
  )

  expect_s3_class(fit, "parametric_hrf_fit")
  expect_equal(nrow(fit$estimated_parameters), n_vox)
  expect_equal(ncol(fit$estimated_parameters), 3)
})
