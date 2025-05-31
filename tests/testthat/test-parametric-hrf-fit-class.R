library(fmriparametric)
library(testthat)

context("parametric_hrf_fit class")

pars <- matrix(rep(1:3, each = 2), nrow = 2, byrow = TRUE)
amps <- c(1, 2)
obj <- new_parametric_hrf_fit(
  estimated_parameters = pars,
  amplitudes = amps,
  parameter_names = c("tau", "sigma", "rho"),
  metadata = list(n_timepoints = 10)
)

test_that("constructor returns correct class and dimensions", {
  expect_s3_class(obj, "parametric_hrf_fit")
  expect_equal(nrow(obj$estimated_parameters), 2)
  expect_equal(ncol(obj$estimated_parameters), 3)
  expect_equal(length(obj$amplitudes), 2)
})

test_that("helper functions return metadata", {
  expect_equal(n_voxels(obj), 2)
  expect_equal(n_timepoints(obj), 10)
})
