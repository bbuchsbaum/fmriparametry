library(testthat)
library(fmriparametric)

context("estimate_parametric_hrf integration")

# simple numeric inputs
data_obj <- matrix(rnorm(20), nrow = 10, ncol = 2)
event_obj <- matrix(rbinom(10, 1, 0.2), ncol = 1)


test_that("defaults are applied when theta_seed and bounds are NULL", {
  res <- estimate_parametric_hrf(data_obj, event_obj)
  expect_s3_class(res, "parametric_hrf_fit")
  expect_equal(res$metadata$theta_seed, c(6, 2.5, 0.35))
  expect_true(is.list(res$metadata$theta_bounds))
  expect_equal(nrow(res$estimated_parameters), ncol(data_obj))
})

test_that("unsupported hrf model errors", {
  expect_error(estimate_parametric_hrf(data_obj, event_obj, parametric_hrf = "bad"))
})

test_that("theta_seed length validated", {
  expect_error(estimate_parametric_hrf(data_obj, event_obj, theta_seed = c(1,2)))
})

test_that("theta_bounds structure validated", {
  expect_error(estimate_parametric_hrf(data_obj, event_obj, theta_bounds = list(lower = c(1,2,3))))
})

test_that("verbose must be logical", {
  expect_error(estimate_parametric_hrf(data_obj, event_obj, verbose = "yes"))
})
