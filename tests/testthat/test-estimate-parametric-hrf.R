library(testthat)
library(fmriparametric)

context("estimate_parametric_hrf input validation")

# simple placeholder objects
data_obj <- "fmri_data"
event_obj <- "event_model"


test_that("defaults are applied when theta_seed and bounds are NULL", {
  res <- estimate_parametric_hrf(data_obj, event_obj)
  expect_equal(res$theta_seed, c(6, 2.5, 0.35))
  expect_true(is.list(res$theta_bounds))
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
