library(fmriparametric)
library(testthat)

# basic time vector
x <- seq(0, 10, length.out = 5)

# simple tests for wrapper functions

test_that(".lwu_hrf_function returns correct length", {
  y <- .lwu_hrf_function(x, c(6, 2.5, 0.35))
  expect_type(y, "double")
  expect_length(y, length(x))
})

test_that(".lwu_hrf_taylor_basis_function dimensions", {
  mat <- .lwu_hrf_taylor_basis_function(c(6, 2.5, 0.35), x)
  expect_equal(dim(mat), c(length(x), 4))
})

test_that("parameter names and defaults", {
  expect_equal(.lwu_hrf_parameter_names(), c("tau", "sigma", "rho"))
  expect_equal(.lwu_hrf_default_seed(), c(6, 2.5, 0.35))
  b <- .lwu_hrf_default_bounds()
  expect_named(b, c("lower", "upper"))
  expect_equal(length(b$lower), 3)
  expect_equal(length(b$upper), 3)
})
