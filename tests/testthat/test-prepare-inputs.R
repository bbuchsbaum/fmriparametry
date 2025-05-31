library(fmriparametric)
library(testthat)

context("prepare_parametric_inputs")

# create simple BOLD matrix (10 timepoints x 2 voxels)
Y <- matrix(rnorm(20), nrow = 10, ncol = 2)

# simple event design (10x1)
S <- matrix(rbinom(10, 1, 0.2), ncol = 1)

# run without confounds
res <- .prepare_parametric_inputs(Y, S)

test_that("output components have correct dimensions", {
  expect_equal(nrow(res$Y_raw), 10)
  expect_equal(ncol(res$Y_raw), 2)
  expect_equal(dim(res$S_target), c(10, 1))
  expect_equal(length(res$scan_times), 10)
})

# with confound projection
res2 <- .prepare_parametric_inputs(Y, S, confound_formula = ~ poly(scan, 2))

test_that("confound projection reduces mean", {
  expect_true(all(abs(colMeans(res2$Y_proj)) < 1e-8))
  expect_true(all(abs(colMeans(res2$S_target_proj)) < 1e-8))
})
