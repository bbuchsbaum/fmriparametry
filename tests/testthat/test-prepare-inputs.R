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

# event_model supplied as a list
S_list <- list(design_matrix = S)
res_list <- .prepare_parametric_inputs(Y, S_list)

test_that("event_model list input works", {
  expect_equal(res_list$S_target, S)
})

# dataset-style input with masking
dataset <- list(data = Y, scan_times = seq(0, 9))
res_mask <- .prepare_parametric_inputs(dataset, S, mask = c(TRUE, FALSE))

test_that("list dataset input and mask subset", {
  expect_equal(ncol(res_mask$Y_raw), 1)
  expect_equal(res_mask$scan_times, seq(0, 9))
})

# projected data orthogonal to confounds
Z <- model.matrix(~ poly(scan, 2), data = data.frame(scan = 1:10))
res_proj <- .prepare_parametric_inputs(Y, S, confound_formula = ~ poly(scan, 2))

test_that("projected components orthogonal to confounds", {
  expect_true(max(abs(c(crossprod(Z, res_proj$Y_proj)))) < 1e-8)
  expect_true(max(abs(c(crossprod(Z, res_proj$S_target_proj)))) < 1e-8)
})

# defaults for hrf_eval_times use scan spacing
res_default <- .prepare_parametric_inputs(Y, S, hrf_span = 4)

test_that("hrf_eval_times defaults use dt", {
  expect_equal(res_default$hrf_eval_times, seq(0, 4, by = 1))
})

# explicit hrf_eval_times respected
custom_times <- seq(0, 5, by = 0.5)
res_custom <- .prepare_parametric_inputs(Y, S, hrf_eval_times = custom_times)

test_that("hrf_eval_times argument respected", {
  expect_equal(res_custom$hrf_eval_times, custom_times)
})

# zero-event design handled
S0 <- matrix(0, nrow = 10, ncol = 1)
res_zero <- .prepare_parametric_inputs(Y, S0)

test_that("zero events handled", {
  expect_true(all(res_zero$S_target == 0))
  expect_true(all(res_zero$S_target_proj == 0))
})

# wrong event row count errors
S_bad <- matrix(1, nrow = 5, ncol = 1)

test_that("event_model row mismatch errors", {
  expect_error(.prepare_parametric_inputs(Y, S_bad), "wrong number of rows")
})
