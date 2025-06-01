library(testthat)
library(fmriparametric)

context("refinement wrapper functions")

set.seed(1)
n_time <- 40
n_vox <- 3
fmri <- matrix(rnorm(n_time * n_vox), nrow = n_time)
event <- matrix(0, nrow = n_time, ncol = 1)
event[c(10,20,30)] <- 1

hrf_int <- .create_hrf_interface("lwu")
base <- .parametric_engine(fmri, event, seq_len(n_time),
                           seq(0,30,length.out=61),
                           hrf_int, hrf_int$default_seed(),
                           hrf_int$default_bounds())

mod_res <- .refine_moderate_voxels(1:n_vox, fmri, event, base$theta_hat,
                                   hrf_int, seq(0,30,length.out=61))

hard_res <- .refine_hard_voxels(1:n_vox, fmri, event, base$theta_hat,
                                 hrf_int, seq(0,30,length.out=61), max_iter = 2)

test_that("moderate refinement wrapper returns correct shapes", {
  expect_equal(dim(mod_res$theta_refined), dim(base$theta_hat))
  expect_length(mod_res$amplitudes, n_vox)
})

test_that("hard refinement wrapper returns correct shapes", {
  expect_equal(dim(hard_res$theta_refined), dim(base$theta_hat))
  expect_length(hard_res$amplitudes, n_vox)
})
