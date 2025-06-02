library(testthat)
library(fmriparametric)

context("batch convolution helper")

time_len <- 20
signals <- cbind(rep(c(1,0,0), length.out = time_len),
                 rep(c(0,1,0), length.out = time_len))
kernels <- matrix(rnorm(5*2), nrow = 5)

result <- fmriparametric:::.batch_convolution(signals, kernels, time_len)
manual <- fmriparametric:::.fast_batch_convolution(signals[,1], kernels, time_len) +
          fmriparametric:::.fast_batch_convolution(signals[,2], kernels, time_len)

test_that("batch convolution sums across signals", {
  expect_equal(result, manual, tolerance = 1e-10)
})
