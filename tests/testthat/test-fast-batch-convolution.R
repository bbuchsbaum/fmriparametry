context("fast batch convolution")

library(fmriparametric)

# simple reproducible test for .fast_batch_convolution

test_that("fast batch convolution matches direct convolution", {
  set.seed(42)
  signal <- rnorm(40)
  kernels <- matrix(rnorm(61 * 3), nrow = 61, ncol = 3)

  # direct convolution using base R convolve
  direct <- matrix(0, nrow = length(signal), ncol = 3)
  for (j in seq_len(ncol(kernels))) {
    conv_result <- convolve(signal, rev(kernels[, j]), type = "open")
    direct[, j] <- conv_result[seq_len(length(signal))]
  }

  # our helper
  fast <- .fast_batch_convolution(signal, kernels, length(signal))

  expect_equal(fast, direct, tolerance = 1e-8)
})


test_that("cached QR solve produces same result", {
  set.seed(123)
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  Y <- matrix(rnorm(40), nrow = 10, ncol = 4)

  .clear_qr_cache()
  res1 <- .cached_qr_solve(X, Y, "test")
  res2 <- .cached_qr_solve(X, Y, "test")

  expect_equal(res1, res2)

  # result should match base qr.solve
  qr_obj <- qr(X)
  base_res <- qr.solve(qr_obj, Y)
  expect_equal(res1, base_res)
})
