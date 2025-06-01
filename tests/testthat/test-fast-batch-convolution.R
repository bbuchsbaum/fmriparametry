library(testthat)

context("fast batch convolution")

test_that("FFT helper matches direct convolution", {
  signal <- rep(c(1, 0, 0), length.out = 50)
  kernels <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)
  n_out <- length(signal)

  result_fft <- .fast_batch_convolution(signal, kernels, n_out)

  direct <- sapply(seq_len(ncol(kernels)), function(j) {
    conv <- convolve(signal, rev(kernels[, j]), type = "open")
    conv[seq_len(n_out)]
  })

  expect_equal(result_fft, direct, tolerance = 1e-10)
})

