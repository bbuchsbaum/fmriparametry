test_that(".smart_convolution handles long output_length without NA", {
  signal <- c(1, 2, 3)
  kernels <- matrix(c(1, 2), ncol = 1)
  out_len <- 6
  result <- .smart_convolution(signal, kernels, out_len)
  expect_equal(nrow(result), out_len)
  expect_false(anyNA(result))
})
