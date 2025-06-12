# Edge case tests for convolution functions
library(fmriparametric)
library(testthat)

test_that(".fast_batch_convolution handles edge cases", {
  # Case 1: Empty signal
  kernels <- matrix(rnorm(10), ncol = 1)
  signal_empty <- numeric(0)
  
  result_empty <- fmriparametric:::.fast_batch_convolution(
    signal_empty, kernels, 10
  )
  expect_equal(dim(result_empty), c(10, 1))
  expect_true(all(result_empty == 0))
  
  # Case 2: Single time point
  signal_single <- 1
  result_single <- fmriparametric:::.fast_batch_convolution(
    signal_single, kernels, 1
  )
  expect_equal(dim(result_single), c(1, 1))
  expect_equal(result_single[1, 1], kernels[1, 1])
  
  # Case 3: Very long signals (test FFT path)
  n_long <- 1000
  signal_long <- c(rep(0, 400), rep(1, 10), rep(0, 590))
  kernels_long <- matrix(exp(-seq(0, 5, length.out = 50)), ncol = 1)
  
  result_long <- fmriparametric:::.fast_batch_convolution(
    signal_long, kernels_long, n_long
  )
  
  expect_equal(dim(result_long), c(n_long, 1))
  expect_true(all(is.finite(result_long)))
  
  # Case 4: Multiple kernels
  kernels_multi <- matrix(rnorm(30), ncol = 3)
  signal <- c(0, 0, 1, 0, 0)
  
  result_multi <- fmriparametric:::.fast_batch_convolution(
    signal, kernels_multi, 5
  )
  
  expect_equal(dim(result_multi), c(5, 3))
  expect_true(all(is.finite(result_multi)))
})

test_that(".batch_convolution handles edge cases", {
  # Case 1: Zero signals
  signals_zero <- matrix(0, nrow = 20, ncol = 1)
  kernels <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)
  
  result_zero <- fmriparametric:::.batch_convolution(
    signals_zero, kernels, 20
  )
  
  expect_equal(dim(result_zero), c(20, 3))
  expect_true(all(result_zero == 0))
  
  # Case 2: Delta function stimulus
  signals_delta <- matrix(0, nrow = 50, ncol = 1)
  signals_delta[25, 1] <- 1
  
  result_delta <- fmriparametric:::.batch_convolution(
    signals_delta, kernels, 50
  )
  
  expect_equal(dim(result_delta), c(50, 3))
  # Response should be shifted kernel starting at position 25
  for (j in 1:3) {
    # Check first few values match (accounting for convolution length)
    n_check <- min(10, 50 - 24)  # Don't go past array bounds
    expect_equal(
      result_delta[25:(24 + n_check - 1), j], 
      kernels[1:(n_check - 1), j]
    )
  }
})

test_that("FFT and C++ convolution produce identical results", {
  set.seed(123)
  
  # Generate test data
  kernel_length <- 20
  kernels <- matrix(exp(-seq(0, 2, length.out = kernel_length)), ncol = 1)
  signal <- numeric(256)
  signal[c(10, 30, 50, 70, 90)] <- 1  # Sparse events
  
  # Force FFT path (large problem)
  result_fft <- fmriparametric:::.fast_batch_convolution(
    signal, kernels, 256
  )
  
  # Force C++ path (small problem)
  signal_short <- signal[1:50]
  result_cpp <- fmriparametric:::.fast_batch_convolution(
    signal_short, kernels, 50
  )
  
  # Compare overlapping region
  expect_equal(
    result_fft[1:50, 1],
    result_cpp[, 1],
    tolerance = 1e-10
  )
})

test_that("Convolution preserves signal energy appropriately", {
  # Test that convolution with unit impulse preserves kernel
  n_time <- 50
  kernel <- exp(-seq(0, 3, length.out = 20))
  kernel <- kernel / sum(kernel)  # Normalize
  
  # Unit impulse
  signal <- numeric(n_time)
  signal[1] <- 1
  
  result <- fmriparametric:::.fast_batch_convolution(
    signal,
    matrix(kernel, ncol = 1),
    n_time
  )
  
  # First part of result should match kernel
  expect_equal(result[1:20, 1], kernel, tolerance = 1e-12)
  
  # Total energy should be preserved
  expect_equal(sum(result), 1, tolerance = 1e-12)
})

test_that("Convolution handles boundary conditions correctly", {
  # Test edge effects
  n_time <- 30
  kernel <- c(0.5, 0.3, 0.2)
  signal <- rep(1, n_time)  # Constant signal
  
  result <- fmriparametric:::.fast_batch_convolution(
    signal,
    matrix(kernel, ncol = 1),
    n_time
  )
  
  # After initial transient, convolution of constant signal with normalized kernel
  # should equal the constant times kernel sum
  kernel_sum <- sum(kernel)
  expect_equal(result[3:n_time, 1], rep(kernel_sum, n_time - 2), tolerance = 1e-12)
  
  # Check initial transient
  expect_equal(result[1, 1], kernel[1])
  expect_equal(result[2, 1], kernel[1] + kernel[2])
})

test_that(".batch_convolution sums multiple signals correctly", {
  # Test the summing behavior
  n_time <- 20
  kernels <- matrix(c(1, 0.5, 0.25), nrow = 3, ncol = 2)
  
  # Two signals
  signals <- matrix(0, nrow = n_time, ncol = 2)
  signals[5, 1] <- 1
  signals[10, 2] <- 1
  
  result <- fmriparametric:::.batch_convolution(
    signals, kernels, n_time
  )
  
  expect_equal(dim(result), c(n_time, 2))
  
  # Check that responses appear at correct times
  expect_equal(result[5:7, 1], kernels[, 1])
  expect_equal(result[10:12, 1], kernels[, 1])
})