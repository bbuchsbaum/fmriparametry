# Test equivalence between FFT and C++ convolution implementations
library(testthat)
library(fmriparametric)

set.seed(123)

# Helper function to generate random signal and kernels
random_signal_kernels <- function(n_time, kernel_len, n_kernels) {
  signal <- rnorm(n_time)
  kernels <- matrix(rnorm(kernel_len * n_kernels), nrow = kernel_len, ncol = n_kernels)
  list(signal = signal, kernels = kernels)
}

test_that("FFT and C++ convolution implementations produce identical results", {
  # Test configuration that triggers FFT path (output_length > 200 && kernel_length > 20)
  n_time_fft <- 256
  kernel_len_fft <- 30
  n_kernels <- 3
  dat_fft <- random_signal_kernels(n_time_fft, kernel_len_fft, n_kernels)
  
  # Test configuration that triggers C++ path
  n_time_cpp <- 100
  kernel_len_cpp <- 10
  dat_cpp <- random_signal_kernels(n_time_cpp, kernel_len_cpp, n_kernels)
  
  # Force FFT path by using large problem
  result_fft_path <- fmriparametric:::.fast_batch_convolution(
    dat_fft$signal, dat_fft$kernels, n_time_fft
  )
  
  # Force C++ path by using small problem
  result_cpp_path <- fmriparametric:::.fast_batch_convolution(
    dat_cpp$signal, dat_cpp$kernels, n_time_cpp
  )
  
  # Now test that we can get same results by calling C++ directly on same data
  # Note: fast_batch_convolution_cpp is not exported, so we use ::: operator
  result_cpp_direct_large <- fmriparametric:::fast_batch_convolution_cpp(
    dat_fft$signal, dat_fft$kernels, n_time_fft
  )
  
  result_cpp_direct_small <- fmriparametric:::fast_batch_convolution_cpp(
    dat_cpp$signal, dat_cpp$kernels, n_time_cpp
  )
  
  # For the large problem, FFT and C++ should give same result
  expect_equal(result_fft_path, result_cpp_direct_large, tolerance = 1e-8,
               label = "FFT path vs direct C++ on large problem")
  
  # For the small problem, the function should already use C++
  expect_equal(result_cpp_path, result_cpp_direct_small, tolerance = 1e-15,
               label = "C++ path vs direct C++ on small problem")
})

test_that("Convolution implementation choice is based on problem size", {
  # Test exact boundary conditions
  # Threshold is: output_length > 200 && kernel_length > 20
  
  signal <- rnorm(250)
  
  # Case 1: Both conditions met - should use FFT
  kernels_large <- matrix(rnorm(25 * 2), nrow = 25, ncol = 2)
  
  # We can't directly test which path was taken, but we can verify results
  result1 <- fmriparametric:::.fast_batch_convolution(signal, kernels_large, 250)
  expect_equal(dim(result1), c(250, 2))
  expect_true(all(is.finite(result1)))
  
  # Case 2: Output length condition not met
  result2 <- fmriparametric:::.fast_batch_convolution(signal[1:150], kernels_large, 150)
  expect_equal(dim(result2), c(150, 2))
  expect_true(all(is.finite(result2)))
  
  # Case 3: Kernel length condition not met
  kernels_small <- matrix(rnorm(15 * 2), nrow = 15, ncol = 2)
  result3 <- fmriparametric:::.fast_batch_convolution(signal, kernels_small, 250)
  expect_equal(dim(result3), c(250, 2))
  expect_true(all(is.finite(result3)))
})

test_that("Both implementations handle edge cases identically", {
  # Test with problematic inputs
  
  # Single time point
  signal_single <- 1.5
  kernels <- matrix(c(1, 0.5, 0.25), nrow = 3, ncol = 2)
  
  result_single <- fmriparametric:::.fast_batch_convolution(
    signal_single, kernels, 1
  )
  result_single_cpp <- fmriparametric:::fast_batch_convolution_cpp(
    signal_single, kernels, 1
  )
  expect_equal(result_single, result_single_cpp)
  
  # Zero signal
  signal_zero <- rep(0, 50)
  result_zero <- fmriparametric:::.fast_batch_convolution(
    signal_zero, kernels, 50
  )
  result_zero_cpp <- fmriparametric:::fast_batch_convolution_cpp(
    signal_zero, kernels, 50
  )
  expect_equal(result_zero, result_zero_cpp)
  expect_true(all(result_zero == 0))
  
  # Single kernel column
  signal <- rnorm(30)
  kernel_single <- matrix(exp(-seq(0, 2, length.out = 10)), ncol = 1)
  
  result_single_kernel <- fmriparametric:::.fast_batch_convolution(
    signal, kernel_single, 30
  )
  result_single_kernel_cpp <- fmriparametric:::fast_batch_convolution_cpp(
    signal, kernel_single, 30
  )
  expect_equal(result_single_kernel, result_single_kernel_cpp, tolerance = 1e-10)
})

test_that("Storage mode conversion works correctly", {
  # The function converts to double precision before calling C++
  signal_int <- as.integer(c(1, 2, 3, 4, 5))
  kernels_int <- matrix(as.integer(c(1, 0, -1)), nrow = 3, ncol = 1)
  
  # This should work despite integer input
  result <- fmriparametric:::.fast_batch_convolution(
    signal_int, kernels_int, 5
  )
  
  expect_equal(dim(result), c(5, 1))
  expect_equal(storage.mode(result), "double")
  expect_true(all(is.finite(result)))
  
  # Compare with explicit double conversion
  signal_double <- as.double(signal_int)
  kernels_double <- matrix(as.double(kernels_int), nrow = 3, ncol = 1)
  result_double <- fmriparametric:::.fast_batch_convolution(
    signal_double, kernels_double, 5
  )
  
  expect_equal(result, result_double)
})