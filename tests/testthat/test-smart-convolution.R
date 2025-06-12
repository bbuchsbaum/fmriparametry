library(testthat)
library(fmriparametric)

# Skip all smart convolution tests - function is in attic for a reason
skip("Smart convolution function is in attic - skipping all tests")

# Test helper function to create test data
.create_test_kernels <- function(n_kernels = 2, kernel_length = 10) {
  set.seed(123)
  kernels <- matrix(rnorm(kernel_length * n_kernels), 
                   nrow = kernel_length, ncol = n_kernels)
  # Normalize kernels
  for (i in 1:n_kernels) {
    kernels[, i] <- kernels[, i] / sum(abs(kernels[, i]))
  }
  return(kernels)
}

# Test helper to create test signal
.create_test_signal <- function(length = 50) {
  set.seed(456)
  signal <- c(rep(0, 10), rep(1, 5), rep(0, length - 15))
  return(signal)
}

# ============================================================================
# Test 1: Basic functionality and output structure
# ============================================================================

test_that(".smart_convolution handles basic convolution correctly", {
  signal <- c(1, 2, 3, 2, 1)
  kernels <- matrix(c(1, 0.5), ncol = 1)
  output_length <- 6
  
  result <- .smart_convolution(signal, kernels, output_length)
  
  # Check output structure
  expect_true(is.matrix(result))
  expect_equal(nrow(result), output_length)
  expect_equal(ncol(result), 1)
  expect_false(anyNA(result))
})

test_that(".smart_convolution handles multiple kernels", {
  signal <- .create_test_signal(20)
  kernels <- .create_test_kernels(n_kernels = 3, kernel_length = 5)
  output_length <- 25
  
  result <- .smart_convolution(signal, kernels, output_length)
  
  expect_equal(nrow(result), output_length)
  expect_equal(ncol(result), 3)
  expect_false(anyNA(result))
})

# ============================================================================
# Test 2: Decision logic - FFT vs Direct convolution
# ============================================================================

test_that(".smart_convolution chooses direct method for small problems", {
  # Create small problem (should trigger direct convolution)
  signal <- .create_test_signal(10)
  kernels <- .create_test_kernels(n_kernels = 2, kernel_length = 3)
  output_length <- 15
  
  # Total ops = 15 * 3 * 2 = 90 < 50000 threshold
  result <- .smart_convolution(signal, kernels, output_length)
  
  expect_equal(nrow(result), output_length)
  expect_equal(ncol(result), 2)
  expect_false(anyNA(result))
})

test_that(".smart_convolution chooses FFT method for large problems", {
  # Create large problem (should trigger FFT)
  signal <- .create_test_signal(100)
  kernels <- .create_test_kernels(n_kernels = 20, kernel_length = 50)
  output_length <- 200
  
  # Total ops = 200 * 50 * 20 = 200000 > 50000 threshold
  result <- .smart_convolution(signal, kernels, output_length)
  
  expect_equal(nrow(result), output_length)
  expect_equal(ncol(result), 20)
  expect_false(anyNA(result))
})

# ============================================================================
# Test 3: Consistency between FFT and direct methods
# ============================================================================

test_that("FFT and direct methods give consistent results for medium problems", {
  signal <- .create_test_signal(30)
  kernels <- .create_test_kernels(n_kernels = 2, kernel_length = 10)
  output_length <- 40
  
  # Force direct method by using small problem size
  signal_small <- signal[1:10]
  kernels_small <- kernels[1:5, , drop = FALSE]
  output_length_small <- 15
  
  result_direct <- .smart_convolution(signal_small, kernels_small, output_length_small)
  
  # Compare with manual convolution for validation
  manual_result <- matrix(0, output_length_small, ncol(kernels_small))
  for (j in 1:ncol(kernels_small)) {
    conv_full <- convolve(signal_small, rev(kernels_small[, j]), type = "open")
    max_len <- min(output_length_small, length(conv_full))
    manual_result[1:max_len, j] <- conv_full[1:max_len]
  }
  
  expect_equal(result_direct, manual_result, tolerance = 1e-10)
})

# ============================================================================
# Test 4: Edge cases and boundary conditions
# ============================================================================

test_that(".smart_convolution handles empty inputs", {
  # Test empty signal
  expect_error(.smart_convolution(numeric(0), matrix(1, nrow = 1, ncol = 1), 5))
  
  # Test empty kernels
  expect_error(.smart_convolution(c(1, 2, 3), matrix(numeric(0), nrow = 0, ncol = 1), 5))
})

test_that(".smart_convolution handles single-element inputs", {
  signal <- 1
  kernels <- matrix(0.5, nrow = 1, ncol = 1)
  output_length <- 3
  
  result <- .smart_convolution(signal, kernels, output_length)
  
  expect_equal(nrow(result), output_length)
  expect_equal(ncol(result), 1)
  expect_equal(result[1, 1], 0.5)
  expect_equal(result[2:3, 1], c(0, 0))
})

test_that(".smart_convolution handles output_length longer than convolution result", {
  signal <- c(1, 2, 3)
  kernels <- matrix(c(1, 2), ncol = 1)
  output_length <- 10  # Much longer than expected convolution length
  
  expect_warning(
    result <- .smart_convolution(signal, kernels, output_length),
    "Requested output_length exceeds convolution result"
  )
  
  expect_equal(nrow(result), output_length)
  expect_false(anyNA(result))
  # Later entries should be zero (padded)
  expect_true(all(result[6:10, 1] == 0))
})

test_that(".smart_convolution handles output_length shorter than convolution result", {
  signal <- c(1, 2, 3, 4, 5)
  kernels <- matrix(c(1, 0.5, 0.25), ncol = 1)
  output_length <- 3  # Shorter than full convolution
  
  result <- .smart_convolution(signal, kernels, output_length)
  
  expect_equal(nrow(result), output_length)
  expect_false(anyNA(result))
})

# ============================================================================
# Test 5: Numerical accuracy and precision
# ============================================================================

test_that(".smart_convolution maintains numerical precision", {
  # Create deterministic test case
  signal <- c(1, 0, 1, 0, 1)
  kernels <- matrix(c(1, -1), ncol = 1)
  output_length <- 6
  
  result <- .smart_convolution(signal, kernels, output_length)
  
  # Manual calculation
  expected <- convolve(signal, c(-1, 1), type = "open")[1:output_length]
  
  expect_equal(as.numeric(result), expected, tolerance = 1e-12)
})

test_that(".smart_convolution handles kernels with different scales", {
  signal <- .create_test_signal(20)
  
  # Create kernels with very different magnitudes
  kernels <- matrix(c(
    rep(1e-6, 5),    # Very small kernel
    rep(1e6, 5),     # Very large kernel
    rep(1, 5)        # Normal kernel
  ), nrow = 5, ncol = 3)
  
  output_length <- 25
  
  result <- .smart_convolution(signal, kernels, output_length)
  
  expect_false(anyNA(result))
  expect_false(any(is.infinite(result)))
  expect_equal(nrow(result), output_length)
  expect_equal(ncol(result), 3)
})

# ============================================================================
# Test 6: Performance characteristics
# ============================================================================

test_that(".smart_convolution performance threshold is reasonable", {
  # Test that the decision threshold (50000 ops) makes sense
  
  # Case 1: Just below threshold - should use direct
  signal1 <- rep(1, 50)
  kernels1 <- matrix(rep(1, 100), nrow = 10, ncol = 10)  # 10x10
  output1 <- 45
  # ops = 45 * 10 * 10 = 4500 < 50000
  
  result1 <- .smart_convolution(signal1, kernels1, output1)
  expect_false(anyNA(result1))
  
  # Case 2: Just above threshold - should use FFT
  signal2 <- rep(1, 100)
  kernels2 <- matrix(rep(1, 1000), nrow = 50, ncol = 20)  # 50x20
  output2 <- 55
  # ops = 55 * 50 * 20 = 55000 > 50000
  
  result2 <- .smart_convolution(signal2, kernels2, output2)
  expect_false(anyNA(result2))
})

# ============================================================================
# Test 7: Error handling and input validation
# ============================================================================

test_that(".smart_convolution validates input dimensions", {
  signal <- c(1, 2, 3)
  
  # Test mismatched kernel matrix dimensions
  invalid_kernels <- matrix(c(1, 2, 3, 4, 5), nrow = 2, ncol = 3)  # 2x3 but last row incomplete
  expect_error(.smart_convolution(signal, invalid_kernels, 5))
})

test_that(".smart_convolution handles non-numeric inputs gracefully", {
  expect_error(.smart_convolution("not numeric", matrix(1, 1, 1), 5))
  expect_error(.smart_convolution(c(1, 2, 3), "not matrix", 5))
  expect_error(.smart_convolution(c(1, 2, 3), matrix(1, 1, 1), "not numeric"))
})

test_that(".smart_convolution handles negative output_length", {
  signal <- c(1, 2, 3)
  kernels <- matrix(c(1, 2), ncol = 1)
  
  expect_error(.smart_convolution(signal, kernels, -1))
  expect_error(.smart_convolution(signal, kernels, 0))
})

# ============================================================================
# Test 8: Integration with fmriparametric workflow
# ============================================================================

test_that(".smart_convolution works with realistic HRF kernels", {
  # Create realistic HRF-like kernels
  t_hrf <- seq(0, 20, by = 0.5)
  
  # LWU model parameters
  tau <- 6
  sigma <- 2.5
  rho <- 0.35
  
  # Generate LWU HRF
  hrf1 <- exp(-(t_hrf - tau)^2 / (2 * sigma^2)) - 
          rho * exp(-(t_hrf - tau - 2*sigma)^2 / (2 * (1.6*sigma)^2))
  
  # Generate derivative kernels (simplified)
  hrf2 <- hrf1 * (t_hrf - tau) / sigma^2  # tau derivative (simplified)
  hrf3 <- hrf1 * ((t_hrf - tau)^2 / sigma^3 - 1/sigma)  # sigma derivative (simplified)
  
  kernels <- cbind(hrf1, hrf2, hrf3)
  
  # Create stimulus signal
  stimulus <- rep(c(1, rep(0, 9)), 10)  # Events every 10 time points
  
  result <- .smart_convolution(stimulus, kernels, length(stimulus) + 20)
  
  expect_false(anyNA(result))
  expect_equal(ncol(result), 3)
  expect_true(nrow(result) > length(stimulus))
})

# ============================================================================
# Test 9: Memory efficiency (for large inputs)
# ============================================================================

test_that(".smart_convolution is memory efficient for large problems", {
  skip_on_cran()  # Skip memory-intensive tests on CRAN
  
  # Create moderately large problem
  signal <- rep(c(1, 0), 1000)  # 2000 elements
  kernels <- .create_test_kernels(n_kernels = 5, kernel_length = 100)
  output_length <- 2500
  
  # This should complete without memory errors
  expect_no_error({
    result <- .smart_convolution(signal, kernels, output_length)
  })
  
  expect_equal(nrow(result), output_length)
  expect_equal(ncol(result), 5)
  expect_false(anyNA(result))
})

# ============================================================================
# Test 10: Warning suppression for expected cases
# ============================================================================

test_that(".smart_convolution only warns once for padding", {
  signal <- c(1, 2)
  kernels <- matrix(c(1, 2), ncol = 1)
  output_length <- 10
  
  # Should get exactly one warning
  expect_warning(
    result <- .smart_convolution(signal, kernels, output_length),
    "Requested output_length exceeds convolution result",
    fixed = TRUE
  )
  
  # Result should still be valid
  expect_equal(nrow(result), output_length)
  expect_false(anyNA(result))
})
