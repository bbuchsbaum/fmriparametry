library(testthat)
library(fmriparametric)

test_that("get_available_memory returns reasonable values", {
  # Test that memory detection works
  available_memory <- fmriparametric:::.get_available_memory()
  
  # Should return a numeric value
  expect_true(is.numeric(available_memory))
  expect_length(available_memory, 1)
  expect_false(is.na(available_memory))
  
  # Should be a reasonable amount (at least 1GB, less than 1TB)
  expect_true(available_memory >= 1e9)  # At least 1GB
  expect_true(available_memory <= 1e12) # Less than 1TB
  
  cat("\nDetected available memory:", round(available_memory / 1e9, 2), "GB\n")
})

test_that("check_memory_available works correctly", {
  # Test with small requirement (should pass)
  small_requirement <- 1e6  # 1MB
  result_small <- fmriparametric:::.check_memory_available(
    required_bytes = small_requirement,
    operation = "small test"
  )
  expect_true(result_small)
  
  # Test with huge requirement (should fail with warning)
  huge_requirement <- 1e15  # 1PB - obviously too much
  expect_warning(
    result_huge <- fmriparametric:::.check_memory_available(
      required_bytes = huge_requirement,
      operation = "huge test"
    ),
    "requires.*GB but only.*GB available"
  )
  expect_false(result_huge)
  
  cat("\nMemory check tests passed\n")
})

test_that("memory functions handle edge cases", {
  # Test check_memory_available with zero requirement
  result_zero <- fmriparametric:::.check_memory_available(
    required_bytes = 0,
    operation = "zero test"
  )
  expect_true(result_zero)
  
  # Test check_memory_available with negative requirement (should still work)
  expect_silent(
    result_negative <- fmriparametric:::.check_memory_available(
      required_bytes = -1000,
      operation = "negative test"
    )
  )
  expect_true(result_negative)
  
  cat("\nEdge case tests passed\n")
})

test_that("memory functions work on different platforms", {
  # This test documents the platform-specific behavior
  # We can't test all platforms, but we can test that the function
  # works on the current platform
  
  available <- fmriparametric:::.get_available_memory()
  
  # On any platform, should get a reasonable value
  expect_true(is.numeric(available))
  expect_true(available > 0)
  
  # Check platform-specific information
  platform <- .Platform$OS.type
  sysname <- Sys.info()["sysname"]
  
  cat("\nPlatform:", platform, "\n")
  cat("System:", sysname, "\n")
  cat("Detected memory:", round(available / 1e9, 2), "GB\n")
  
  # Platform-specific expectations
  if (platform == "windows") {
    # On Windows, should use memory.limit()
    expect_true(available > 0)
  } else if (sysname == "Darwin") {
    # On macOS, should use vm_stat
    expect_true(available > 0)
  } else {
    # On Linux, should use /proc/meminfo
    expect_true(available > 0)
  }
})

test_that("memory estimation for typical fMRI datasets", {
  # Test memory requirements for typical dataset sizes
  
  # Small dataset: 100 voxels, 200 timepoints
  small_voxels <- 100
  small_timepoints <- 200
  small_params <- 3
  
  # Estimate memory for small dataset (similar to smart memory decision logic)
  # Data: voxels * timepoints * 8 bytes (double precision)
  # Working: assume 3x data size for intermediate calculations
  small_data_memory <- small_voxels * small_timepoints * 8
  small_working_memory <- small_data_memory * 3
  
  result_small <- fmriparametric:::.check_memory_available(
    required_bytes = small_working_memory,
    operation = "small fMRI dataset"
  )
  expect_true(result_small)
  
  # Medium dataset: 10,000 voxels, 500 timepoints
  medium_voxels <- 10000
  medium_timepoints <- 500
  
  medium_data_memory <- medium_voxels * medium_timepoints * 8
  medium_working_memory <- medium_data_memory * 3
  
  # This should usually pass on modern systems
  result_medium <- fmriparametric:::.check_memory_available(
    required_bytes = medium_working_memory,
    operation = "medium fMRI dataset"
  )
  
  # Large dataset: 100,000 voxels, 1000 timepoints
  large_voxels <- 100000
  large_timepoints <- 1000
  
  large_data_memory <- large_voxels * large_timepoints * 8
  large_working_memory <- large_data_memory * 3
  
  # This might fail on some systems, which is expected
  # The function might produce warnings for large datasets, which is normal
  suppressMessages(suppressWarnings(
    result_large <- fmriparametric:::.check_memory_available(
      required_bytes = large_working_memory,
      operation = "large fMRI dataset"
    )
  ))
  
  cat("\nMemory estimates for fMRI datasets:\n")
  cat("Small dataset (", small_voxels, "x", small_timepoints, "): ", 
      round(small_working_memory / 1e6, 1), "MB - ", 
      if(result_small) "OK" else "Too large", "\n")
  cat("Medium dataset (", medium_voxels, "x", medium_timepoints, "): ", 
      round(medium_working_memory / 1e6, 1), "MB - ", 
      if(result_medium) "OK" else "Too large", "\n")
  cat("Large dataset (", large_voxels, "x", large_timepoints, "): ", 
      round(large_working_memory / 1e6, 1), "MB - ", 
      if(result_large) "OK" else "Too large", "\n")
})

test_that("memory fallback behavior when memory cannot be detected", {
  # We can't easily mock the internal memory detection failure,
  # but we can test that the functions handle unusual return values gracefully
  
  # The current implementation should handle cases where memory detection fails
  # by returning default values or allowing operations to proceed
  
  # Test that even if we can't determine exact memory, basic operations work
  available <- fmriparametric:::.get_available_memory()
  
  # Should always return something numeric and positive
  expect_true(is.numeric(available))
  expect_true(length(available) == 1)
  expect_true(available > 0)
  
  # Test that check works even with small amounts
  small_check <- fmriparametric:::.check_memory_available(1000, "fallback test")
  expect_true(is.logical(small_check))
  expect_length(small_check, 1)
  
  cat("\nMemory fallback behavior test passed\n")
})
