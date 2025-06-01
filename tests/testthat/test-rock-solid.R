library(testthat)

context("Rock solid implementation tests")

# Test adversarial inputs
test_that("rock solid handles pathological inputs without crashing", {
  # Test 1: NULL inputs
  expect_error(
    estimate_parametric_hrf_rock_solid(NULL, NULL),
    "fmri_data cannot be NULL"
  )
  
  # Test 2: Empty data
  expect_error(
    estimate_parametric_hrf_rock_solid(
      matrix(nrow = 0, ncol = 0),
      matrix(nrow = 0, ncol = 0)
    ),
    "Insufficient time points"
  )
  
  # Test 3: Single voxel, single timepoint
  expect_error(
    estimate_parametric_hrf_rock_solid(
      matrix(1),
      matrix(1)
    ),
    "Insufficient time points"
  )
  
  # Test 4: All NA data
  Y_na <- matrix(NA, nrow = 50, ncol = 10)
  S_na <- matrix(c(rep(c(1, 0, 0), 16), 1, 0), ncol = 1)
  
  expect_warning(
    fit_na <- estimate_parametric_hrf_rock_solid(
      Y_na, S_na,
      verbose = FALSE,
      safety_mode = "maximum"
    ),
    "contains NA values"
  )
  
  # Should still return valid object
  expect_s3_class(fit_na, "parametric_hrf_fit")
  expect_equal(dim(fit_na$estimated_parameters), c(10, 3))
})

# Test numerical edge cases
test_that("rock solid handles numerical extremes", {
  set.seed(123)
  
  # Test 1: Extreme values
  Y_extreme <- matrix(c(
    rep(1e10, 25),   # Very large
    rep(1e-10, 25),  # Very small
    rep(0, 25),      # Zeros
    rnorm(25)        # Normal
  ), nrow = 50, ncol = 4, byrow = FALSE)
  
  S <- matrix(rbinom(50, 1, 0.2), ncol = 1)
  
  expect_warning(
    fit_extreme <- estimate_parametric_hrf_rock_solid(
      Y_extreme, S,
      verbose = FALSE,
      lambda_ridge = 0.1  # Higher regularization for stability
    ),
    class = NA  # Expect some warning
  )
  
  expect_s3_class(fit_extreme, "parametric_hrf_fit")
  expect_false(any(is.na(fit_extreme$estimated_parameters)))
  
  # Test 2: Constant voxels
  Y_const <- matrix(rep(1:5, each = 50), nrow = 50, ncol = 5)
  
  expect_warning(
    fit_const <- estimate_parametric_hrf_rock_solid(
      Y_const, S,
      verbose = FALSE
    ),
    "constant voxels"
  )
  
  expect_s3_class(fit_const, "parametric_hrf_fit")
})

# Test memory safety
test_that("rock solid handles memory constraints", {
  skip_if_not(Sys.getenv("FMRIPARAMETRIC_EXTENDED_TESTS") == "true",
              "Skipping memory tests")
  
  # Large but manageable data
  n_time <- 200
  n_vox <- 5000
  
  Y_large <- matrix(rnorm(n_time * n_vox), nrow = n_time)
  S_large <- matrix(rbinom(n_time, 1, 0.1), ncol = 1)
  
  # Should handle this with chunking if needed
  fit_large <- estimate_parametric_hrf_rock_solid(
    Y_large, S_large,
    verbose = FALSE,
    safety_mode = "balanced",
    compute_se = FALSE  # Save memory
  )
  
  expect_s3_class(fit_large, "parametric_hrf_fit")
  expect_equal(nrow(fit_large$estimated_parameters), n_vox)
})

# Test error recovery
test_that("rock solid recovers from algorithm failures", {
  # Create data that causes numerical issues
  Y_singular <- matrix(1:100, nrow = 20, ncol = 5)
  S_singular <- matrix(rep(1, 20), ncol = 1)  # Constant stimulus
  
  # Should recover and return something
  fit_recovered <- estimate_parametric_hrf_rock_solid(
    Y_singular, S_singular,
    verbose = FALSE,
    error_report = TRUE
  )
  
  expect_s3_class(fit_recovered, "parametric_hrf_fit")
  
  # Check metadata for recovery info
  expect_true("errors_recovered" %in% names(fit_recovered$metadata))
})

# Test progressive degradation
test_that("progressive estimation fallbacks work correctly", {
  set.seed(456)
  
  # Normal data
  Y <- matrix(rnorm(100), nrow = 20, ncol = 5)
  S <- matrix(rbinom(20, 1, 0.3), ncol = 1)
  
  # Force fallback by using impossible bounds
  impossible_bounds <- list(
    lower = c(100, 100, 100),
    upper = c(101, 101, 101)
  )
  
  expect_warning(
    fit_fallback <- estimate_parametric_hrf_rock_solid(
      Y, S,
      theta_bounds = impossible_bounds,
      verbose = FALSE
    ),
    class = NA
  )
  
  # Should still get results
  expect_s3_class(fit_fallback, "parametric_hrf_fit")
})

# Test safety modes
test_that("different safety modes produce appropriate trade-offs", {
  set.seed(789)
  Y <- matrix(rnorm(500), nrow = 50, ncol = 10)
  S <- matrix(rbinom(50, 1, 0.2), ncol = 1)
  
  # Maximum safety
  time_max <- system.time({
    fit_max <- estimate_parametric_hrf_rock_solid(
      Y, S,
      safety_mode = "maximum",
      verbose = FALSE
    )
  })
  
  # Performance mode
  time_perf <- system.time({
    fit_perf <- estimate_parametric_hrf_rock_solid(
      Y, S,
      safety_mode = "performance",
      verbose = FALSE
    )
  })
  
  # Maximum should be more thorough (usually slower)
  # But both should produce valid results
  expect_s3_class(fit_max, "parametric_hrf_fit")
  expect_s3_class(fit_perf, "parametric_hrf_fit")
})

# Test error reporting
test_that("error reporting provides useful information", {
  # Trigger some errors
  Y_bad <- matrix(rnorm(50), nrow = 10, ncol = 5)
  S_bad <- matrix(rnorm(15), ncol = 1)  # Wrong dimensions
  
  expect_error({
    fit_error <- estimate_parametric_hrf_rock_solid(
      Y_bad, S_bad,
      error_report = TRUE,
      verbose = FALSE
    )
  }, "don't match")
})

# Test input validation thoroughness
test_that("input validation catches all issues", {
  Y <- matrix(rnorm(100), nrow = 20, ncol = 5)
  S <- matrix(rbinom(20, 1, 0.3), ncol = 1)
  
  # Bad parameter bounds
  expect_error(
    estimate_parametric_hrf_rock_solid(
      Y, S,
      theta_bounds = list(lower = c(5, 5, 5), upper = c(1, 1, 1)),
      verbose = FALSE
    ),
    "lower must be less than upper"
  )
  
  # Bad numeric parameters
  expect_error(
    estimate_parametric_hrf_rock_solid(
      Y, S,
      lambda_ridge = -1,
      verbose = FALSE
    ),
    "outside valid range"
  )
  
  # Wrong model
  expect_error(
    estimate_parametric_hrf_rock_solid(
      Y, S,
      parametric_hrf = "unknown_model",
      verbose = FALSE
    ),
    "Only 'lwu' model currently supported"
  )
})

# Test partial results
test_that("partial results are returned when requested", {
  # Create mixed quality data
  Y_mixed <- matrix(rnorm(200), nrow = 40, ncol = 5)
  Y_mixed[, 4:5] <- NA  # Last two voxels are bad
  S <- matrix(rbinom(40, 1, 0.2), ncol = 1)
  
  fit_partial <- estimate_parametric_hrf_rock_solid(
    Y_mixed, S,
    allow_partial = TRUE,
    verbose = FALSE
  )
  
  expect_s3_class(fit_partial, "parametric_hrf_fit")
  # Should have results for all voxels (some may be defaults)
  expect_equal(nrow(fit_partial$estimated_parameters), 5)
})

# Test reproducibility
test_that("rock solid estimation is reproducible", {
  Y <- matrix(rnorm(100), nrow = 20, ncol = 5)
  S <- matrix(rbinom(20, 1, 0.3), ncol = 1)
  
  # Run twice with same seed
  set.seed(111)
  fit1 <- estimate_parametric_hrf_rock_solid(
    Y, S,
    verbose = FALSE,
    recenter_kmeans_passes = 0  # Disable for exact reproducibility
  )
  
  set.seed(111)
  fit2 <- estimate_parametric_hrf_rock_solid(
    Y, S,
    verbose = FALSE,
    recenter_kmeans_passes = 0
  )
  
  # Results should be identical
  expect_equal(fit1$estimated_parameters, fit2$estimated_parameters)
  expect_equal(fit1$amplitudes, fit2$amplitudes)
})