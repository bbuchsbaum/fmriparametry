# Edge case tests for HRF interface functions
library(fmriparametric)
library(testthat)

test_that(".lwu_hrf_function clamps parameters to bounds", {
  t <- seq(0, 30, by = 1)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  
  # Test parameters below bounds
  params_low <- c(-1, 0.01, -0.1)
  res_low <- fmriparametric:::.lwu_hrf_function(t, params_low, bounds)
  
  # Expected: parameters should be clamped to bounds
  expected_low <- fmrihrf::hrf_lwu(
    t = t,
    tau = bounds$lower[1],      # 0
    sigma = 0.051,               # Special handling for sigma
    rho = bounds$lower[3],       # 0
    normalize = "none"
  )
  expect_equal(res_low, expected_low)
  
  # Test parameters above bounds
  params_high <- c(25, 11, 2)
  res_high <- fmriparametric:::.lwu_hrf_function(t, params_high, bounds)
  
  expected_high <- fmrihrf::hrf_lwu(
    t = t,
    tau = bounds$upper[1],       # 20
    sigma = bounds$upper[2],     # 10
    rho = bounds$upper[3],       # 1.5
    normalize = "none"
  )
  expect_equal(res_high, expected_high)
})

test_that(".lwu_hrf_function forces sigma > 0.05", {
  t <- seq(0, 2, length.out = 3)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  
  # Test exact boundary value for sigma
  res <- fmriparametric:::.lwu_hrf_function(t, c(6, 0.05, 0.3), bounds)
  
  expected <- fmrihrf::hrf_lwu(
    t = t,
    tau = 6,
    sigma = 0.051,  # Should be forced to 0.051
    rho = 0.3,
    normalize = "none"
  )
  expect_equal(res, expected)
  
  # Test value below 0.05
  res2 <- fmriparametric:::.lwu_hrf_function(t, c(6, 0.01, 0.3), bounds)
  expect_equal(res2, expected)  # Should give same result
})

test_that(".lwu_hrf_function validates parameter length", {
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  
  # Wrong length should trigger an error
  expect_error(
    fmriparametric:::.lwu_hrf_function(seq(0, 1), c(1, 2), bounds),
    "length\\(params_vector\\) == 3"
  )
  
  # Also test with too many parameters
  expect_error(
    fmriparametric:::.lwu_hrf_function(seq(0, 1), c(1, 2, 3, 4), bounds),
    "length\\(params_vector\\) == 3"
  )
})

test_that(".lwu_hrf_taylor_basis_function clamps parameters and returns matrix", {
  t_hrf <- seq(0, 30, by = 1)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  
  # Test with out-of-bounds parameters
  params <- c(-2, 0.01, 2)
  basis <- fmriparametric:::.lwu_hrf_taylor_basis_function(params, t_hrf, bounds)
  
  # Check return type
  expect_true(is.matrix(basis))
  expect_equal(ncol(basis), 4)  # HRF + 3 derivatives
  expect_equal(nrow(basis), length(t_hrf))
  expect_equal(storage.mode(basis), "double")
  
  # Verify parameters were clamped appropriately
  # The function clamps to bounds, then adds epsilon to avoid exact boundaries
  trimmed <- pmax(bounds$lower, pmin(params, bounds$upper))
  eps <- c(0.01, 0.01, 0.01)
  trimmed <- pmax(bounds$lower + eps, pmin(trimmed, bounds$upper - eps))
  names(trimmed) <- c("tau", "sigma", "rho")
  
  expected <- fmrihrf::hrf_basis_lwu(
    theta0 = trimmed,
    t = t_hrf,
    normalize_primary = "none"
  )
  
  if (!is.matrix(expected)) {
    expected <- matrix(expected, ncol = 4)
  }
  storage.mode(expected) <- "double"
  
  expect_equal(basis, expected)
})

test_that(".lwu_hrf_taylor_basis_function validates parameter length", {
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  
  expect_error(
    fmriparametric:::.lwu_hrf_taylor_basis_function(c(1, 2), seq(0, 1), bounds),
    "length\\(params_vector0\\) == 3"
  )
})

test_that("HRF interface through registry works correctly", {
  # Test the actual interface as used in the package
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  
  t <- seq(0, 10, by = 0.5)
  params <- c(6, 3, 0.5)
  
  # Test HRF function through interface
  hrf_vals <- hrf_interface$hrf_function(t, params)
  expect_equal(length(hrf_vals), length(t))
  expect_true(all(is.finite(hrf_vals)))
  
  # Test Taylor basis through interface
  basis_vals <- hrf_interface$taylor_basis(params, t)
  expect_true(is.matrix(basis_vals))
  expect_equal(dim(basis_vals), c(length(t), 4))
  expect_true(all(is.finite(basis_vals)))
  
  # Test other interface components
  expect_equal(hrf_interface$parameter_names, c("tau", "sigma", "rho"))
  expect_equal(hrf_interface$default_seed(), c(6, 2.5, 0.35))
  expect_equal(names(hrf_interface$default_bounds()), c("lower", "upper"))
})

test_that("Boundary parameters produce valid HRF values", {
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  bounds <- hrf_interface$default_bounds()
  t <- seq(0, 30, length.out = 100)
  
  # Test at exact lower bounds (adjusted for sigma)
  params_lower <- c(bounds$lower[1], 0.1, bounds$lower[3])
  hrf_lower <- hrf_interface$hrf_function(t, params_lower)
  expect_true(all(is.finite(hrf_lower)))
  expect_false(all(hrf_lower == 0))  # Should produce non-zero HRF
  
  # Test at exact upper bounds
  params_upper <- bounds$upper
  hrf_upper <- hrf_interface$hrf_function(t, params_upper)
  expect_true(all(is.finite(hrf_upper)))
  expect_false(all(hrf_upper == 0))
  
  # Test mixed boundary conditions
  params_mixed <- c(bounds$lower[1], 5, bounds$upper[3])
  hrf_mixed <- hrf_interface$hrf_function(t, params_mixed)
  expect_true(all(is.finite(hrf_mixed)))
})

test_that("HRF interface handles invalid model names", {
  expect_error(
    fmriparametric:::.create_hrf_interface("invalid_model"),
    "HRF model 'invalid_model' is not registered"
  )
})