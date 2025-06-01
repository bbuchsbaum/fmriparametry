library(testthat)

context("Sprint 2 S3 methods")

# Create test fit object with Sprint 2 features
create_test_fit_v2 <- function(n_vox = 10, n_time = 50) {
  set.seed(999)
  
  params <- matrix(runif(n_vox * 3, 
                         min = c(3, 1, 0), 
                         max = c(9, 4, 1)), 
                   nrow = n_vox, ncol = 3)
  colnames(params) <- c("tau", "sigma", "rho")
  
  amps <- runif(n_vox, 0.5, 2.5)
  r2 <- runif(n_vox, 0.1, 0.9)
  residuals <- matrix(rnorm(n_time * n_vox, sd = 0.1), nrow = n_time)
  
  ses <- matrix(runif(n_vox * 3, 0.1, 0.5), nrow = n_vox, ncol = 3)
  colnames(ses) <- c("tau", "sigma", "rho")
  
  conv_info <- list(
    trajectory = list(c(6, 2.5, 0.35), c(5.5, 2.3, 0.33), c(5.2, 2.2, 0.32)),
    n_iterations = 3,
    final_global_theta = c(5.2, 2.2, 0.32),
    converged = TRUE
  )
  
  new_parametric_hrf_fit(
    estimated_parameters = params,
    amplitudes = amps,
    parameter_names = c("tau", "sigma", "rho"),
    hrf_model = "lwu",
    r_squared = r2,
    residuals = residuals,
    parameter_ses = ses,
    convergence_info = conv_info,
    metadata = list(
      n_voxels = n_vox,
      n_timepoints = n_time,
      theta_seed = c(6, 2.5, 0.35),
      recenter_global_passes = 3
    )
  )
}

# Test enhanced print method
test_that("print method shows Sprint 2 information", {
  fit <- create_test_fit_v2()
  
  output <- capture.output(print(fit))
  
  expect_true(any(grepl("Mean R²:", output)))
  expect_true(any(grepl("Global iterations:", output)))
  expect_true(any(grepl("Parametric HRF Fit", output)))
})

# Test enhanced summary method
test_that("summary method includes Sprint 2 fields", {
  fit <- create_test_fit_v2()
  
  summ <- summary(fit)
  
  # Check structure
  expect_s3_class(summ, "summary.parametric_hrf_fit")
  expect_true(!is.null(summ$r_squared_summary))
  expect_true(!is.null(summ$parameter_se_summary))
  expect_true(!is.null(summ$convergence_info))
  
  # Check convergence info
  expect_equal(summ$convergence_info$n_iterations, 3)
  expect_true(summ$convergence_info$converged)
  expect_equal(length(summ$convergence_info$delta_theta), 3)
  
  # Check summaries have correct structure
  expect_equal(dim(summ$parameter_summary), c(6, 3))
  expect_equal(dim(summ$parameter_se_summary), c(6, 3))
  expect_length(summ$r_squared_summary, 6)
  
  # Test print method for summary
  output <- capture.output(print(summ))
  expect_true(any(grepl("R-squared Summary:", output)))
  expect_true(any(grepl("Parameter Standard Errors:", output)))
  expect_true(any(grepl("Convergence Information:", output)))
})

# Test enhanced coef method
test_that("coef method supports type argument", {
  fit <- create_test_fit_v2()
  
  # Default: parameters
  params <- coef(fit)
  expect_equal(dim(params), c(10, 3))
  expect_equal(colnames(params), c("tau", "sigma", "rho"))
  
  # Amplitudes
  amps <- coef(fit, type = "amplitude")
  expect_length(amps, 10)
  expect_true(all(amps > 0))
  
  # Standard errors
  ses <- coef(fit, type = "se")
  expect_equal(dim(ses), c(10, 3))
  expect_true(all(ses > 0))
  
  # Test with missing SEs
  fit_no_se <- fit
  fit_no_se$parameter_ses <- NULL
  expect_warning(se_null <- coef(fit_no_se, type = "se"))
  expect_null(se_null)
})

# Test fitted method
test_that("fitted method works correctly", {
  fit <- create_test_fit_v2(n_vox = 5, n_time = 20)
  
  # Create fake Y data
  Y_proj <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  
  # Get fitted values
  fitted_vals <- fitted(fit, Y_proj = Y_proj)
  
  expect_equal(dim(fitted_vals), dim(Y_proj))
  
  # fitted + residuals = Y
  expect_equal(fitted_vals + fit$residuals, Y_proj, tolerance = 1e-10)
  
  # Error without Y_proj
  expect_error(fitted(fit), "Y_proj required")
  
  # Error without residuals
  fit_no_resid <- fit
  fit_no_resid$residuals <- NULL
  expect_error(fitted(fit_no_resid, Y_proj), "Cannot compute fitted values")
})

# Test residuals method
test_that("residuals method works correctly", {
  fit <- create_test_fit_v2(n_vox = 5, n_time = 20)
  
  resid <- residuals(fit)
  expect_equal(dim(resid), c(20, 5))
  
  # Warning when no residuals
  fit_no_resid <- fit
  fit_no_resid$residuals <- NULL
  expect_warning(resid_null <- residuals(fit_no_resid))
  expect_null(resid_null)
})

# Test vcov method
test_that("vcov method works for individual voxels", {
  fit <- create_test_fit_v2()
  
  # Default: first voxel
  vcov1 <- vcov(fit)
  expect_equal(dim(vcov1), c(3, 3))
  expect_true(all(diag(vcov1) > 0))
  expect_equal(diag(vcov1), fit$parameter_ses[1, ]^2)
  
  # Specific voxel
  vcov5 <- vcov(fit, voxel_index = 5)
  expect_equal(diag(vcov5), fit$parameter_ses[5, ]^2)
  
  # Invalid voxel
  expect_error(vcov(fit, voxel_index = 100))
  expect_error(vcov(fit, voxel_index = 0))
  
  # No SEs
  fit_no_se <- fit
  fit_no_se$parameter_ses <- NULL
  expect_warning(vcov_null <- vcov(fit_no_se))
  expect_null(vcov_null)
})

# Test plot method
test_that("plot method produces output without errors", {
  fit <- create_test_fit_v2(n_vox = 20)
  
  # Default plot (median HRF)
  expect_silent(plot(fit, hrf_time_max = 20))
  
  # Multiple voxels
  expect_silent(plot(fit, voxel_indices = c(1, 5, 10)))
  
  # Parameter distributions
  expect_silent(plot(fit, type = "parameters"))
  
  # Many voxels (should warn)
  expect_warning(plot(fit, voxel_indices = 1:30))
})

# Test backward compatibility
test_that("v2 methods work with v1 objects", {
  # Create v1-style object (no Sprint 2 fields)
  params <- matrix(runif(15), nrow = 5, ncol = 3)
  colnames(params) <- c("tau", "sigma", "rho")
  
  fit_v1 <- structure(
    list(
      estimated_parameters = params,
      amplitudes = runif(5),
      parameter_names = c("tau", "sigma", "rho"),
      hrf_model = "lwu",
      convergence = list(),
      metadata = list(n_voxels = 5, n_timepoints = 50)
    ),
    class = "parametric_hrf_fit"
  )
  
  # Should work without errors
  expect_output(print(fit_v1))
  expect_s3_class(summary(fit_v1), "summary.parametric_hrf_fit")
  expect_equal(dim(coef(fit_v1)), c(5, 3))
  
  # Should handle missing fields gracefully
  summ_v1 <- summary(fit_v1)
  expect_null(summ_v1$r_squared_summary)
  expect_null(summ_v1$parameter_se_summary)
})

# Test is_v2_fit helper
test_that("is_v2_fit correctly identifies Sprint 2 objects", {
  source(test_path("../../R/parametric-hrf-fit-class-v2.R"))
  
  fit_v2 <- create_test_fit_v2()
  expect_true(is_v2_fit(fit_v2))
  
  # Remove one v2 field at a time
  fit_partial <- fit_v2
  fit_partial$r_squared <- NULL
  expect_true(is_v2_fit(fit_partial))  # Still has other v2 fields
  
  fit_partial$residuals <- NULL
  expect_true(is_v2_fit(fit_partial))  # Still has SEs
  
  fit_partial$parameter_ses <- NULL
  expect_false(is_v2_fit(fit_partial))  # No v2 fields left
})