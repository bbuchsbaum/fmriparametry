test_that("S3 methods work for parametric_hrf_fit objects", {
  # Create a fit object
  set.seed(123)
  fmri_data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  event_design <- matrix(rbinom(20, 1, 0.2), ncol = 1)
  
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    verbose = FALSE
  )
  
  # Test print method
  expect_output(print(fit), "Parametric HRF Fit")
  expect_output(print(fit), "Model: lwu")
  expect_output(print(fit), "Voxels: 5")
  
  # Test coef method
  coefficients <- coef(fit)
  expect_equal(dim(coefficients), c(5, 3))
  expect_true(is.matrix(coefficients))
  expect_equal(rownames(coefficients), paste0("Voxel_", 1:5))
  expect_equal(colnames(coefficients), c("tau", "sigma", "rho"))
  
  # Test summary method
  summary_obj <- summary(fit)
  expect_s3_class(summary_obj, "summary_parametric_hrf_fit")
  expect_output(print(summary_obj), "Summary of Parametric HRF Fit")
  
  # Test fitted method
  # Need to pass Y_proj since fit doesn't store fitted values
  fitted_values <- fitted(fit, Y_proj = fmri_data)
  expect_true(is.matrix(fitted_values))
  expect_equal(dim(fitted_values), c(20, 5))
  
  # Test residuals method
  resid_values <- residuals(fit)
  expect_true(is.matrix(resid_values))
  expect_equal(dim(resid_values), c(20, 5))
  
  # Check that fitted + residuals = original data (approximately)
  reconstructed <- fitted_values + resid_values
  expect_equal(reconstructed, fmri_data, tolerance = 1e-10)
})

test_that("plot method works for parametric_hrf_fit", {
  skip_if_not_installed("ggplot2")
  
  set.seed(456)
  fmri_data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  event_design <- matrix(rbinom(20, 1, 0.2), ncol = 1)
  
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    verbose = FALSE
  )
  
  # Test that plot returns a ggplot object
  p <- plot(fit)
  expect_s3_class(p, "ggplot")
  
  # Test plotting specific voxels
  p_subset <- plot(fit, voxels = c(1, 3))
  expect_s3_class(p_subset, "ggplot")
  
  # Test plotting with custom time range
  p_custom <- plot(fit, time_range = c(0, 20))
  expect_s3_class(p_custom, "ggplot")
})