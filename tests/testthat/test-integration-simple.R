# Integration tests for core functionality

test_that("estimate_parametric_hrf integration with minimal data", {
  set.seed(42)
  
  # Very simple test data
  n_time <- 50
  n_vox <- 5
  
  # Simple synthetic fMRI data
  fmri_data <- matrix(rnorm(n_time * n_vox, mean = 100, sd = 10), 
                      nrow = n_time, ncol = n_vox)
  
  # Simple event model - just a few events
  event_model <- matrix(0, nrow = n_time, ncol = 1)
  event_model[c(10, 25, 40), 1] <- 1
  
  # Should run without errors
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  # Basic structure checks
  expect_s3_class(fit, "parametric_hrf_fit")
  expect_true(is.matrix(coef(fit)))
  expect_equal(nrow(coef(fit)), n_vox)
  expect_equal(ncol(coef(fit)), 3)
  expect_equal(colnames(coef(fit)), c("tau", "sigma", "rho"))
  
  # Parameters should be within LWU bounds
  params <- coef(fit)
  expect_true(all(params[, "tau"] >= 0 & params[, "tau"] <= 20))
  expect_true(all(params[, "sigma"] >= 0.05 & params[, "sigma"] <= 10))
  expect_true(all(params[, "rho"] >= 0 & params[, "rho"] <= 1.5))
})

test_that("S3 methods work correctly", {
  set.seed(123)
  
  # Generate test data
  fmri_data <- matrix(rnorm(30 * 3), nrow = 30, ncol = 3)
  event_model <- matrix(rbinom(30, 1, 0.2), ncol = 1)
  
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  # Test print method
  expect_output(print(fit), "Parametric HRF Fit")
  
  # Test coef method  
  params <- coef(fit)
  expect_true(is.matrix(params))
  expect_equal(dim(params), c(3, 3))
  
  # Test summary method
  summary_fit <- summary(fit)
  expect_s3_class(summary_fit, "summary.parametric_hrf_fit")
  expect_output(print(summary_fit), "Parameter Estimates")
})