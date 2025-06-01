# Comprehensive workflow validation tests
# These tests validate the complete end-to-end functionality

test_that("Basic workflow with simulated data works", {
  # Load the package to ensure functions are available
  library(fmriparametric)
  set.seed(42)
  
  # Generate realistic synthetic data
  n_time <- 200
  n_vox <- 10
  tr <- 2.0
  
  # True HRF parameters (LWU model)
  true_params <- c(tau = 6, sigma = 2.5, rho = 0.35)
  
  # Generate event times
  scan_times <- seq(0, (n_time - 1) * tr, by = tr)
  event_times <- c(20, 60, 120, 180)  # Events at these times
  
  # Create event model matrix
  event_model <- matrix(0, nrow = n_time, ncol = 1)
  event_indices <- round(event_times / tr) + 1
  event_indices <- event_indices[event_indices <= n_time]
  event_model[event_indices, 1] <- 1
  
  # Generate true HRF
  hrf_times <- seq(0, 30, by = tr)
  true_hrf <- fmrireg::hrf_lwu(hrf_times, 
                               tau = true_params[1], 
                               sigma = true_params[2], 
                               rho = true_params[3])
  
  # Convolve with events to get expected BOLD
  expected_bold <- rep(0, n_time)
  for (i in seq_along(event_indices)) {
    start_idx <- event_indices[i]
    end_idx <- min(start_idx + length(true_hrf) - 1, n_time)
    response_length <- end_idx - start_idx + 1
    expected_bold[start_idx:end_idx] <- expected_bold[start_idx:end_idx] + 
      2.0 * true_hrf[1:response_length]  # amplitude = 2.0
  }
  
  # Add noise and create multi-voxel data
  noise_level <- 0.1
  fmri_data <- matrix(expected_bold, nrow = n_time, ncol = n_vox) + 
    matrix(rnorm(n_time * n_vox, sd = noise_level), nrow = n_time, ncol = n_vox)
  
  # Test 1: Basic estimation
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    hrf_eval_times = hrf_times,
    verbose = FALSE
  )
  
  # Validate structure
  expect_s3_class(fit, "parametric_hrf_fit")
  expect_true(is.matrix(coef(fit)))
  expect_equal(nrow(coef(fit)), n_vox)
  expect_equal(ncol(coef(fit)), 3)  # tau, sigma, rho
  expect_equal(colnames(coef(fit)), c("tau", "sigma", "rho"))
  
  # Validate parameter recovery (should be close to true values)
  mean_params <- colMeans(coef(fit))
  expect_true(abs(mean_params["tau"] - true_params["tau"]) < 1.0)
  expect_true(abs(mean_params["sigma"] - true_params["sigma"]) < 0.5)
  expect_true(abs(mean_params["rho"] - true_params["rho"]) < 0.2)
  
  # Test 2: S3 methods work
  expect_output(print(fit), "Parametric HRF Fit")
  
  summary_fit <- summary(fit)
  expect_s3_class(summary_fit, "summary.parametric_hrf_fit")
  expect_output(print(summary_fit), "Parameter Statistics")
  
  # Test 3: Standard errors are reasonable
  if (!is.null(fit$standard_errors)) {
    expect_true(all(is.finite(fit$standard_errors)))
    expect_true(all(fit$standard_errors > 0))
  }
})

test_that("Workflow with fmrireg objects works", {
  skip_if_not_installed("fmrireg")

  set.seed(123)

  # Create small fmrireg dataset and event model
  n_time <- 100
  n_vox <- 5

  bold_mat <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  fmri_ds <- fmrireg::matrix_dataset(bold_mat, TR = 2, run_length = n_time)

  ev <- fmrireg::event_model(
    onset = c(20, 40, 60, 80),
    blockids = rep(1, 4),
    durations = rep(0, 4)
  )

  fit <- estimate_parametric_hrf(
    fmri_data = fmri_ds,
    event_model = ev,
    parametric_hrf = "lwu",
    verbose = FALSE,
    compute_se = FALSE
  )

  expect_s3_class(fit, "parametric_hrf_fit")
  expect_equal(nrow(coef(fit)), n_vox)
})

test_that("Parameter bounds are enforced", {
  set.seed(456)
  
  # Create simple test data
  n_time <- 50
  n_vox <- 3
  
  fmri_data <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  event_model <- matrix(rbinom(n_time, 1, 0.2), ncol = 1)
  
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  params <- coef(fit)
  
  # Check LWU bounds are enforced
  expect_true(all(params[, "tau"] >= 0))
  expect_true(all(params[, "tau"] <= 20))
  expect_true(all(params[, "sigma"] >= 0.05))
  expect_true(all(params[, "sigma"] <= 10))
  expect_true(all(params[, "rho"] >= 0))
  expect_true(all(params[, "rho"] <= 1.5))
})

test_that("Error handling for invalid inputs", {
  # Test invalid fmri_data
  expect_error(
    estimate_parametric_hrf(
      fmri_data = "not_a_matrix",
      event_model = matrix(1, 10, 1),
      parametric_hrf = "lwu"
    ),
    "fmri_data.*matrix"
  )
  
  # Test dimension mismatch
  expect_error(
    estimate_parametric_hrf(
      fmri_data = matrix(1, 10, 5),
      event_model = matrix(1, 20, 1),  # Wrong number of rows
      parametric_hrf = "lwu"
    ),
    "dimension"
  )
  
  # Test invalid HRF model
  expect_error(
    estimate_parametric_hrf(
      fmri_data = matrix(1, 10, 5),
      event_model = matrix(1, 10, 1),
      parametric_hrf = "invalid_model"
    ),
    "Unsupported.*HRF.*model"
  )
})

test_that("Performance is reasonable for medium datasets", {
  skip_on_cran()  # Skip performance tests on CRAN
  
  set.seed(789)
  
  # Medium-sized dataset
  n_time <- 500   # ~16 minutes of data at TR=2s
  n_vox <- 100    # 100 voxels
  
  fmri_data <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  event_model <- matrix(rbinom(n_time, 1, 0.05), ncol = 1)
  
  # Should complete within reasonable time
  start_time <- Sys.time()
  
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    verbose = FALSE
  )
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Should complete within 30 seconds for this size
  expect_true(elapsed < 30, 
              info = sprintf("Took %.2f seconds for %d timepoints x %d voxels", 
                           elapsed, n_time, n_vox))
  
  # Results should be valid
  expect_s3_class(fit, "parametric_hrf_fit")
  expect_equal(dim(coef(fit)), c(n_vox, 3))
})