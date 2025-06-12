library(fmriparametric)

test_that("estimate_parametric_hrf handles pathological data patterns", {
  set.seed(789)
  
  # Test 1: Perfectly periodic signal (could cause aliasing)
  n_time <- 100
  n_vox <- 3
  t <- seq(0, n_time - 1)
  fmri_periodic <- matrix(
    sin(2 * pi * t / 10) + rnorm(n_time * n_vox, sd = 0.1),
    nrow = n_time, ncol = n_vox
  )
  
  # Events at same frequency
  event_periodic <- matrix(0, nrow = n_time, ncol = 1)
  event_periodic[seq(5, n_time, by = 10), 1] <- 1
  
  result <- estimate_parametric_hrf(
    fmri_data = fmri_periodic,
    event_model = event_periodic,
    verbose = FALSE
  )
  
  expect_s3_class(result, "parametric_hrf_fit")
  expect_true(all(is.finite(result$estimated_parameters)))
  
  # Test 2: Step function data
  fmri_step <- matrix(0, nrow = n_time, ncol = 2)
  fmri_step[1:50, ] <- 10
  fmri_step[51:100, ] <- -10
  fmri_step <- fmri_step + rnorm(n_time * 2, sd = 0.5)
  
  event_step <- matrix(0, nrow = n_time, ncol = 1)
  event_step[c(25, 75), 1] <- 1
  
  result <- estimate_parametric_hrf(
    fmri_data = fmri_step,
    event_model = event_step,
    verbose = FALSE
  )
  
  expect_true(all(is.finite(result$estimated_parameters)))
  
  # Test 3: Extreme outliers
  fmri_outliers <- matrix(rnorm(n_time * 3), nrow = n_time, ncol = 3)
  # Add outliers
  fmri_outliers[c(10, 50, 90), ] <- 1000
  
  event_normal <- matrix(rbinom(n_time, 1, 0.1), ncol = 1)
  
  result <- estimate_parametric_hrf(
    fmri_data = fmri_outliers,
    event_model = event_normal,
    verbose = FALSE
  )
  
  expect_s3_class(result, "parametric_hrf_fit")
  expect_true(all(is.finite(result$estimated_parameters)))
})

test_that("estimate_parametric_hrf handles extreme parameter scenarios", {
  set.seed(456)
  n_time <- 80
  n_vox <- 5
  
  # Test 1: Very narrow bounds (forces all voxels to similar solutions)
  fmri_data <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  events <- matrix(rbinom(n_time, 1, 0.1), ncol = 1)
  
  narrow_bounds <- list(
    lower = c(5.5, 2.0, 0.3),
    upper = c(6.5, 3.0, 0.4)
  )
  
  result <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = events,
    theta_bounds = narrow_bounds,
    verbose = FALSE
  )
  
  # Check all parameters are within narrow bounds
  params <- result$estimated_parameters
  expect_true(all(params[,1] >= 5.5 & params[,1] <= 6.5))
  expect_true(all(params[,2] >= 2.0 & params[,2] <= 3.0))
  expect_true(all(params[,3] >= 0.3 & params[,3] <= 0.4))
  
  # Test 2: Seed outside bounds (should be adjusted)
  bad_seed <- c(100, 100, 100)  # Way outside bounds
  
  result <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = events,
    theta_seed = bad_seed,
    theta_bounds = fmriparametric:::.lwu_hrf_default_bounds(),
    verbose = FALSE
  )
  
  expect_s3_class(result, "parametric_hrf_fit")
  # Parameters should be within bounds despite bad seed
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  params <- result$estimated_parameters
  expect_true(all(params[,1] >= bounds$lower[1] & params[,1] <= bounds$upper[1]))
})

test_that("estimate_parametric_hrf handles data structure edge cases", {
  # Test 1: Single time point (should error)
  expect_error(
    estimate_parametric_hrf(
      fmri_data = matrix(1:5, nrow = 1),
      event_model = matrix(1, nrow = 1, ncol = 1)
    ),
    "scan_times|time points|n_time"
  )
  
  # Test 2: More event columns than voxels
  fmri_data <- matrix(rnorm(100 * 2), nrow = 100, ncol = 2)
  events_many <- matrix(rbinom(100 * 5, 1, 0.1), nrow = 100, ncol = 5)
  
  result <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = events_many,
    verbose = FALSE
  )
  
  expect_equal(nrow(result$estimated_parameters), 2)  # One per voxel
  
  # Test 3: Very long time series
  n_time_long <- 1000
  fmri_long <- matrix(rnorm(n_time_long * 2), nrow = n_time_long, ncol = 2)
  events_long <- matrix(rbinom(n_time_long, 1, 0.05), ncol = 1)
  
  result <- estimate_parametric_hrf(
    fmri_data = fmri_long,
    event_model = events_long,
    verbose = FALSE
  )
  
  expect_s3_class(result, "parametric_hrf_fit")
  expect_equal(length(result$residuals), n_time_long * 2)
})

test_that("estimate_parametric_hrf handles numerical precision limits", {
  set.seed(999)
  n_time <- 60
  
  # Test 1: Data near machine precision
  fmri_tiny <- matrix(rnorm(n_time * 2, sd = .Machine$double.eps), 
                      nrow = n_time, ncol = 2)
  events <- matrix(rbinom(n_time, 1, 0.1), ncol = 1)
  
  # Should handle gracefully with regularization
  result <- estimate_parametric_hrf(
    fmri_data = fmri_tiny,
    event_model = events,
    lambda_ridge = 0.1,  # Higher regularization for stability
    verbose = FALSE
  )
  
  expect_true(all(is.finite(result$estimated_parameters)))
  
  # Test 2: Mixed scales
  fmri_mixed <- cbind(
    rnorm(n_time, sd = 1e-10),
    rnorm(n_time, sd = 1e10)
  )
  
  result <- estimate_parametric_hrf(
    fmri_data = fmri_mixed,
    event_model = events,
    verbose = FALSE
  )
  
  expect_s3_class(result, "parametric_hrf_fit")
  expect_true(all(is.finite(result$estimated_parameters)))
  
  # Test 3: Constant columns mixed with variable
  fmri_const_mix <- cbind(
    rep(5.0, n_time),  # Constant
    rnorm(n_time),     # Variable
    rep(-3.0, n_time)  # Constant
  )
  
  result <- estimate_parametric_hrf(
    fmri_data = fmri_const_mix,
    event_model = events,
    baseline_model = "intercept",
    verbose = FALSE
  )
  
  # R² for constant voxels should reflect perfect baseline fit
  expect_true(result$r_squared[1] >= 0)  # Constant can be fit by intercept
  expect_true(result$r_squared[3] >= 0)
})

test_that("Integration with refinement handles edge cases", {
  # Enable global refinement for this test
  old_opt <- getOption("fmriparametric.refine_global")
  options(fmriparametric.refine_global = TRUE)
  on.exit(options(fmriparametric.refine_global = old_opt))
  
  set.seed(321)
  n_time <- 100
  n_vox <- 10
  
  # Create mix of easy and hard voxels
  true_params <- rbind(
    matrix(rep(c(6, 2.5, 0.35), 5), nrow = 5, byrow = TRUE),  # Easy
    matrix(rep(c(2, 0.5, 1.2), 5), nrow = 5, byrow = TRUE)    # Hard (extreme)
  )
  
  # Generate data
  fmri_data <- matrix(0, nrow = n_time, ncol = n_vox)
  events <- matrix(0, nrow = n_time, ncol = 1)
  events[seq(10, 90, by = 20), 1] <- 1
  
  hrf_times <- seq(0, 30, by = 0.5)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  
  for (v in 1:n_vox) {
    hrf <- fmriparametric:::.lwu_hrf_function(hrf_times, true_params[v,], bounds)
    conv_full <- stats::convolve(events[,1], rev(hrf), type = "open")
    fmri_data[, v] <- conv_full[1:n_time] + rnorm(n_time, sd = 0.5)
  }
  
  # Run with refinement
  result <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = events,
    refine_iters = 2,
    verbose = FALSE
  )
  
  expect_s3_class(result, "parametric_hrf_fit")
  expect_true("refinement_info" %in% names(result))
  
  # Check that hard voxels were identified
  if (!is.null(result$refinement_info$queue_labels)) {
    queue_labels <- result$refinement_info$queue_labels
    expect_true(any(queue_labels %in% c("moderate", "hard_GN")))
  }
  
  # Parameters should be reasonable
  expect_true(all(is.finite(result$estimated_parameters)))
  params <- result$estimated_parameters
  expect_true(all(params[,1] >= bounds$lower[1] & params[,1] <= bounds$upper[1]))
  expect_true(all(params[,2] >= bounds$lower[2] & params[,2] <= bounds$upper[2]))
  expect_true(all(params[,3] >= bounds$lower[3] & params[,3] <= bounds$upper[3]))
})

test_that("S3 methods handle edge cases", {
  # Create minimal fit object
  n_vox <- 2
  n_time <- 50
  
  minimal_fit <- structure(
    list(
      estimated_parameters = matrix(c(6, 2.5, 0.35, 7, 3, 0.4), 
                                   nrow = n_vox, byrow = TRUE),
      amplitudes = c(1.5, -0.5),
      parameter_names = c("tau", "sigma", "rho"),
      hrf_model = "lwu",
      r_squared = c(0.8, 0.3),
      residuals = matrix(rnorm(n_time * n_vox), nrow = n_time),
      parameter_ses = matrix(0.1, nrow = n_vox, ncol = 3),
      convergence_info = list(converged = TRUE),
      metadata = list(
        n_voxels = n_vox,
        n_timepoints = n_time,
        hrf_model = "lwu"
      )
    ),
    class = "parametric_hrf_fit"
  )
  
  # Test print with minimal object
  expect_output(print(minimal_fit), "parametric_hrf_fit")
  
  # Test coef
  coeffs <- coef(minimal_fit)
  expect_equal(dim(coeffs), c(n_vox, 3))
  
  # Test summary
  summ <- summary(minimal_fit)
  expect_s3_class(summ, "summary.parametric_hrf_fit")
  
  # Test with missing components
  incomplete_fit <- minimal_fit
  incomplete_fit$parameter_ses <- NULL
  
  # Should still work
  expect_output(print(incomplete_fit), "parametric_hrf_fit")
  summ2 <- summary(incomplete_fit)
  expect_s3_class(summ2, "summary.parametric_hrf_fit")
})