library(fmriparametric)

test_that(".validate_fmri_data handles edge cases", {
  # Test 1: NULL input
  expect_error(
    fmriparametric:::.validate_fmri_data(NULL),
    "fmri_data cannot be NULL"
  )
  
  # Test 2: Empty matrix
  expect_error(
    fmriparametric:::.validate_fmri_data(matrix(nrow = 0, ncol = 0)),
    "Input arrays cannot be empty"
  )
  
  # Test 3: Too few time points
  expect_error(
    fmriparametric:::.validate_fmri_data(matrix(1:5, nrow = 5, ncol = 1)),
    "Insufficient time points"
  )
  
  # Test 4: No voxels
  expect_error(
    fmriparametric:::.validate_fmri_data(matrix(1:20, nrow = 20, ncol = 0)),
    "No voxels found"
  )
  
  # Test 5: Non-numeric data
  expect_error(
    fmriparametric:::.validate_fmri_data(matrix(letters[1:20], nrow = 10, ncol = 2)),
    "must contain numeric values"
  )
  
  # Test 6: Data with NAs (> 50%)
  data_mostly_na <- matrix(NA, nrow = 20, ncol = 5)
  data_mostly_na[1:5, 1] <- rnorm(5)  # Only 5% non-NA
  expect_error(
    fmriparametric:::.validate_fmri_data(data_mostly_na),
    "More than 50% of data is NA"
  )
  
  # Test 7: Data with some NAs (10-50%)
  data_some_na <- matrix(rnorm(100), nrow = 20, ncol = 5)
  data_some_na[1:30] <- NA  # 30% NA
  expect_warning(
    result <- fmriparametric:::.validate_fmri_data(data_some_na),
    "30% of data contains NA values"
  )
  expect_equal(result$n_na, 30)
  
  # Test 8: Data with infinite values
  data_inf <- matrix(rnorm(100), nrow = 20, ncol = 5)
  data_inf[10, 2] <- Inf
  data_inf[15, 3] <- -Inf
  expect_error(
    fmriparametric:::.validate_fmri_data(data_inf),
    "Data contains 2 infinite values"
  )
  
  # Test 9: Constant voxels
  data_const <- matrix(rnorm(60), nrow = 20, ncol = 3)
  data_const[, 2] <- 5  # Constant column
  expect_warning(
    result <- fmriparametric:::.validate_fmri_data(data_const),
    "Found 1 constant voxels"
  )
  expect_equal(result$n_constant, 1)
  
  # Test 10: All-zero voxels
  data_zero <- matrix(rnorm(60), nrow = 20, ncol = 3)
  data_zero[, 3] <- 0
  expect_warning(
    result <- fmriparametric:::.validate_fmri_data(data_zero),
    "Found 1 all-zero voxels"
  )
  expect_equal(result$n_zero, 1)
  
  # Test 11: Data frame input
  df_input <- data.frame(V1 = rnorm(20), V2 = rnorm(20))
  expect_warning(
    result <- fmriparametric:::.validate_fmri_data(df_input),
    "Converting data.frame to matrix"
  )
  expect_true(is.matrix(result$data))
  
  # Test 12: Mock fmri_dataset object
  mock_dataset <- structure(
    list(data = matrix(rnorm(100), nrow = 20, ncol = 5)),
    class = "fmri_dataset"
  )
  result <- fmriparametric:::.validate_fmri_data(mock_dataset)
  expect_equal(dim(result$data), c(20, 5))
  expect_equal(result$type, "fmri_dataset")
})

test_that(".validate_event_model handles edge cases", {
  n_time <- 50
  
  # Test 1: NULL input
  expect_error(
    fmriparametric:::.validate_event_model(NULL, n_time),
    "event_model cannot be NULL"
  )
  
  # Test 2: Dimension mismatch
  events_wrong <- matrix(1:40, nrow = 40, ncol = 1)
  expect_error(
    fmriparametric:::.validate_event_model(events_wrong, n_time),
    "don't match fmri_data time points"
  )
  
  # Test 3: No events (all zeros)
  events_empty <- matrix(0, nrow = n_time, ncol = 1)
  expect_warning(
    result <- fmriparametric:::.validate_event_model(events_empty, n_time),
    "No events detected"
  )
  expect_equal(result$n_events, 0)
  
  # Test 4: Very high event density
  events_dense <- matrix(1, nrow = n_time, ncol = 1)  # All ones
  expect_warning(
    result <- fmriparametric:::.validate_event_model(events_dense, n_time),
    "Very high event density"
  )
  expect_equal(result$event_density, 1.0)
  
  # Test 5: Very low event density
  events_sparse <- matrix(0, nrow = n_time, ncol = 1)
  events_sparse[25] <- 1  # Single event
  expect_warning(
    result <- fmriparametric:::.validate_event_model(events_sparse, n_time),
    "Very low event density"
  )
  expect_equal(result$n_events, 1)
  
  # Test 6: Mock event_model object
  mock_event <- structure(
    list(terms = list(matrix(rbinom(n_time, 1, 0.1), ncol = 1))),
    class = "event_model"
  )
  result <- fmriparametric:::.validate_event_model(mock_event, n_time)
  expect_equal(nrow(result$design), n_time)
  expect_equal(result$type, "event_model")
  
  # Test 7: Numeric vector (should be coerced to matrix)
  events_vec <- rbinom(n_time, 1, 0.2)
  result <- fmriparametric:::.validate_event_model(events_vec, n_time)
  expect_true(is.matrix(result$design))
  expect_equal(dim(result$design), c(n_time, 1))
})

test_that(".validate_theta_bounds handles edge cases", {
  n_params <- 3
  param_names <- c("tau", "sigma", "rho")
  
  # Test 1: NULL is valid (use defaults)
  result <- fmriparametric:::.validate_theta_bounds(NULL, n_params)
  expect_null(result)
  
  # Test 2: Not a list
  expect_error(
    fmriparametric:::.validate_theta_bounds(c(1, 2, 3), n_params),
    "theta_bounds must be a list"
  )
  
  # Test 3: Missing required elements
  bounds_incomplete <- list(lower = c(0, 0.1, 0))
  expect_error(
    fmriparametric:::.validate_theta_bounds(bounds_incomplete, n_params),
    "missing required elements: upper"
  )
  
  # Test 4: Wrong dimensions
  bounds_wrong_dim <- list(
    lower = c(0, 0.1),  # Only 2 values
    upper = c(10, 5, 1)  # 3 values
  )
  expect_error(
    fmriparametric:::.validate_theta_bounds(bounds_wrong_dim, n_params),
    "dimensions incorrect"
  )
  
  # Test 5: Non-numeric bounds
  bounds_non_numeric <- list(
    lower = c("0", "0.1", "0"),
    upper = c(10, 5, 1)
  )
  expect_error(
    fmriparametric:::.validate_theta_bounds(bounds_non_numeric, n_params),
    "must contain numeric values"
  )
  
  # Test 6: Lower >= upper
  bounds_invalid <- list(
    lower = c(0, 5, 0),
    upper = c(10, 2, 1)  # upper[2] < lower[2]
  )
  expect_error(
    fmriparametric:::.validate_theta_bounds(bounds_invalid, n_params, param_names),
    "lower must be less than upper.*sigma"
  )
  
  # Test 7: Non-physiological LWU bounds (warnings)
  bounds_weird <- list(
    lower = c(-5, 0.05, -1),  # Negative tau and rho
    upper = c(50, 30, 3)      # Very large values
  )
  expect_warning(
    result <- fmriparametric:::.validate_theta_bounds(bounds_weird, n_params, param_names),
    "tau lower bound < 0"
  )
  
  # Test 8: Valid bounds
  bounds_valid <- list(
    lower = c(0, 0.1, 0),
    upper = c(20, 10, 1.5)
  )
  result <- fmriparametric:::.validate_theta_bounds(bounds_valid, n_params, param_names)
  expect_equal(result, bounds_valid)
})

test_that(".validate_numeric_param handles edge cases", {
  # Test 1: NULL with allow_null = TRUE
  result <- fmriparametric:::.validate_numeric_param(
    NULL, "test_param", allow_null = TRUE, default = 5
  )
  expect_equal(result, 5)
  
  # Test 2: NULL with allow_null = FALSE
  expect_error(
    fmriparametric:::.validate_numeric_param(
      NULL, "test_param", allow_null = FALSE
    ),
    "test_param cannot be NULL"
  )
  
  # Test 3: Non-numeric input
  expect_error(
    fmriparametric:::.validate_numeric_param("5", "test_param"),
    "must be a single numeric value"
  )
  
  # Test 4: Vector input
  expect_error(
    fmriparametric:::.validate_numeric_param(c(1, 2, 3), "test_param"),
    "must be a single numeric value.*length 3"
  )
  
  # Test 5: NA value
  expect_error(
    fmriparametric:::.validate_numeric_param(NA, "test_param"),
    "test_param cannot be NA"
  )
  
  # Test 6: Infinite value
  expect_error(
    fmriparametric:::.validate_numeric_param(Inf, "test_param"),
    "test_param cannot be infinite"
  )
  
  # Test 7: Value outside range
  expect_error(
    fmriparametric:::.validate_numeric_param(
      -5, "test_param", min_val = 0, max_val = 10
    ),
    "outside valid range.*0.*10"
  )
  
  # Test 8: Valid value at boundary
  result <- fmriparametric:::.validate_numeric_param(
    0, "test_param", min_val = 0, max_val = 10
  )
  expect_equal(result, 0)
  
  # Test 9: Valid value
  result <- fmriparametric:::.validate_numeric_param(
    5.5, "test_param", min_val = 0, max_val = 10
  )
  expect_equal(result, 5.5)
})

test_that(".rock_solid_validate_inputs performs comprehensive validation", {
  # Create valid test data
  fmri_data <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  event_model <- matrix(rbinom(100, 1, 0.1), ncol = 1)
  
  # Test 1: All valid inputs
  result <- fmriparametric:::.rock_solid_validate_inputs(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "lwu",
    theta_seed = c(6, 2.5, 0.35),
    theta_bounds = list(lower = c(0, 0.1, 0), upper = c(20, 10, 1.5)),
    hrf_span = 30,
    lambda_ridge = 0.01,
    recenter_global_passes = 3,
    recenter_epsilon = 0.01,
    r2_threshold = 0.1,
    mask = NULL,
    verbose = FALSE
  )
  
  expect_type(result, "list")
  expect_equal(result$fmri_data$n_time, 100)
  expect_equal(result$fmri_data$n_vox, 5)
  expect_equal(result$parametric_hrf, "lwu")
  
  # Test 2: Invalid HRF model
  expect_error(
    fmriparametric:::.rock_solid_validate_inputs(
      fmri_data = fmri_data,
      event_model = event_model,
      parametric_hrf = "invalid_model",
      theta_seed = c(6, 2.5, 0.35),
      theta_bounds = NULL,
      hrf_span = 30,
      lambda_ridge = 0.01,
      recenter_global_passes = 3,
      recenter_epsilon = 0.01,
      r2_threshold = 0.1,
      mask = NULL,
      verbose = FALSE
    ),
    "Only 'lwu' model currently supported"
  )
  
  # Test 3: Parameter out of range
  expect_error(
    fmriparametric:::.rock_solid_validate_inputs(
      fmri_data = fmri_data,
      event_model = event_model,
      parametric_hrf = "lwu",
      theta_seed = c(6, 2.5, 0.35),
      theta_bounds = NULL,
      hrf_span = 100,  # Too large
      lambda_ridge = 0.01,
      recenter_global_passes = 3,
      recenter_epsilon = 0.01,
      r2_threshold = 0.1,
      mask = NULL,
      verbose = FALSE
    ),
    "outside valid range.*5.*60"
  )
  
  # Test 4: Defaults are applied
  result <- fmriparametric:::.rock_solid_validate_inputs(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = "LWU",  # Should be lowercased
    theta_seed = NULL,
    theta_bounds = NULL,
    hrf_span = NULL,  # Should use default
    lambda_ridge = NULL,
    recenter_global_passes = NULL,
    recenter_epsilon = NULL,
    r2_threshold = NULL,
    mask = NULL,
    verbose = TRUE
  )
  
  expect_equal(result$parametric_hrf, "lwu")
  expect_equal(result$hrf_span, 30)
  expect_equal(result$lambda_ridge, 0.01)
  expect_equal(result$recenter_global_passes, 3)
  expect_true(result$verbose)
})