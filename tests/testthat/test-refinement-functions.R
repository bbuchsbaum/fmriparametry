library(testthat)
library(fmriparametric)

test_that("refinement queue classification works correctly", {
  set.seed(123)
  n_vox <- 50
  
  # Create realistic R-squared values
  r2_voxel <- c(
    runif(15, 0.8, 0.95),    # Easy voxels (high R²)
    runif(20, 0.4, 0.7),     # Moderate voxels 
    runif(10, 0.1, 0.3),     # Hard voxels (low R²)
    runif(5, 0.0, 0.1)       # Very hard voxels
  )
  
  # Test queue classification without standard errors
  queue_result <- fmriparametric:::.classify_refinement_queue(
    r2_voxel = r2_voxel,
    se_theta_hat_voxel = NULL,
    refinement_opts = list(
      apply_refinement = TRUE,
      r2_threshold_hard = 0.3,
      r2_threshold_moderate = 0.7
    )
  )
  
  # Check structure
  expect_type(queue_result, "list")
  expect_true("queue_labels" %in% names(queue_result))
  expect_true("queue_summary" %in% names(queue_result))
  expect_true("refinement_needed" %in% names(queue_result))
  
  # Check queue assignments
  expect_equal(length(queue_result$queue_labels), n_vox)
  expect_true(queue_result$refinement_needed)
  
  # Should have some hard cases
  expect_true(any(queue_result$queue_labels == "hard_GN"))
  # Should have some moderate cases  
  expect_true(any(queue_result$queue_labels == "moderate_local_recenter"))
  # Should have some easy cases
  expect_true(any(queue_result$queue_labels == "easy"))
  
  cat("\nQueue classification test passed\n")
  print(queue_result$queue_summary)
})

test_that("refinement queue handles edge cases", {
  # Test with all good voxels
  r2_good <- rep(0.9, 10)
  result_good <- fmriparametric:::.classify_refinement_queue(
    r2_voxel = r2_good,
    refinement_opts = list(
      apply_refinement = TRUE,
      r2_threshold_hard = 0.3,
      r2_threshold_moderate = 0.7
    )
  )
  expect_true(all(result_good$queue_labels == "easy"))
  expect_false(result_good$refinement_needed)
  
  # Test with all bad voxels
  r2_bad <- rep(0.1, 10)
  result_bad <- fmriparametric:::.classify_refinement_queue(
    r2_voxel = r2_bad,
    refinement_opts = list(
      apply_refinement = TRUE,
      r2_threshold_hard = 0.3,
      r2_threshold_moderate = 0.7
    )
  )
  expect_true(all(result_bad$queue_labels == "hard_GN"))
  expect_true(result_bad$refinement_needed)
  
  # Test with refinement disabled
  result_disabled <- fmriparametric:::.classify_refinement_queue(
    r2_voxel = r2_bad,
    refinement_opts = list(apply_refinement = FALSE)
  )
  expect_true(all(result_disabled$queue_labels == "easy"))
  expect_false(result_disabled$refinement_needed)
  
  cat("\nEdge case tests passed\n")
})

test_that("moderate voxel refinement function works", {
  set.seed(456)
  n_time <- 30
  n_vox <- 5
  
  # Create test data
  Y_proj <- matrix(rnorm(n_time * n_vox, 100, 10), nrow = n_time)
  S_target_proj <- matrix(rbinom(n_time, 1, 0.2), ncol = 1)
  hrf_eval_times <- seq(0, 20, by = 0.5)
  
  # Create proper HRF interface with bounds
  base_interface <- fmriparametric:::.get_hrf_interface("lwu")
  theta_bounds <- base_interface$default_bounds()
  hrf_interface <- list(
    hrf_function = function(t, params_vector, ...) {
      base_interface$hrf_function(t, params_vector, bounds = theta_bounds, ...)
    },
    taylor_basis = function(params_vector0, t_hrf_eval, ...) {
      base_interface$taylor_basis(params_vector0, t_hrf_eval, bounds = theta_bounds, ...)
    },
    parameter_names = base_interface$parameter_names,
    default_seed = base_interface$default_seed,
    default_bounds = base_interface$default_bounds
  )
  
  # Initialize parameters
  theta_current <- matrix(rep(hrf_interface$default_seed(), n_vox),
                          nrow = n_vox, byrow = TRUE)
  r_squared <- runif(n_vox, 0.4, 0.7)  # Moderate R² values
  
  # Test moderate voxel refinement
  expect_silent(
    result_moderate <- fmriparametric:::.refine_moderate_voxels_grouped(
      voxel_idx = 1:n_vox,
      Y_proj = Y_proj,
      S_target_proj = S_target_proj,
      theta_current = theta_current,
      r_squared = r_squared,
      hrf_interface = hrf_interface,
      hrf_eval_times = hrf_eval_times,
      theta_bounds = theta_bounds
    )
  )
  
  # Check output structure
  expect_type(result_moderate, "list")
  expect_true("theta_refined" %in% names(result_moderate))
  expect_true("amplitudes" %in% names(result_moderate))
  expect_true("n_improved" %in% names(result_moderate))
  
  # Check dimensions
  expect_equal(nrow(result_moderate$theta_refined), n_vox)
  expect_equal(ncol(result_moderate$theta_refined), length(hrf_interface$parameter_names))
  expect_equal(length(result_moderate$amplitudes), n_vox)
  expect_true(is.numeric(result_moderate$n_improved))
  
  cat("\nModerate voxel refinement test passed\n")
  cat("Number improved:", result_moderate$n_improved, "\n")
})

test_that("Gauss-Newton refinement handles empty input", {
  # Test with no hard voxels
  n_vox <- 3
  n_params <- 3
  theta_hat <- matrix(runif(n_vox * n_params), nrow = n_vox)
  r2_voxel <- rep(0.8, n_vox)  # All good R²
  Y_proj <- matrix(rnorm(20 * n_vox), nrow = 20)
  S_target_proj <- matrix(rbinom(20, 1, 0.2), ncol = 1)
  scan_times <- seq(0, 38, by = 2)
  hrf_eval_times <- seq(0, 20, by = 0.5)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  theta_bounds <- hrf_interface$default_bounds()
  queue_labels <- rep("easy", n_vox)  # No hard voxels
  
  expect_silent(
    result_gn <- fmriparametric:::.gauss_newton_refinement(
      theta_hat_voxel = theta_hat,
      r2_voxel = r2_voxel,
      Y_proj = Y_proj,
      S_target_proj = S_target_proj,
      scan_times = scan_times,
      hrf_eval_times = hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_bounds = theta_bounds,
      queue_labels = queue_labels,
      max_iter_gn = 2,
      verbose = FALSE
    )
  )
  
  # Should handle empty input gracefully
  expect_type(result_gn, "list")
  expect_equal(result_gn$n_refined, 0)
  expect_equal(result_gn$n_converged, 0)
  expect_equal(result_gn$n_improved, 0)
  
  cat("\nGauss-Newton empty input test passed\n")
})

test_that("amplitude computation function works", {
  set.seed(789)
  n_time <- 25
  n_vox <- 3
  
  # Create test data
  Y_proj <- matrix(rnorm(n_time * n_vox, 50, 8), nrow = n_time)
  S_target_proj <- matrix(rbinom(n_time, 1, 0.3), ncol = 1)
  hrf_eval_times <- seq(0, 15, by = 0.5)
  
  # Create proper HRF interface with bounds
  base_interface <- fmriparametric:::.get_hrf_interface("lwu")
  theta_bounds <- base_interface$default_bounds()
  hrf_interface <- list(
    hrf_function = function(t, params_vector, ...) {
      base_interface$hrf_function(t, params_vector, bounds = theta_bounds, ...)
    },
    taylor_basis = function(params_vector0, t_hrf_eval, ...) {
      base_interface$taylor_basis(params_vector0, t_hrf_eval, bounds = theta_bounds, ...)
    },
    parameter_names = base_interface$parameter_names,
    default_seed = base_interface$default_seed,
    default_bounds = base_interface$default_bounds
  )
  
  # Parameters
  theta_hat <- matrix(rep(c(6, 3, 0.3), n_vox), nrow = n_vox, byrow = TRUE)
  
  # Test amplitude computation
  expect_silent(
    amplitudes <- fmriparametric:::.compute_amplitudes_for_voxels(
      voxel_idx = 1:n_vox,
      theta_hat = theta_hat,
      Y_proj = Y_proj,
      S_target_proj = S_target_proj,
      hrf_interface = hrf_interface,
      hrf_eval_times = hrf_eval_times
    )
  )
  
  # Check output
  expect_length(amplitudes, n_vox)
  expect_true(is.numeric(amplitudes))
  expect_false(any(is.na(amplitudes)))
  
  cat("\nAmplitude computation test passed\n")
  cat("Computed amplitudes:", round(amplitudes, 3), "\n")
})

test_that("refinement print function works", {
  # Create a mock queue result
  queue_result <- list(
    queue_summary = c(easy = 30, moderate_local_recenter = 15, hard_GN = 5),
    queue_proportions = c(easy = 0.6, moderate_local_recenter = 0.3, hard_GN = 0.1),
    queue_details = list(
      easy = list(n = 30, proportion = 0.6, r2_mean = 0.85),
      moderate_local_recenter = list(n = 15, proportion = 0.3, r2_mean = 0.55),
      hard_GN = list(n = 5, proportion = 0.1, r2_mean = 0.25)
    ),
    refinement_needed = TRUE,
    classification_criteria = list(
      r2_thresholds = c(hard = 0.3, moderate = 0.7),
      se_thresholds = NULL,
      se_available = FALSE
    )
  )
  
  # Test print function (should run without error)
  expect_output(
    fmriparametric:::.print_refinement_summary(queue_result, verbose = TRUE),
    "Refinement Queue Classification"
  )
  
  # Test silent mode
  expect_silent(
    fmriparametric:::.print_refinement_summary(queue_result, verbose = FALSE)
  )
  
  cat("\nRefinement print function test passed\n")
})
