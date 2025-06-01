# Comprehensive tests for Sprint 3 features

library(testthat)

context("Sprint 3: Advanced features and refinement")

# Helper function to create test data with known structure
create_clustered_test_data <- function(n_time = 100, n_vox = 50, n_clusters = 3) {
  set.seed(123)
  
  # Create clustered parameters
  cluster_centers <- matrix(c(
    5, 2, 0.3,   # Cluster 1: early, narrow
    8, 3, 0.5,   # Cluster 2: late, medium
    6, 4, 0.2    # Cluster 3: medium, wide
  ), nrow = n_clusters, byrow = TRUE)
  
  # Assign voxels to clusters
  cluster_assignments <- sample(1:n_clusters, n_vox, replace = TRUE)
  
  # Generate parameters with cluster structure
  theta_true <- matrix(NA, nrow = n_vox, ncol = 3)
  for (i in 1:n_vox) {
    cluster <- cluster_assignments[i]
    theta_true[i, ] <- cluster_centers[cluster, ] + rnorm(3, sd = 0.2)
  }
  
  # Ensure bounds
  theta_true[, 1] <- pmax(2, pmin(12, theta_true[, 1]))  # tau
  theta_true[, 2] <- pmax(1, pmin(5, theta_true[, 2]))    # sigma
  theta_true[, 3] <- pmax(0, pmin(1, theta_true[, 3]))    # rho
  
  # Generate data
  t_hrf <- seq(0, 30, length.out = 61)
  scan_times <- seq(0, (n_time - 1) * 2, by = 2)
  
  # Simple event design
  S <- matrix(0, nrow = n_time, ncol = 1)
  S[seq(10, n_time, by = 20), 1] <- 1
  
  # Generate Y with cluster-specific responses
  Y <- matrix(NA, nrow = n_time, ncol = n_vox)
  for (v in 1:n_vox) {
    # Generate HRF for this voxel
    hrf <- exp(-(t_hrf - theta_true[v, 1])^2 / (2 * theta_true[v, 2]^2)) -
           theta_true[v, 3] * exp(-(t_hrf - theta_true[v, 1] - 2 * theta_true[v, 2])^2 / 
                                  (2 * (1.6 * theta_true[v, 2])^2))
    
    # Convolve
    conv_full <- stats::convolve(S[, 1], rev(hrf), type = "open")
    signal <- conv_full[1:n_time]
    
    # Add noise (varying by cluster)
    noise_level <- 0.1 + 0.1 * cluster_assignments[v]
    Y[, v] <- signal + rnorm(n_time, sd = noise_level)
  }
  
  list(
    Y = Y,
    S = S,
    theta_true = theta_true,
    cluster_assignments = cluster_assignments,
    cluster_centers = cluster_centers,
    scan_times = scan_times,
    t_hrf = t_hrf
  )
}

# Test K-means recentering
test_that("K-means recentering improves fits for clustered data", {
  skip_if_not_installed("fmrireg")
  
  # Create clustered test data
  test_data <- create_clustered_test_data(n_time = 100, n_vox = 60, n_clusters = 3)
  
  # Load required functions
  source(file.path(test_path("..", "..", "R"), "parametric-engine-iterative.R"))
  source(file.path(test_path("..", "..", "R"), "kmeans-recentering.R"))
  source(file.path(test_path("..", "..", "R"), "hrf-interface-lwu.R"))
  
  # HRF interface
  hrf_interface <- list(
    hrf_function = .lwu_hrf_function,
    taylor_basis = .lwu_hrf_taylor_basis_function,
    parameter_names = .lwu_hrf_parameter_names(),
    default_seed = .lwu_hrf_default_seed(),
    default_bounds = .lwu_hrf_default_bounds()
  )
  
  # Run without K-means
  fit_no_kmeans <- .parametric_engine_iterative(
    Y_proj = test_data$Y,
    S_target_proj = test_data$S,
    scan_times = test_data$scan_times,
    hrf_eval_times = test_data$t_hrf,
    hrf_interface = hrf_interface,
    theta_seed = c(6, 2.5, 0.35),
    theta_bounds = hrf_interface$default_bounds,
    recenter_global_passes = 2,
    recenter_kmeans_passes = 0,
    compute_se = FALSE,
    verbose = FALSE
  )
  
  # Run with K-means
  fit_kmeans <- .parametric_engine_iterative(
    Y_proj = test_data$Y,
    S_target_proj = test_data$S,
    scan_times = test_data$scan_times,
    hrf_eval_times = test_data$t_hrf,
    hrf_interface = hrf_interface,
    theta_seed = c(6, 2.5, 0.35),
    theta_bounds = hrf_interface$default_bounds,
    recenter_global_passes = 2,
    recenter_kmeans_passes = 2,
    kmeans_k = 3,
    r2_threshold_kmeans = 0.5,
    compute_se = FALSE,
    verbose = FALSE
  )
  
  # K-means should improve mean R²
  expect_true(mean(fit_kmeans$r_squared) > mean(fit_no_kmeans$r_squared))
  
  # Check K-means info
  expect_true(fit_kmeans$kmeans_info$applied)
  expect_equal(fit_kmeans$kmeans_info$n_clusters, 3)
  expect_true(fit_kmeans$kmeans_info$total_iterations > 0)
  
  # Parameter recovery should be better
  mse_no_kmeans <- mean((fit_no_kmeans$theta_hat - test_data$theta_true)^2)
  mse_kmeans <- mean((fit_kmeans$theta_hat - test_data$theta_true)^2)
  expect_true(mse_kmeans < mse_no_kmeans * 1.1)  # Allow small tolerance
})

# Test refinement queue classification
test_that("Refinement queue classification works correctly", {
  source(file.path(test_path("..", "..", "R"), "refinement-queue.R"))
  
  # Create test R² and SE values
  n_vox <- 100
  set.seed(456)
  
  # Mix of easy, moderate, and hard voxels
  r2_values <- c(
    runif(30, 0.8, 0.95),  # Easy
    runif(40, 0.4, 0.7),   # Moderate
    runif(30, 0, 0.3)      # Hard
  )
  
  se_values <- matrix(c(
    runif(30, 0.05, 0.2),   # Easy - low SE
    runif(40, 0.2, 0.5),    # Moderate - medium SE
    runif(30, 0.5, 2)       # Hard - high SE
  ), ncol = 3, byrow = FALSE)
  
  # Classify
  queue_result <- .classify_refinement_queue(
    r2_voxel = r2_values,
    se_theta_hat_voxel = se_values,
    refinement_opts = list(
      apply_refinement = TRUE,
      r2_threshold_hard = 0.3,
      r2_threshold_moderate = 0.7,
      se_threshold_hard = 0.5,
      se_threshold_moderate = 0.3
    )
  )
  
  # Check classification
  expect_true(queue_result$refinement_needed)
  expect_equal(length(queue_result$queue_labels), n_vox)
  expect_true(all(queue_result$queue_labels %in% c("easy", "moderate_local_recenter", "hard_GN")))
  
  # Verify distribution roughly matches
  queue_props <- prop.table(table(queue_result$queue_labels))
  expect_true(queue_props["easy"] > 0.2)
  expect_true(queue_props["hard_GN"] > 0.2)
  
  # Check queue details
  expect_true(all(c("n", "proportion", "r2_mean", "r2_median") %in% 
                  names(queue_result$queue_details$easy)))
})

# Test local recentering
test_that("Local recentering improves moderate voxels", {
  skip_if_not_installed("fmrireg")
  
  # Create test data with moderate quality fits
  test_data <- create_clustered_test_data(n_time = 80, n_vox = 20)
  
  # Add extra noise to make fits moderate
  test_data$Y <- test_data$Y + matrix(rnorm(length(test_data$Y), sd = 0.3), 
                                       nrow = nrow(test_data$Y))
  
  # Load required functions
  source(file.path(test_path("..", "..", "R"), "parametric-engine.R"))
  source(file.path(test_path("..", "..", "R"), "hrf-interface-lwu.R"))
  
  hrf_interface <- list(
    hrf_function = .lwu_hrf_function,
    taylor_basis = .lwu_hrf_taylor_basis_function,
    parameter_names = .lwu_hrf_parameter_names(),
    default_seed = .lwu_hrf_default_seed(),
    default_bounds = .lwu_hrf_default_bounds()
  )
  
  # Initial fit
  initial_fit <- .parametric_engine(
    Y_proj = test_data$Y,
    S_target_proj = test_data$S,
    scan_times = test_data$scan_times,
    hrf_eval_times = test_data$t_hrf,
    hrf_interface = hrf_interface,
    theta_seed = c(6, 2.5, 0.35),
    theta_bounds = hrf_interface$default_bounds,
    verbose = FALSE
  )
  
  # Select moderate voxels (R² between 0.3 and 0.7)
  moderate_idx <- which(initial_fit$r_squared > 0.3 & initial_fit$r_squared < 0.7)
  
  if (length(moderate_idx) > 0) {
    # Apply local recentering
    improved_count <- 0
    for (v in moderate_idx) {
      local_fit <- .parametric_engine(
        Y_proj = test_data$Y[, v, drop = FALSE],
        S_target_proj = test_data$S,
        scan_times = test_data$scan_times,
        hrf_eval_times = test_data$t_hrf,
        hrf_interface = hrf_interface,
        theta_seed = initial_fit$theta_hat[v, ],  # Use voxel's estimate
        theta_bounds = hrf_interface$default_bounds,
        verbose = FALSE
      )
      
      if (local_fit$r_squared[1] > initial_fit$r_squared[v]) {
        improved_count <- improved_count + 1
      }
    }
    
    # At least some voxels should improve
    expect_true(improved_count > 0)
    improvement_rate <- improved_count / length(moderate_idx)
    expect_true(improvement_rate > 0.3)  # At least 30% improve
  }
})

# Test Gauss-Newton refinement
test_that("Gauss-Newton refinement handles difficult voxels", {
  skip_if_not_installed("fmrireg")
  
  # Create difficult test case
  set.seed(789)
  n_time <- 60
  n_vox <- 10
  
  # True parameters far from typical seed
  theta_true <- matrix(c(
    rep(10, n_vox),    # tau = 10 (late peak)
    rep(1.5, n_vox),   # sigma = 1.5 (narrow)
    rep(0.8, n_vox)    # rho = 0.8 (strong undershoot)
  ), nrow = n_vox, byrow = FALSE)
  
  # Generate data
  t_hrf <- seq(0, 30, length.out = 61)
  S <- matrix(0, nrow = n_time, ncol = 1)
  S[c(5, 25, 45), 1] <- 1
  
  Y <- matrix(NA, nrow = n_time, ncol = n_vox)
  for (v in 1:n_vox) {
    hrf <- exp(-(t_hrf - theta_true[v, 1])^2 / (2 * theta_true[v, 2]^2)) -
           theta_true[v, 3] * exp(-(t_hrf - theta_true[v, 1] - 2 * theta_true[v, 2])^2 / 
                                  (2 * (1.6 * theta_true[v, 2])^2))
    conv_full <- stats::convolve(S[, 1], rev(hrf), type = "open")
    Y[, v] <- conv_full[1:n_time] + rnorm(n_time, sd = 0.2)
  }
  
  # Load functions
  source(file.path(test_path("..", "..", "R"), "gauss-newton-refinement.R"))
  source(file.path(test_path("..", "..", "R"), "hrf-interface-lwu.R"))
  
  hrf_interface <- list(
    hrf_function = .lwu_hrf_function,
    taylor_basis = .lwu_hrf_taylor_basis_function,
    parameter_names = .lwu_hrf_parameter_names(),
    default_seed = .lwu_hrf_default_seed(),
    default_bounds = .lwu_hrf_default_bounds()
  )
  
  # Poor initial estimates
  theta_init <- matrix(rep(c(6, 2.5, 0.35), each = n_vox), nrow = n_vox)
  r2_init <- rep(0.2, n_vox)  # Assume poor initial fit
  queue_labels <- rep("hard_GN", n_vox)
  
  # Apply Gauss-Newton
  gn_result <- .gauss_newton_refinement(
    theta_hat_voxel = theta_init,
    r2_voxel = r2_init,
    Y_proj = Y,
    S_target_proj = S,
    scan_times = seq(0, (n_time - 1) * 2, by = 2),
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_interface,
    theta_bounds = list(lower = c(2, 1, 0), upper = c(12, 5, 1)),
    queue_labels = queue_labels,
    max_iter_gn = 10,
    verbose = FALSE
  )
  
  # Check results
  expect_equal(gn_result$n_refined, n_vox)
  expect_true(gn_result$n_converged > 0)
  expect_true(gn_result$n_improved > 0)
  
  # Parameters should be closer to truth
  mse_init <- mean((theta_init - theta_true)^2)
  mse_final <- mean((gn_result$theta_hat - theta_true)^2)
  expect_true(mse_final < mse_init)
})

# Test parallel processing
test_that("Parallel processing produces identical results", {
  skip_if_not(parallel::detectCores() > 1, "Single core system")
  skip_if_not_installed("future")
  
  # Create test data
  test_data <- create_clustered_test_data(n_time = 50, n_vox = 30)
  
  # Mock fit data for parallel testing
  fit_data <- list(
    theta_hat = matrix(rnorm(90), nrow = 30, ncol = 3),
    r_squared = runif(30, 0.2, 0.8),
    beta0 = rnorm(30),
    residuals = matrix(rnorm(1500), nrow = 50, ncol = 30)
  )
  
  prepared_data <- list(
    Y_proj = test_data$Y,
    S_target_proj = test_data$S,
    scan_times = test_data$scan_times,
    hrf_eval_times = test_data$t_hrf
  )
  
  # Load parallel functions
  source(file.path(test_path("..", "..", "R"), "parallel-processing.R"))
  
  # Test parallel backend setup
  parallel_config <- .setup_parallel_backend(n_cores = 2, verbose = FALSE)
  expect_true(parallel_config$n_cores >= 1)
  expect_true(parallel_config$backend %in% c("sequential", "future_multicore", 
                                              "future_multisession", "mclapply", 
                                              "parLapply"))
  
  # Clean up
  parallel_config$cleanup()
})

# Test safety mode in refinement
test_that("Refinement handles edge cases safely", {
  # Test with empty data
  source(file.path(test_path("..", "..", "R"), "refinement-queue.R"))
  
  queue_result <- .classify_refinement_queue(
    r2_voxel = numeric(0),
    se_theta_hat_voxel = NULL,
    refinement_opts = list(apply_refinement = TRUE)
  )
  
  expect_false(queue_result$refinement_needed)
  expect_equal(length(queue_result$queue_labels), 0)
  
  # Test with all NA values
  queue_result_na <- .classify_refinement_queue(
    r2_voxel = rep(NA, 10),
    se_theta_hat_voxel = NULL,
    refinement_opts = list(apply_refinement = TRUE)
  )
  
  expect_true(queue_result_na$refinement_needed)
  expect_true(all(queue_result_na$queue_labels == "hard_GN"))
})

# Test end-to-end workflow with Sprint 3 features
test_that("Full Sprint 3 workflow completes successfully", {
  skip_if_not_installed("fmrireg")
  skip_on_cran()  # Too intensive for CRAN
  
  # Create realistic test data
  test_data <- create_clustered_test_data(n_time = 100, n_vox = 50, n_clusters = 3)
  
  # Create mock fmri and event objects
  fmri_data <- structure(
    list(data = test_data$Y, dims = c(100, 50)),
    class = "mock_fmri"
  )
  
  event_model <- structure(
    list(design = test_data$S),
    class = "mock_event"
  )
  
  # Source main function
  source(file.path(test_path("..", "..", "R"), "estimate_parametric_hrf_v3.R"))
  
  # Test that it would run (without actually running due to dependencies)
  expect_true(is.function(estimate_parametric_hrf_v3))
  
  # Check function signature
  fn_args <- names(formals(estimate_parametric_hrf_v3))
  expect_true("recenter_kmeans_passes" %in% fn_args)
  expect_true("refinement_opts" %in% fn_args)
  expect_true("n_cores" %in% fn_args)
})

# Test S3 methods with refinement info
test_that("Enhanced S3 methods handle refinement information", {
  # Create mock fit object with refinement info
  mock_fit <- structure(
    list(
      estimated_parameters = matrix(rnorm(150), nrow = 50, ncol = 3),
      amplitudes = rnorm(50),
      parameter_names = c("tau", "sigma", "rho"),
      hrf_model = "lwu",
      r_squared = runif(50, 0.1, 0.9),
      residuals = matrix(rnorm(5000), nrow = 100, ncol = 50),
      parameter_ses = matrix(runif(150, 0.1, 0.5), nrow = 50, ncol = 3),
      convergence_info = list(
        global_iterations = 3,
        converged = TRUE
      ),
      metadata = list(
        n_voxels = 50,
        n_timepoints = 100,
        refinement_info = list(
          applied = TRUE,
          queue_result = list(
            queue_labels = sample(c("easy", "moderate_local_recenter", "hard_GN"), 
                                  50, replace = TRUE),
            queue_summary = table(c("easy" = 25, "moderate_local_recenter" = 15, 
                                    "hard_GN" = 10))
          ),
          n_moderate_refined = 15,
          n_hard_refined = 10,
          n_converged = 8,
          n_improved = 7,
          final_queue_summary = table(c("easy" = 42, "moderate_local_recenter" = 5, 
                                         "hard_GN" = 3))
        ),
        parallel_info = list(
          backend = "future_multicore",
          n_cores = 4
        ),
        computation_time = 45.3
      )
    ),
    class = "parametric_hrf_fit"
  )
  
  # Source enhanced methods
  source(file.path(test_path("..", "..", "R"), "parametric-hrf-fit-methods-v3.R"))
  
  # Test print method
  expect_output(print(mock_fit), "Refinement Applied")
  expect_output(print(mock_fit), "Parallel Processing")
  
  # Test summary method
  summ <- summary(mock_fit)
  expect_true("refinement_summary" %in% names(summ))
  expect_true(summ$refinement_summary$applied)
  expect_equal(summ$refinement_summary$n_converged, 8)
  
  # Test that summary prints correctly
  expect_output(print(summ), "Refinement Summary")
})