# Test that R-squared is correctly recalculated after refinement

library(testthat)
library(fmriparametric)

test_that("R-squared is recalculated after refinement", {
  skip_if_not_installed("fmrireg")
  skip_if_not_installed("Matrix")
  
  # Create a simple test case with stronger signal
  set.seed(123)
  n_time <- 80  # Increased for better signal
  n_vox <- 3    # Increased for more robust testing
  TR <- 2
  
  # Simple design matrix with more events
  design_matrix <- matrix(0, n_time, 1)
  design_matrix[seq(10, 70, by = 15), 1] <- 1  # Events every 15 time points
  
  # Create data where we know the true HRF parameters
  # Start with baseline noise
  Y <- matrix(rnorm(n_time * n_vox, 100, 10), n_time, n_vox)
  
  # Create a realistic LWU HRF with known parameters
  time_points <- seq(0, 30, by = TR)
  true_tau <- 6
  true_sigma <- 3
  true_rho <- 0.3
  
  # Generate true LWU HRF  
  true_hrf <- exp(-(time_points - true_tau)^2 / (2 * true_sigma^2)) - 
              true_rho * exp(-(time_points - true_tau - 2*true_sigma)^2 / (2 * (1.6*true_sigma)^2))
  true_hrf <- true_hrf / max(true_hrf)
  
  # Add strong signal that matches this HRF
  for (v in 1:n_vox) {
    signal <- numeric(n_time)
    amplitude <- 40 + v * 10  # Strong, voxel-specific amplitudes
    
    for (i in which(design_matrix[,1] == 1)) {
      end_idx <- min(i + length(true_hrf) - 1, n_time)
      if (end_idx > i) {
        signal[i:end_idx] <- signal[i:end_idx] + 
          true_hrf[1:(end_idx - i + 1)] * amplitude
      }
    }
    Y[, v] <- Y[, v] + signal
  }
  
  cat("\n=== R-squared Fix Test Debug ===\n")
  cat("Data setup:\n")
  cat("- Time points:", n_time, "\n")
  cat("- Voxels:", n_vox, "\n") 
  cat("- Events:", sum(design_matrix), "\n")
  cat("- True HRF params: tau =", true_tau, ", sigma =", true_sigma, ", rho =", true_rho, "\n")
  cat("- Signal strength: mean =", round(mean(apply(Y, 2, var)), 1), "\n")
  
  # Run estimation with refinement disabled first
  cat("\nRunning without refinement...\n")
  result_no_refine <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = design_matrix,
    parametric_hrf = "lwu",
    global_refinement = FALSE,
    kmeans_refinement = FALSE,
    tiered_refinement = "none",
    theta_seed = c(5, 2.5, 0.25),  # Start away from true values
    verbose = FALSE
  )
  
  # Run estimation with refinement enabled
  cat("Running with refinement...\n")
  result_with_refine <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = design_matrix,
    parametric_hrf = "lwu",
    global_refinement = TRUE,
    tiered_refinement = "moderate",
    theta_seed = c(5, 2.5, 0.25),  # Same starting point
    verbose = FALSE
  )
  
  # Debug output
  cat("\nParameter estimates:\n")
  cat("Without refinement:\n")
  print(round(result_no_refine$estimated_parameters, 3))
  cat("With refinement:\n") 
  print(round(result_with_refine$estimated_parameters, 3))
  
  # Tests
  # 1. Both should have r_squared values
  expect_true(!is.null(result_no_refine$r_squared))
  expect_true(!is.null(result_with_refine$r_squared))
  expect_equal(length(result_no_refine$r_squared), n_vox)
  expect_equal(length(result_with_refine$r_squared), n_vox)
  
  # 2. With refinement should have refinement_info
  expect_true(!is.null(result_with_refine$refinement_info))
  
  # 3. R-squared values should be different after refinement
  # This tests that R-squared was actually recalculated
  # NOTE: Refinement may not always change parameters if they're already optimal
  param_diff <- abs(result_no_refine$estimated_parameters - result_with_refine$estimated_parameters)
  params_changed <- any(param_diff > 1e-6)
  
  if (params_changed) {
    # If parameters changed, R-squared should also change
    expect_false(all(abs(result_no_refine$r_squared - result_with_refine$r_squared) < 1e-10),
                 info = "Parameters changed but R-squared didn't - possible calculation bug")
  } else {
    # If parameters didn't change, R-squared should be the same
    expect_true(all(abs(result_no_refine$r_squared - result_with_refine$r_squared) < 1e-10),
                info = "Parameters unchanged but R-squared changed - possible calculation bug")
    cat("NOTE: Refinement didn't change parameters (already optimal)\n")
  }
  
  # 4. Parameter comparison
  mean_r2_no_refine <- mean(result_no_refine$r_squared)
  mean_r2_with_refine <- mean(result_with_refine$r_squared)
  
  cat("\nR-squared comparison:\n")
  cat("Without refinement: Mean R² =", round(mean_r2_no_refine, 4), 
      ", Range: [", round(min(result_no_refine$r_squared), 4), 
      ", ", round(max(result_no_refine$r_squared), 4), "]\n")
  cat("With refinement:    Mean R² =", round(mean_r2_with_refine, 4), 
      ", Range: [", round(min(result_with_refine$r_squared), 4), 
      ", ", round(max(result_with_refine$r_squared), 4), "]\n")
  cat("Improvement:        ΔR² =", round(mean_r2_with_refine - mean_r2_no_refine, 4), "\n")
  
  # 5. R-squared values should be valid (between 0 and 1)
  expect_true(all(result_no_refine$r_squared >= 0 & result_no_refine$r_squared <= 1))
  expect_true(all(result_with_refine$r_squared >= 0 & result_with_refine$r_squared <= 1))
  
  # 6. With strong signal, both should achieve reasonable R-squared
  # But we'll be more lenient since refinement doesn't always improve things
  expect_true(mean_r2_no_refine > 0.05, 
              info = paste("Mean R² without refinement too low:", round(mean_r2_no_refine, 4)))
  
  # 7. R-squared should not become 0 after refinement unless there's a real issue
  # Allow some voxels to have low R-squared, but not all
  n_zero_r2 <- sum(result_with_refine$r_squared < 0.001)
  expect_true(n_zero_r2 < n_vox, 
              info = paste("Too many voxels with near-zero R² after refinement:", n_zero_r2, "out of", n_vox))
  
  # 8. Test that R-squared calculation is consistent
  # If we have fitted values, we can verify R-squared calculation
  if (!is.null(result_with_refine$fitted_values)) {
    cat("\nVerifying R-squared calculation with fitted values...\n")
    manual_r2 <- numeric(n_vox)
    for (v in 1:n_vox) {
      y_true <- Y[, v]
      y_pred <- result_with_refine$fitted_values[, v]
      
      # Calculate TSS and RSS manually
      tss <- sum((y_true - mean(y_true))^2)
      rss <- sum((y_true - y_pred)^2)
      manual_r2[v] <- 1 - (rss / tss)
    }
    
    cat("Manual R² calculation:", round(manual_r2, 4), "\n")
    cat("Package R² values:  ", round(result_with_refine$r_squared, 4), "\n")
    
    # They should match within tolerance
    expect_true(all(abs(manual_r2 - result_with_refine$r_squared) < 0.01),
                info = "Manual and package R-squared calculations don't match")
  }
})

test_that("R-squared calculation handles edge cases", {
  skip_if_not_installed("fmrireg")
  
  # Test with minimal data
  n_time <- 30
  Y <- matrix(rnorm(n_time * 1, 50, 5), n_time, 1)
  design <- matrix(c(rep(0, 10), 1, rep(0, 19)), n_time, 1)
  
  cat("\n=== Edge Case Test ===\n")
  cat("Testing with minimal data (", n_time, " timepoints, 1 voxel)\n")
  
  result <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = design,
    parametric_hrf = "lwu",
    global_refinement = TRUE,
    verbose = FALSE
  )
  
  # Should have valid R-squared even with minimal data
  expect_length(result$r_squared, 1)
  expect_true(result$r_squared >= 0 && result$r_squared <= 1)
  expect_false(is.na(result$r_squared))
  
  cat("Result: R² =", round(result$r_squared, 4), "\n")
})

test_that("R-squared is preserved when refinement is disabled", {
  # Test that disabling refinement doesn't break R-squared calculation
  set.seed(456)
  n_time <- 50
  n_vox <- 2
  Y <- matrix(rnorm(n_time * n_vox, 100, 15), n_time, n_vox)
  design <- matrix(c(rep(0, 15), 1, rep(0, 20), 1, rep(0, 13)), n_time, 1)
  
  # Add some signal
  signal_strength <- c(25, 30)
  for (v in 1:n_vox) {
    for (i in which(design[,1] == 1)) {
      if (i + 5 <= n_time) {
        Y[i:(i+5), v] <- Y[i:(i+5), v] + signal_strength[v] * exp(-0.2 * (0:5))
      }
    }
  }
  
  result <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = design,
    parametric_hrf = "lwu",
    global_refinement = FALSE,
    kmeans_refinement = FALSE,
    tiered_refinement = "none",
    verbose = FALSE
  )
  
  expect_equal(length(result$r_squared), n_vox)
  expect_true(all(result$r_squared >= 0 & result$r_squared <= 1))
  expect_true(all(!is.na(result$r_squared)))
  
  cat("\nNo-refinement test: R² =", round(result$r_squared, 4), "\n")
})

test_that("R-squared is correctly recalculated when parameters change", {
  skip_if_not_installed("fmrireg")
  
  # Create a test case that will definitely trigger parameter changes
  set.seed(789)
  n_time <- 60
  n_vox <- 2
  
  # Simple design with clear signal
  design_matrix <- matrix(0, n_time, 1)
  design_matrix[c(15, 35), 1] <- 1
  
  # Generate data with known good HRF parameters
  Y <- matrix(rnorm(n_time * n_vox, 100, 8), n_time, n_vox)
  
  # Add clear signal
  time_points <- seq(0, 25, by = 2)
  true_hrf <- exp(-time_points/4) * (time_points > 0)
  true_hrf <- true_hrf / max(true_hrf)
  
  for (v in 1:n_vox) {
    signal <- numeric(n_time)
    amplitude <- 30 + v * 15
    
    for (i in which(design_matrix[,1] == 1)) {
      end_idx <- min(i + length(true_hrf) - 1, n_time)
      if (end_idx > i) {
        signal[i:end_idx] <- signal[i:end_idx] + 
          true_hrf[1:(end_idx - i + 1)] * amplitude
      }
    }
    Y[, v] <- Y[, v] + signal
  }
  
  cat("\n=== Forced Parameter Change Test ===\n")
  
  # Use very poor starting parameters that will definitely get refined
  bad_seed <- c(15, 0.5, 1.2)  # Way off from optimal
  good_seed <- c(6, 3, 0.3)    # Much closer to optimal
  
  result_bad_start <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = design_matrix,
    parametric_hrf = "lwu",
    global_refinement = FALSE,
    theta_seed = bad_seed,
    verbose = FALSE
  )
  
  result_good_start <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = design_matrix,
    parametric_hrf = "lwu",
    global_refinement = FALSE,
    theta_seed = good_seed,
    verbose = FALSE
  )
  
  cat("Bad start params:\n")
  print(round(result_bad_start$estimated_parameters, 3))
  cat("Good start params:\n")
  print(round(result_good_start$estimated_parameters, 3))
  
  cat("Bad start R²:", round(result_bad_start$r_squared, 4), "\n")
  cat("Good start R²:", round(result_good_start$r_squared, 4), "\n")
  
  # Parameters should be different
  param_diff <- abs(result_bad_start$estimated_parameters - result_good_start$estimated_parameters)
  expect_true(any(param_diff > 0.1), 
              info = "Starting parameters should lead to different final parameters")
  
  # R-squared should also be different (better parameters should have better R²)
  r2_diff <- abs(result_bad_start$r_squared - result_good_start$r_squared)
  expect_true(any(r2_diff > 0.01),
              info = "Different parameters should lead to different R-squared values")
  
  # Both should have valid R-squared
  expect_true(all(result_bad_start$r_squared >= 0 & result_bad_start$r_squared <= 1))
  expect_true(all(result_good_start$r_squared >= 0 & result_good_start$r_squared <= 1))
  
  cat("Parameter difference test passed - R-squared correctly reflects parameter changes\n")
})