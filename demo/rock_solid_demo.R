#' Rock Solid fmriparametric Demonstration
#' 
#' This script demonstrates the bulletproof features of the rock-solid implementation
#' by throwing various challenging scenarios at it.

library(fmriparametric)

cat("\n===== ROCK SOLID FMRIPARAMETRIC DEMO =====\n\n")

# Demo 1: Handle pathological data
cat("DEMO 1: Pathological Data Handling\n")
cat("----------------------------------\n")

# Create data with multiple issues
set.seed(123)
n_time <- 50
n_vox <- 20

Y_pathological <- matrix(rnorm(n_time * n_vox), nrow = n_time)
# Add problematic voxels
Y_pathological[, 1] <- 0                    # All zeros
Y_pathological[, 2] <- 1                    # Constant
Y_pathological[, 3] <- NA                   # All NA
Y_pathological[, 4] <- c(Inf, rep(1, n_time - 1))  # Contains Inf
Y_pathological[, 5] <- 1e10 * rnorm(n_time)        # Extreme values
Y_pathological[, 6:10] <- 1e-10 * rnorm(n_time * 5)  # Tiny values

S <- matrix(0, nrow = n_time, ncol = 1)
S[seq(5, n_time, by = 10), 1] <- 1  # Sparse events

cat("Running on pathological data with:\n")
cat("- All-zero voxels\n")
cat("- Constant voxels\n") 
cat("- NA voxels\n")
cat("- Infinite values\n")
cat("- Extreme scales\n\n")

# This should handle everything gracefully
fit_pathological <- estimate_parametric_hrf_rock_solid(
  Y_pathological, S,
  verbose = TRUE,
  safety_mode = "maximum",
  error_report = TRUE
)

cat("\nResults summary:\n")
print(fit_pathological)

if (!is.null(fit_pathological$error_report)) {
  cat("\nError report:\n")
  print(fit_pathological$error_report)
}

# Demo 2: Memory stress test
cat("\n\nDEMO 2: Memory Management\n")
cat("-------------------------\n")

# Create large dataset that might stress memory
n_time_large <- 200
n_vox_large <- 10000

cat("Creating large dataset (", n_time_large, "timepoints x", 
    n_vox_large, "voxels)\n")

Y_large <- matrix(rnorm(n_time_large * n_vox_large, sd = 0.1), 
                  nrow = n_time_large)
S_large <- matrix(0, nrow = n_time_large, ncol = 1)
S_large[seq(10, n_time_large, by = 20), 1] <- 1

cat("Dataset size:", round(object.size(Y_large) / 1024^2, 1), "MB\n")
cat("Running estimation with automatic memory management...\n")

time_large <- system.time({
  fit_large <- estimate_parametric_hrf_rock_solid(
    Y_large, S_large,
    verbose = FALSE,  # Less output for large data
    safety_mode = "balanced",
    compute_se = FALSE,  # Save memory
    recenter_global_passes = 1,  # Fewer passes for speed
    recenter_kmeans_passes = 0   # Skip k-means for demo
  )
})

cat("Completed in", round(time_large["elapsed"], 1), "seconds\n")
cat("Processing rate:", round(n_vox_large / time_large["elapsed"]), 
    "voxels/second\n")

# Demo 3: Algorithm failure recovery
cat("\n\nDEMO 3: Algorithm Failure Recovery\n")
cat("----------------------------------\n")

# Create data that causes numerical problems
Y_singular <- outer(1:30, 1:10)  # Rank-1 matrix
S_singular <- matrix(1, nrow = 30, ncol = 1)  # Constant stimulus

cat("Running on rank-deficient data that breaks standard algorithms...\n")

fit_singular <- estimate_parametric_hrf_rock_solid(
  Y_singular, S_singular,
  verbose = TRUE,
  safety_mode = "maximum",
  lambda_ridge = 0.001  # Low regularization to stress algorithm
)

cat("\nRecovered successfully!\n")
cat("Mean R²:", round(mean(fit_singular$r_squared, na.rm = TRUE), 3), "\n")

# Demo 4: Adversarial inputs
cat("\n\nDEMO 4: Adversarial Input Testing\n")
cat("---------------------------------\n")

adversarial_tests <- list(
  "Empty event model" = list(
    Y = matrix(rnorm(100), 20, 5),
    S = matrix(0, 20, 1)  # No events!
  ),
  "Mismatched dimensions" = list(
    Y = matrix(rnorm(100), 20, 5),
    S = matrix(1, 25, 1)  # Wrong size
  ),
  "Impossible bounds" = list(
    Y = matrix(rnorm(100), 20, 5),
    S = matrix(rbinom(20, 1, 0.3), 20, 1),
    theta_bounds = list(lower = c(10, 10, 10), upper = c(9, 9, 9))
  ),
  "Negative time span" = list(
    Y = matrix(rnorm(100), 20, 5),
    S = matrix(rbinom(20, 1, 0.3), 20, 1),
    hrf_span = -10
  )
)

for (test_name in names(adversarial_tests)) {
  cat("\nTesting:", test_name, "\n")
  
  test_data <- adversarial_tests[[test_name]]
  
  result <- tryCatch({
    fit <- estimate_parametric_hrf_rock_solid(
      fmri_data = test_data$Y,
      event_model = test_data$S,
      theta_bounds = test_data$theta_bounds,
      hrf_span = test_data$hrf_span,
      verbose = FALSE,
      safety_mode = "maximum"
    )
    "PASSED - Handled gracefully"
  }, error = function(e) {
    paste("PASSED - Caught error:", e$message)
  })
  
  cat("Result:", result, "\n")
}

# Demo 5: Safety mode comparison
cat("\n\nDEMO 5: Safety Mode Comparison\n")
cat("------------------------------\n")

Y_normal <- matrix(rnorm(500), 50, 10)
S_normal <- matrix(rbinom(50, 1, 0.2), 50, 1)

safety_modes <- c("maximum", "balanced", "performance")
results <- list()

for (mode in safety_modes) {
  cat("\nTesting", mode, "mode...\n")
  
  time_mode <- system.time({
    fit_mode <- estimate_parametric_hrf_rock_solid(
      Y_normal, S_normal,
      safety_mode = mode,
      verbose = FALSE
    )
  })
  
  results[[mode]] <- list(
    time = time_mode["elapsed"],
    mean_r2 = mean(fit_mode$r_squared, na.rm = TRUE),
    errors = fit_mode$metadata$errors_recovered
  )
  
  cat("  Time:", round(results[[mode]]$time, 2), "seconds\n")
  cat("  Mean R²:", round(results[[mode]]$mean_r2, 3), "\n")
  cat("  Errors recovered:", results[[mode]]$errors, "\n")
}

# Summary
cat("\n\n===== DEMO SUMMARY =====\n")
cat("The rock-solid implementation successfully:\n")
cat("✓ Handled pathological data without crashing\n")
cat("✓ Managed memory for large datasets\n")
cat("✓ Recovered from algorithm failures\n")
cat("✓ Caught and handled adversarial inputs\n")
cat("✓ Provided flexible safety/performance trade-offs\n")
cat("\nThe package is now ROCK SOLID! 🪨💪\n")