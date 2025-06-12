# Debug script for R-squared = 0 issue
# Following Gemini's systematic approach

library(fmriparametric)
library(fmrireg)

cat("=== Debugging R-squared = 0 issue ===\n\n")

# Set debug mode
options(fmriparametric.debug = TRUE)

# 1. Create simple test data with known signal
set.seed(42)
n_time <- 100
n_vox <- 3

# Create event design with clear structure
events <- c(10, 30, 50, 70, 90)
design_matrix <- matrix(0, n_time, 1)
design_matrix[events, 1] <- 1

cat("Event design:\n")
cat("Events at time points:", events, "\n")
cat("Design matrix sum:", sum(design_matrix), "\n\n")

# Create synthetic BOLD data with strong signal
Y <- matrix(rnorm(n_time * n_vox, mean = 100, sd = 10), n_time, n_vox)

# Add a clear HRF-like response to events
true_hrf <- fmrireg::hrf_spmg1(seq(0, 20, by = 1))
for (event_time in events) {
  end_time <- min(event_time + length(true_hrf) - 1, n_time)
  hrf_length <- end_time - event_time + 1
  Y[event_time:end_time, ] <- Y[event_time:end_time, ] + 
    matrix(true_hrf[1:hrf_length] * 50, nrow = hrf_length, ncol = n_vox)
}

cat("Data properties:\n")
cat("Y dimensions:", dim(Y), "\n")
cat("Y range:", range(Y), "\n")
cat("Y variance by voxel:", apply(Y, 2, var), "\n\n")

# 2. Run estimation with minimal configuration
cat("Running parametric HRF estimation...\n")
cat("================================\n")

result <- tryCatch({
  estimate_parametric_hrf(
    fmri_data = Y,
    event_model = design_matrix,
    parametric_hrf = "lwu",
    theta_seed = c(6, 1, 0.5),  # tau, sigma, rho
    global_refinement = FALSE,   # Disable refinement to focus on core engine
    tiered_refinement = "none",
    compute_se = FALSE,
    verbose = TRUE
  )
}, error = function(e) {
  cat("\nERROR in estimation:", e$message, "\n")
  NULL
})

if (!is.null(result)) {
  cat("\n\n=== RESULTS ===\n")
  cat("R-squared values:", result$r_squared, "\n")
  cat("Mean R-squared:", mean(result$r_squared), "\n")
  cat("Parameter estimates:\n")
  print(result$estimated_parameters)
  cat("Amplitudes:", result$amplitudes, "\n")
} else {
  cat("\nEstimation failed!\n")
}

# 3. Test the prepare_parametric_inputs directly
cat("\n\n=== Testing prepare_parametric_inputs ===\n")
inputs <- fmriparametric:::.prepare_parametric_inputs(
  fmri_data = Y,
  event_model = design_matrix,
  hrf_span = 30
)

cat("Prepared inputs:\n")
cat("Y_proj dim:", dim(inputs$Y_proj), "\n")
cat("S_target dim:", dim(inputs$S_target), "\n")
cat("S_target_proj dim:", dim(inputs$S_target_proj), "\n")
cat("S_target sum:", sum(abs(inputs$S_target)), "\n")
cat("S_target_proj sum:", sum(abs(inputs$S_target_proj)), "\n")

# Check if S_target_proj is all zeros
if (all(abs(inputs$S_target_proj) < 1e-10)) {
  cat("\nWARNING: S_target_proj is all zeros! This would cause R² = 0\n")
  cat("Checking S_target (before projection):\n")
  print(head(inputs$S_target, 20))
}

# Reset debug mode
options(fmriparametric.debug = FALSE)