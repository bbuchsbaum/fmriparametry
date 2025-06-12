# Test the refactored estimation function

library(fmriparametric)

cat("=== Testing Refactored Estimation Function ===\n\n")

# Create simple test data
set.seed(42)
n_time <- 50
n_vox <- 3

# Events
events <- c(10, 25, 40)
design_matrix <- matrix(0, n_time, 1)
design_matrix[events, 1] <- 1

# Data with baseline
Y <- matrix(rnorm(n_time * n_vox, mean = 100, sd = 5), n_time, n_vox)

# Add signal
for (v in 1:n_vox) {
  for (e in events) {
    Y[e:(e+5), v] <- Y[e:(e+5), v] + c(0, 5, 10, 15, 10, 5)
  }
}

cat("Running legacy implementation...\n")
result_legacy <- estimate_parametric_hrf(
  fmri_data = Y,
  event_model = design_matrix,
  parametric_hrf = "lwu",
  global_refinement = FALSE,
  tiered_refinement = "none",
  compute_se = FALSE,
  verbose = FALSE,
  .implementation = "legacy"
)

cat("Running refactored implementation...\n")
result_refactored <- estimate_parametric_hrf(
  fmri_data = Y,
  event_model = design_matrix,
  parametric_hrf = "lwu",
  global_refinement = FALSE,
  tiered_refinement = "none",
  compute_se = FALSE,
  verbose = FALSE,
  .implementation = "refactored"
)

cat("\n=== COMPARISON ===\n")
cat("Legacy R²:", round(result_legacy$r_squared, 3), "\n")
cat("Refactored R²:", round(result_refactored$r_squared, 3), "\n")
cat("Difference:", max(abs(result_legacy$r_squared - result_refactored$r_squared)), "\n")

cat("\nLegacy parameters:\n")
print(round(result_legacy$estimated_parameters, 3))
cat("\nRefactored parameters:\n")
print(round(result_refactored$estimated_parameters, 3))

cat("\nParameter difference:\n")
print(round(result_legacy$estimated_parameters - result_refactored$estimated_parameters, 6))

# Test compare mode
cat("\n\n=== Testing Compare Mode ===\n")
result_compare <- estimate_parametric_hrf(
  fmri_data = Y,
  event_model = design_matrix,
  parametric_hrf = "lwu",
  global_refinement = FALSE,
  tiered_refinement = "none",
  compute_se = FALSE,
  verbose = TRUE,
  .implementation = "compare"
)