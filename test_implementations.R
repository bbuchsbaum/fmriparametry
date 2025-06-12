# Test script to verify implementations work
library(devtools)
load_all()

# Create simple test data
set.seed(42)
n_time <- 50
n_vox <- 3

# Events
events <- c(10, 25, 40)
event_model <- matrix(0, n_time, 1)
event_model[events, 1] <- 1

# Data with baseline
Y <- matrix(rnorm(n_time * n_vox, mean = 100, sd = 5), n_time, n_vox)

# Add signal
for (v in 1:n_vox) {
  for (e in events) {
    if (e + 5 <= n_time) {
      Y[e:(e+5), v] <- Y[e:(e+5), v] + c(0, 5, 10, 15, 10, 5)
    }
  }
}

cat("Testing legacy implementation...\n")
result_legacy <- estimate_parametric_hrf(
  fmri_data = Y,
  event_model = event_model,
  parametric_hrf = "lwu",
  global_refinement = FALSE,
  tiered_refinement = "none",
  compute_se = FALSE,
  verbose = FALSE,
  .implementation = "legacy"
)
cat("Legacy R²:", round(result_legacy$r_squared, 3), "\n")

cat("\nTesting refactored implementation...\n")
result_refactored <- estimate_parametric_hrf(
  fmri_data = Y,
  event_model = event_model,
  parametric_hrf = "lwu",
  global_refinement = FALSE,
  tiered_refinement = "none",
  compute_se = FALSE,
  verbose = FALSE,
  .implementation = "refactored"
)
cat("Refactored R²:", round(result_refactored$r_squared, 3), "\n")

cat("\nComparing results...\n")
cat("R² match:", all.equal(result_legacy$r_squared, result_refactored$r_squared), "\n")
cat("Parameters match:", all.equal(result_legacy$estimated_parameters, result_refactored$estimated_parameters), "\n")