# Simple test without refinement to verify R² fix

library(fmriparametric)

cat("=== Simple Test (No Refinement) ===\n\n")

# Create simple data
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

cat("Y mean:", mean(Y), "\n")
cat("Y range:", range(Y), "\n\n")

# Run without refinement
result <- estimate_parametric_hrf(
  fmri_data = Y,
  event_model = design_matrix,
  parametric_hrf = "lwu",
  global_refinement = FALSE,
  tiered_refinement = "none",
  compute_se = FALSE,
  verbose = TRUE
)

cat("\n\nRESULTS:\n")
cat("R-squared:", round(result$r_squared, 3), "\n")
cat("Mean R-squared:", round(mean(result$r_squared), 3), "\n")

# Check fit object
cat("\nFit object structure:\n")
str(result, max.level = 1)