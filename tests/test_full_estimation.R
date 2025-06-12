# Test full estimation pipeline with the R² fix

library(fmriparametric)
library(fmrireg)

cat("=== Testing Full Estimation Pipeline ===\n\n")

# Create test data
set.seed(123)
n_time <- 100
n_vox <- 5

# Event design
events <- seq(10, 90, by = 20)
design_matrix <- matrix(0, n_time, 1)
design_matrix[events, 1] <- 1

# Create data with baseline and signal
Y <- matrix(rnorm(n_time * n_vox, mean = 1000, sd = 50), n_time, n_vox)

# Add HRF responses
true_hrf <- fmrireg::hrf_spmg1(seq(0, 20, by = 1))
for (event in events) {
  end_idx <- min(event + length(true_hrf) - 1, n_time)
  hrf_len <- end_idx - event + 1
  Y[event:end_idx, ] <- Y[event:end_idx, ] + 
    matrix(true_hrf[1:hrf_len] * 100, nrow = hrf_len, ncol = n_vox)
}

cat("Data properties:\n")
cat("Y mean:", mean(Y), "\n")
cat("Y range:", range(Y), "\n")
cat("Events at:", events, "\n\n")

# Run estimation
cat("Running parametric HRF estimation...\n")
result <- estimate_parametric_hrf(
  fmri_data = Y,
  event_model = design_matrix,
  parametric_hrf = "lwu",
  theta_seed = c(6, 1, 0.5),
  global_refinement = TRUE,
  tiered_refinement = "moderate",
  compute_se = TRUE,
  verbose = FALSE
)

cat("\nRESULTS:\n")
cat("R-squared values:", round(result$r_squared, 3), "\n")
cat("Mean R-squared:", round(mean(result$r_squared), 3), "\n")
cat("Min R-squared:", round(min(result$r_squared), 3), "\n")
cat("Max R-squared:", round(max(result$r_squared), 3), "\n\n")

cat("Parameter estimates:\n")
print(round(result$estimated_parameters, 2))

cat("\nAmplitudes:", round(result$amplitudes, 1), "\n")

# Check if refinement improved fits
if (!is.null(result$refinement_info)) {
  cat("\nRefinement info:\n")
  cat("Easy voxels:", result$refinement_info$n_easy, "\n")
  cat("Moderate voxels:", result$refinement_info$n_moderate, "\n")
  cat("Hard voxels:", result$refinement_info$n_hard, "\n")
}

# Test summary method
cat("\n")
summary(result)