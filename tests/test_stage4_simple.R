# Simplified test for Stage 4 refinement
library(fmriparametric)

cat("\n=== Testing Stage 4 refinement with mock data ===\n")

# Create very simple mock data
set.seed(42)
n_time <- 100
n_vox <- 20

# Create mock prepared data
Y_proj <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
S_target_proj <- matrix(c(rep(0, 10), rep(1, 5), rep(0, 85)), n_time, 1)

prepared_data <- list(
  inputs = list(
    Y_proj = Y_proj,
    S_target_proj = S_target_proj,
    hrf_eval_times = seq(0, 30, by = 1)
  ),
  n_vox = n_vox,
  n_time = n_time
)

# Create HRF interface
hrf_interface <- fmriparametric:::.get_hrf_interface("lwu")

# Create initial results from "Stage 3"
theta_current <- matrix(c(5, 2, 0.3), n_vox, 3, byrow = TRUE)
colnames(theta_current) <- hrf_interface$parameter_names

# Add some variation to make classification interesting
theta_current[1:5, 1] <- 4  # Some different lag values
theta_current[6:10, 2] <- 1.5  # Some different width values

amplitudes <- runif(n_vox, 0.5, 2)
r_squared <- c(
  runif(5, 0.8, 0.95),   # 5 easy voxels
  runif(10, 0.4, 0.65),  # 10 moderate voxels  
  runif(5, 0.1, 0.25)    # 5 hard voxels
)

# Mock SE values
se_theta <- matrix(runif(n_vox * 3, 0.1, 0.8), n_vox, 3)

# Create refined results from Stage 3
refined_results <- list(
  theta_current = theta_current,
  amplitudes = amplitudes,
  r_squared = r_squared,
  convergence_info = list(),
  se_theta = se_theta
)

# Create config
config <- list(
  tiered_refinement = "aggressive",
  refinement_thresholds = list(
    r2_easy = 0.7,
    r2_hard = 0.3,
    se_low = 0.3,
    se_high = 0.7,
    gauss_newton_maxiter = 3
  ),
  lambda_ridge = 0.01,
  theta_bounds = hrf_interface$default_bounds(),
  baseline_model = "intercept",
  compute_se = TRUE,
  verbose = TRUE
)

# Test Stage 4 directly
cat("\nCalling Stage 4...\n")
result <- fmriparametric:::.stage4_tiered_refinement(
  prepared_data = prepared_data,
  refined_results = refined_results,
  hrf_interface = hrf_interface,
  config = config
)

cat("\nStage 4 Results:\n")
cat("- Moderate refined:", result$refinement_info$moderate_refined, "\n")
cat("- Moderate improved:", result$refinement_info$moderate_improved, "\n")
cat("- Hard refined:", result$refinement_info$hard_refined, "\n")
cat("- Converged:", result$refinement_info$n_converged, "\n")
cat("- Improved:", result$refinement_info$n_improved, "\n")
cat("- Mean R² before:", mean(refined_results$r_squared), "\n")
cat("- Mean R² after:", mean(result$r_squared), "\n")

cat("\n=== Test completed ===\n")