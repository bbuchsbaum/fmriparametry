# Simplified debug test for R-squared = 0 issue
# Testing core components directly

library(fmriparametric)
library(fmrihrf)
library(Matrix)

cat("=== Direct Component Testing for R² = 0 Issue ===\n\n")

# 1. Create simple test data
set.seed(42)
n_time <- 50
n_vox <- 2

# Simple event design 
events <- c(10, 25, 40)
S_target <- matrix(0, n_time, 1)
S_target[events, 1] <- 1

cat("Event design - events at:", events, "\n")
cat("S_target sum:", sum(S_target), "\n\n")

# Create data with clear signal
Y <- matrix(rnorm(n_time * n_vox, 100, 5), n_time, n_vox)

# Add HRF response
hrf_times <- seq(0, 20, by = 1)
true_hrf <- fmrihrf::hrf_spmg1(hrf_times)
for (event in events) {
  end_idx <- min(event + length(true_hrf) - 1, n_time)
  hrf_len <- end_idx - event + 1
  Y[event:end_idx, ] <- Y[event:end_idx, ] + 
    matrix(true_hrf[1:hrf_len] * 30, nrow = hrf_len, ncol = n_vox)
}

cat("Y data properties:\n")
cat("  Range:", range(Y), "\n")
cat("  Mean:", mean(Y), "\n")
cat("  Variance by voxel:", apply(Y, 2, var), "\n\n")

# 2. Test HRF interface creation
cat("Testing HRF interface...\n")
hrf_interface <- fmriparametric:::.get_hrf_interface("lwu")
cat("  Parameter names:", hrf_interface$parameter_names, "\n")
cat("  Default seed:", hrf_interface$default_seed(), "\n\n")

# 3. Test Taylor basis generation
cat("Testing Taylor basis generation...\n")
theta_seed <- c(6, 1, 0.5)
hrf_eval_times <- seq(0, 30, by = 2)
taylor_basis <- hrf_interface$taylor_basis(theta_seed, hrf_eval_times)

cat("  Taylor basis dim:", dim(taylor_basis), "\n")
cat("  HRF values (col 1) range:", range(taylor_basis[,1]), "\n")
cat("  HRF sum:", sum(taylor_basis[,1]), "\n\n")

# 4. Test batch convolution
cat("Testing batch convolution...\n")
X_design <- fmriparametric:::.batch_convolution(S_target, taylor_basis, n_time)
cat("  X_design dim:", dim(X_design), "\n")
cat("  X_design[,1] range:", range(X_design[,1]), "\n")
cat("  X_design[,1] variance:", var(X_design[,1]), "\n\n")

# 5. Test parametric engine directly
cat("Testing parametric engine...\n")
options(fmriparametric.debug = TRUE)

engine_result <- fmriparametric:::.parametric_engine(
  Y_proj = Y,
  S_target_proj = S_target,
  hrf_eval_times = hrf_eval_times,
  hrf_interface = hrf_interface,
  theta_seed = theta_seed,
  theta_bounds = list(lower = c(0, 0.1, 0), upper = c(20, 10, 1.5)),
  lambda_ridge = 0.01,
  verbose = TRUE
)

cat("\nEngine results:\n")
cat("  R-squared:", engine_result$r_squared, "\n")
cat("  Beta0 (amplitudes):", engine_result$beta0, "\n")
cat("  Theta estimates:\n")
print(engine_result$theta_hat)

# 6. Manually compute R-squared to verify
cat("\n\nManual R-squared calculation:\n")
fitted <- X_design %*% engine_result$coeffs
residuals <- Y - fitted
ss_res <- colSums(residuals^2)
ss_tot <- colSums(scale(Y, scale = FALSE)^2)
r2_manual <- 1 - ss_res / ss_tot

cat("  Manual R²:", r2_manual, "\n")
cat("  SS_res:", ss_res, "\n")
cat("  SS_tot:", ss_tot, "\n")

# Check if the issue is with the data projection
if (all(engine_result$r_squared < 0.01)) {
  cat("\n\nDEBUG: Low R² detected. Checking data flow:\n")
  cat("1. Original Y variance:", apply(Y, 2, var), "\n")
  cat("2. S_target non-zero elements:", sum(S_target != 0), "\n")
  cat("3. X_design non-zero proportion:", mean(X_design != 0), "\n")
  cat("4. Coefficients range:", range(engine_result$coeffs), "\n")
  cat("5. Fitted values range:", range(fitted), "\n")
}

options(fmriparametric.debug = FALSE)