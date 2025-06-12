# Debug Gauss-Newton directly
library(devtools)
load_all(quiet = TRUE)

# Create simple test case
set.seed(123)
n_vox <- 5
n_time <- 100
n_params <- 3

# Create test data
Y_proj <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
S_target_proj <- matrix(c(rep(0, 20), rep(1, 10), rep(0, 70)), ncol = 1)
theta_hat <- matrix(c(5, 2, 0.5), nrow = n_vox, ncol = n_params, byrow = TRUE)
colnames(theta_hat) <- c("tau", "sigma", "rho")

# Add some noise to parameters
theta_hat <- theta_hat + matrix(rnorm(n_vox * n_params, 0, 0.1), n_vox, n_params)

r2_voxel <- runif(n_vox, 0, 0.3)
hrf_interface <- .get_hrf_interface("lwu")
theta_bounds <- hrf_interface$default_bounds()
queue_labels <- rep("hard_GN", n_vox)

cat("Test data created:\n")
cat("Y_proj dimensions:", dim(Y_proj), "\n")
cat("S_target_proj dimensions:", dim(S_target_proj), "\n")
cat("theta_hat dimensions:", dim(theta_hat), "\n")
cat("theta_hat:\n")
print(round(theta_hat, 3))

# Test Gauss-Newton directly
cat("\nTesting Gauss-Newton refinement...\n")
result <- tryCatch({
  .gauss_newton_refinement(
    theta_hat_voxel = theta_hat,
    r2_voxel = r2_voxel,
    Y_proj = Y_proj,
    S_target_proj = S_target_proj,
    scan_times = seq_len(n_time),
    hrf_eval_times = 0:30,
    hrf_interface = hrf_interface,
    theta_bounds = theta_bounds,
    queue_labels = queue_labels,
    max_iter_gn = 2,
    verbose = TRUE
  )
}, error = function(e) {
  cat("Error in Gauss-Newton:", e$message, "\n")
  cat("\nTraceback:\n")
  print(traceback())
  NULL
})

if (!is.null(result)) {
  cat("\nGauss-Newton succeeded!\n")
  cat("Refined parameters:\n")
  print(round(result$theta_hat[1:5,], 3))
}