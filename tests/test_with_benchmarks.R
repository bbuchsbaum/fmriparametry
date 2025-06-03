# Test fmriparametric with fmrireg benchmark datasets
library(fmrireg)
library(testthat)

# Source our functions
source("R/parametric-engine.R")
source("R/hrf-interface-lwu.R")
source("R/prepare-parametric-inputs.R")
source("R/parametric-hrf-fit-class.R")
source("R/estimate_parametric_hrf.R")

# List available benchmarks
cat("Available benchmark datasets:\n")
print(fmrireg:::list_benchmark_datasets())

# Test 1: High SNR Canonical HRF
cat("\n\n=== TEST 1: High SNR Canonical HRF ===\n")
bm_high <- fmrireg::get_benchmark_dataset("BM_Canonical_HighSNR")
summary_high <- fmrireg:::get_benchmark_summary("BM_Canonical_HighSNR")

cat("Dataset dimensions:", summary_high$dimensions$n_timepoints, "timepoints x", 
    summary_high$dimensions$n_voxels, "voxels\n")
cat("SNR:", summary_high$experimental_design$target_snr, "\n")

# The benchmark provides an SPMG1 HRF, but we need to test our LWU model
# Let's create a simple test with the benchmark structure

# Extract data
Y <- bm_high$data  # Should be timepoints x voxels
event_model <- bm_high$event_model

# Try basic estimation
cat("\nTesting basic parametric engine...\n")
tryCatch({
  # First test just the engine
  hrf_interface <- list(
    hrf_function = .lwu_hrf_function,
    taylor_basis = .lwu_hrf_taylor_basis_function,
    parameter_names = .lwu_hrf_parameter_names(),
    default_seed = .lwu_hrf_default_seed(),
    default_bounds = .lwu_hrf_default_bounds()
  )
  
  # Get design matrix from event model
  S <- fmrireg::construct_design_matrix(
    event_model, 
    hrf = fmrireg::HRF_IDENTITY,  # Use identity to get raw design
    sample_times = seq(0, (nrow(Y)-1)*2, by = 2)
  )
  
  # Test on a few voxels first
  n_test_voxels <- 10
  fit_engine <- .parametric_engine(
    Y_proj = Y[, 1:n_test_voxels],
    S_target_proj = S[, 1, drop = FALSE],  # Use first regressor
    hrf_eval_times = seq(0, 30, length.out = 61),
    hrf_interface = hrf_interface,
    theta_seed = c(6, 2.5, 0.35),
    theta_bounds = list(lower = c(2, 1, 0), upper = c(12, 5, 1)),
    lambda_ridge = 0.01,
    verbose = TRUE
  )
  
  cat("Success! Engine returned estimates for", nrow(fit_engine$theta_hat), "voxels\n")
  cat("Parameter estimates (first voxel):\n")
  cat("  tau (peak time):", round(fit_engine$theta_hat[1, 1], 2), "s\n")
  cat("  sigma (width):", round(fit_engine$theta_hat[1, 2], 2), "s\n")
  cat("  rho (undershoot):", round(fit_engine$theta_hat[1, 3], 2), "\n")
  cat("  R²:", round(fit_engine$r_squared[1], 3), "\n")
  
}, error = function(e) {
  cat("ERROR in parametric engine:", e$message, "\n")
})

# Test 2: Low SNR to check robustness
cat("\n\n=== TEST 2: Low SNR Canonical HRF ===\n")
bm_low <- fmrireg::get_benchmark_dataset("BM_Canonical_LowSNR")
summary_low <- fmrireg:::get_benchmark_summary("BM_Canonical_LowSNR")
cat("SNR:", summary_low$experimental_design$target_snr, "\n")

# Test 3: HRF Variability - Perfect for K-means!
cat("\n\n=== TEST 3: HRF Variability Across Voxels ===\n")
bm_var <- fmrireg::get_benchmark_dataset("BM_HRF_Variability_AcrossVoxels")
summary_var <- fmrireg:::get_benchmark_summary("BM_HRF_Variability_AcrossVoxels")
cat("This dataset has HRF varying across voxel groups - perfect for testing K-means!\n")

# Test the full estimation pipeline
cat("\nTesting full estimation pipeline with estimate_parametric_hrf...\n")
tryCatch({
  fit_full <- estimate_parametric_hrf(
    fmri_data = bm_var$data,
    event_model = bm_var$event_model,
    parametric_hrf = "lwu",
    verbose = TRUE
  )
  
  cat("\nFull pipeline results:\n")
  print(fit_full)
  
  # Check if we can detect the groups
  cat("\nParameter distributions (should show clusters):\n")
  cat("Tau (peak time) range:", range(coef(fit_full)[, "tau"]), "\n")
  cat("Sigma (width) range:", range(coef(fit_full)[, "sigma"]), "\n")
  
}, error = function(e) {
  cat("ERROR in full pipeline:", e$message, "\n")
  cat("Traceback:\n")
  traceback()
})

# Test 4: Complex realistic scenario
cat("\n\n=== TEST 4: Complex Realistic Scenario ===\n")
bm_complex <- fmrireg::get_benchmark_dataset("BM_Complex_Realistic")
summary_complex <- fmrireg:::get_benchmark_summary("BM_Complex_Realistic")
cat("This has 3 HRF groups, 3 conditions, variable durations, and AR(2) noise\n")

# Generate our own test data with known LWU parameters
cat("\n\n=== TEST 5: Synthetic LWU Recovery Test ===\n")
create_lwu_test_data <- function(n_time = 100, n_vox = 50, 
                                 true_params = matrix(c(6, 2.5, 0.35), nrow = 1)) {
  # Create event design
  S <- matrix(0, n_time, 1)
  S[seq(10, n_time, by = 20), 1] <- 1
  
  # Generate HRF
  t_hrf <- seq(0, 30, length.out = 61)
  Y <- matrix(NA, n_time, n_vox)
  
  for (v in 1:n_vox) {
    # Add some variation
    theta_v <- true_params[1, ] + rnorm(3, sd = c(0.5, 0.2, 0.05))
    theta_v <- pmax(c(2, 1, 0), pmin(theta_v, c(12, 5, 1)))  # Bounds
    
    # Generate HRF using our function
    hrf <- .lwu_hrf_function(t_hrf, theta_v)
    
    # Convolve
    conv_full <- stats::convolve(S[, 1], rev(hrf), type = "open")
    signal <- conv_full[1:n_time]
    
    # Add noise
    Y[, v] <- signal + rnorm(n_time, sd = 0.2)
  }
  
  list(Y = Y, S = S, true_params = true_params)
}

test_data <- create_lwu_test_data()

cat("Testing parameter recovery on synthetic LWU data...\n")
tryCatch({
  fit_recovery <- .parametric_engine(
    Y_proj = test_data$Y,
    S_target_proj = test_data$S,
    hrf_eval_times = seq(0, 30, length.out = 61),
    hrf_interface = list(
      hrf_function = .lwu_hrf_function,
      taylor_basis = .lwu_hrf_taylor_basis_function,
      parameter_names = .lwu_hrf_parameter_names(),
      default_seed = .lwu_hrf_default_seed(),
      default_bounds = .lwu_hrf_default_bounds()
    ),
    theta_seed = c(6, 2.5, 0.35),
    theta_bounds = list(lower = c(2, 1, 0), upper = c(12, 5, 1)),
    verbose = TRUE
  )
  
  cat("\nParameter recovery results:\n")
  cat("True parameters:", test_data$true_params[1, ], "\n")
  cat("Mean estimated:", colMeans(fit_recovery$theta_hat), "\n")
  cat("Mean absolute error:", mean(abs(fit_recovery$theta_hat - 
                                       matrix(test_data$true_params[1, ], 
                                              nrow = nrow(fit_recovery$theta_hat), 
                                              ncol = 3, byrow = TRUE))), "\n")
  cat("Mean R²:", mean(fit_recovery$r_squared), "\n")
  
}, error = function(e) {
  cat("ERROR in synthetic recovery:", e$message, "\n")
})

# Performance test
cat("\n\n=== PERFORMANCE TEST ===\n")
cat("Testing speed on different data sizes...\n")

for (n_vox in c(100, 1000)) {
  test_data <- create_lwu_test_data(n_time = 100, n_vox = n_vox)
  
  time_taken <- system.time({
    fit <- .parametric_engine(
      Y_proj = test_data$Y,
      S_target_proj = test_data$S,
      hrf_eval_times = seq(0, 30, length.out = 61),
      hrf_interface = list(
        hrf_function = .lwu_hrf_function,
        taylor_basis = .lwu_hrf_taylor_basis_function,
        parameter_names = .lwu_hrf_parameter_names(),
        default_seed = .lwu_hrf_default_seed(),
        default_bounds = .lwu_hrf_default_bounds()
      ),
      theta_seed = c(6, 2.5, 0.35),
      theta_bounds = list(lower = c(2, 1, 0), upper = c(12, 5, 1)),
      verbose = FALSE
    )
  })
  
  cat(n_vox, "voxels:", round(time_taken["elapsed"], 2), "seconds (",
      round(n_vox / time_taken["elapsed"]), "voxels/sec)\n")
}