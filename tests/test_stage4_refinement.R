# Test script for Stage 4 tiered refinement
library(fmriparametric)
library(testthat)

# Create simple test data
set.seed(123)
n_time <- 100
n_vox <- 50
TR <- 1

# Create synthetic BOLD data with some signal
true_tau <- 5
true_sigma = 2
true_rho = 0.3
true_amp = 2

# Generate HRF
t_hrf <- seq(0, 30, by = TR)
true_hrf <- exp(-(t_hrf - true_tau)^2 / (2 * true_sigma^2)) - 
            true_rho * exp(-(t_hrf - true_tau - 2*true_sigma)^2 / (2 * (1.6*true_sigma)^2))
true_hrf <- true_hrf / max(true_hrf)

# Create stimulus vector
stim_times <- seq(10, 90, by = 20)
stim_vec <- numeric(n_time)
stim_vec[stim_times] <- 1

# Convolve to get signal
signal <- stats::filter(stim_vec, true_hrf, sides = 1)
signal[is.na(signal)] <- 0

# Create data matrix with varying SNR
Y <- matrix(0, n_time, n_vox)
for (v in 1:n_vox) {
  noise_level <- runif(1, 0.2, 1.5)  # Varying noise levels
  Y[, v] <- true_amp * signal + rnorm(n_time, sd = noise_level)
}

# Create event model with HRF term
event_df <- data.frame(
  onset = (stim_times - 1) * TR,
  duration = rep(1, length(stim_times)),
  trial_type = "stimulus",
  block = 1
)

# Use fmrireg if available, otherwise create mock
if (requireNamespace("fmrireg", quietly = TRUE)) {
  event_model <- fmrireg::event_model(
    onset ~ fmrireg::hrf(trial_type),
    data = event_df,
    block = event_df$block,
    sampling_frame = fmrireg::sampling_frame(n_time, TR = TR)
  )
  
  fmri_data <- fmrireg::matrix_dataset(Y, TR = TR, run_length = n_time)
} else {
  # Mock objects
  event_model <- list(
    terms = list(matrix(stim_vec, ncol = 1)),
    sampling_rate = 1/TR
  )
  class(event_model) <- c("event_model", "list")
  
  fmri_data <- list(
    data = Y,
    TR = TR,
    run_length = n_time
  )
  class(fmri_data) <- c("matrix_dataset", "list")
}

# Test with tiered refinement
cat("\n=== Testing estimate_parametric_hrf with tiered refinement ===\n")
result <- estimate_parametric_hrf(
  fmri_data = fmri_data,
  event_model = event_model,
  parametric_hrf = "lwu",
  tiered_refinement = "aggressive",
  refinement_thresholds = list(
    r2_easy = 0.7,
    r2_hard = 0.3,
    se_low = 0.3,
    se_high = 0.7,
    gauss_newton_maxiter = 5
  ),
  global_refinement = TRUE,
  global_passes = 2,
  verbose = TRUE
)

# Check results
cat("\n=== Results Summary ===\n")
cat("Mean R²:", mean(result$r_squared), "\n")
cat("Refinement info:\n")
if (!is.null(result$refinement_info)) {
  print(result$refinement_info)
}

# Test that refinement actually ran
test_that("Tiered refinement was executed", {
  expect_true(!is.null(result$refinement_info))
  expect_true(result$refinement_info$n_easy >= 0)
  expect_true(result$refinement_info$n_moderate >= 0)
  expect_true(result$refinement_info$n_hard >= 0)
})

cat("\n=== Test completed successfully ===\n")