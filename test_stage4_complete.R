library(devtools)
load_all(quiet = TRUE)
set.seed(123)

# Create realistic test data
n_time <- 200
n_vox <- 20
TR <- 2

# True HRF parameters
true_tau <- runif(n_vox, 4, 8)
true_sigma <- runif(n_vox, 1, 3)
true_rho <- runif(n_vox, 0.2, 0.8)
true_amp <- runif(n_vox, 0.5, 2)

# Create events
events <- matrix(0, n_time, 1)
events[seq(20, 180, by=40), 1] <- 1

# Generate data
Y <- matrix(0, n_time, n_vox)
hrf_interface <- .get_hrf_interface('lwu')
for (v in 1:n_vox) {
  hrf <- hrf_interface$hrf_function(0:30, c(true_tau[v], true_sigma[v], true_rho[v]))
  conv_signal <- stats::filter(events[,1], hrf, sides=1)
  conv_signal[is.na(conv_signal)] <- 0
  Y[,v] <- true_amp[v] * conv_signal + rnorm(n_time, 0, 0.2)
}

# Test with tiered refinement
cat("Testing Stage 4 with tiered refinement...\n")
fit <- estimate_parametric_hrf(
  Y, events,
  tiered_refinement = 'moderate',
  global_refinement = TRUE,
  global_passes = 2,
  verbose = TRUE
)

cat('\n=== Results ===\n')
cat('Mean R²:', round(mean(fit$r_squared), 3), '\n')
cat('R² range:', round(range(fit$r_squared), 3), '\n')
if (!is.null(fit$refinement_info)) {
  cat('\nVoxel classification:\n')
  cat('  Easy:', fit$refinement_info$n_easy, '\n')
  cat('  Moderate:', fit$refinement_info$n_moderate, '\n')
  cat('  Hard:', fit$refinement_info$n_hard, '\n')
  if (fit$refinement_info$n_moderate > 0) {
    cat('  Moderate refined:', fit$refinement_info$moderate_refined, '\n')
  }
}

# Compare parameter recovery
param_recovery <- cbind(
  tau_error = mean(abs(fit$estimated_parameters[,'tau'] - true_tau)),
  sigma_error = mean(abs(fit$estimated_parameters[,'sigma'] - true_sigma)),
  rho_error = mean(abs(fit$estimated_parameters[,'rho'] - true_rho))
)
cat('\nParameter recovery (mean absolute error):\n')
print(round(param_recovery, 3))

# Check if refinement improved fit
if (!is.null(fit$convergence_info$global_iterations)) {
  cat('\nGlobal refinement iterations:', fit$convergence_info$global_iterations, '\n')
}

cat('\nStage 4 implementation complete!\n')