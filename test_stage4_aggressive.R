library(devtools)
load_all(quiet = TRUE)
set.seed(123)

# Create realistic test data
n_time <- 200
n_vox <- 30  # More voxels to ensure we get some hard ones
TR <- 2

# True HRF parameters with more variation
true_tau <- runif(n_vox, 3, 10)     # Wider range
true_sigma <- runif(n_vox, 0.5, 4)  # Wider range
true_rho <- runif(n_vox, 0.1, 1.0)   # Wider range
true_amp <- runif(n_vox, 0.2, 3)     # Wider range

# Add some difficult voxels with low SNR
noise_levels <- c(rep(0.2, 20), rep(0.5, 5), rep(1.0, 5))

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
  Y[,v] <- true_amp[v] * conv_signal + rnorm(n_time, 0, noise_levels[v])
}

# Test with AGGRESSIVE tiered refinement (includes hard voxels)
cat("Testing Stage 4 with aggressive refinement (includes Gauss-Newton for hard voxels)...\n")
fit <- estimate_parametric_hrf(
  Y, events,
  tiered_refinement = 'aggressive',  # This will refine hard voxels too
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
  if (!is.null(fit$refinement_info$hard_refined)) {
    cat('  Hard refined:', fit$refinement_info$hard_refined, '\n')
  }
}

# Check convergence info for Gauss-Newton
if (!is.null(fit$convergence_info$gauss_newton)) {
  cat('\nGauss-Newton refinement:\n')
  cat('  Voxels processed:', fit$convergence_info$gauss_newton$n_hard, '\n')
  cat('  Converged:', fit$convergence_info$gauss_newton$n_converged, '\n')
  cat('  Improved:', fit$convergence_info$gauss_newton$n_improved, '\n')
}

# Compare parameter recovery by noise level
cat('\nParameter recovery by noise level:\n')
for (noise in unique(noise_levels)) {
  idx <- which(noise_levels == noise)
  param_err <- c(
    tau = mean(abs(fit$estimated_parameters[idx,'tau'] - true_tau[idx])),
    sigma = mean(abs(fit$estimated_parameters[idx,'sigma'] - true_sigma[idx])),
    rho = mean(abs(fit$estimated_parameters[idx,'rho'] - true_rho[idx]))
  )
  cat(sprintf('  Noise %.1f: tau=%.2f, sigma=%.2f, rho=%.2f\n', 
              noise, param_err['tau'], param_err['sigma'], param_err['rho']))
}

cat('\nStage 4 with aggressive refinement complete!\n')