# Example: Using fmrireg::baseline_model with estimate_parametric_hrf

library(fmriparametric)
library(fmrireg)

# Create sample fMRI data
set.seed(123)
sframe <- sampling_frame(blocklens = c(100, 100), TR = 2.0)
n_voxels <- 50

# Event model - simple block design
event_times <- seq(10, 90, by = 20)
event_onsets <- list(
  event_times,
  event_times + 100  # Second block
)
events <- event_model(
  data = data.frame(
    onset = unlist(event_onsets),
    duration = 10
  ),
  sframe = sframe
)

# True HRF parameters
true_tau <- 5.5
true_sigma <- 1.8
true_rho <- 0.35

# Generate data with drift and noise
time_points <- sframe$samples
y_data <- matrix(0, nrow = length(time_points), ncol = n_voxels)

# Add signal for each voxel
for (v in 1:n_voxels) {
  # Generate HRF
  hrf_times <- seq(0, 30, by = 0.1)
  hrf <- hrf_spmg1(hrf_times, tau = true_tau, sigma = true_sigma, rho = true_rho)
  
  # Convolve with events
  signal <- convolve_events(events, hrf, sframe)
  
  # Add polynomial drift
  drift <- 0.5 * poly(time_points, 3) %*% rnorm(3)
  
  # Add noise
  noise <- rnorm(length(time_points), sd = 0.5)
  
  y_data[, v] <- signal + drift + noise
}

# Example 1: Using baseline_model for polynomial drift
baseline_poly <- baseline_model(
  ~ poly(time, 3),
  sframe = sframe
)

fit_baseline <- estimate_parametric_hrf(
  fmri_data = y_data,
  event_model = events,
  baseline_model = baseline_poly,
  verbose = FALSE
)

# Example 2: Using baseline_model with splines
baseline_spline <- baseline_model(
  ~ ns(time, df = 5),
  sframe = sframe
)

fit_spline <- estimate_parametric_hrf(
  fmri_data = y_data,
  event_model = events,
  baseline_model = baseline_spline,
  verbose = FALSE
)

# Example 3: Using baseline_model with custom nuisance regressors
# (e.g., motion parameters)
motion_params <- matrix(rnorm(length(time_points) * 6, sd = 0.1), ncol = 6)
colnames(motion_params) <- paste0("motion", 1:6)

baseline_nuisance <- baseline_model(
  ~ poly(time, 3) + motion1 + motion2 + motion3 + motion4 + motion5 + motion6,
  sframe = sframe,
  data = as.data.frame(motion_params)
)

fit_nuisance <- estimate_parametric_hrf(
  fmri_data = y_data,
  event_model = events,
  baseline_model = baseline_nuisance,
  verbose = FALSE
)

# Compare results
cat("Parameter recovery with different baseline models:\n")
cat("True parameters: tau =", true_tau, ", sigma =", true_sigma, ", rho =", true_rho, "\n\n")

cat("Polynomial baseline:\n")
params_poly <- colMeans(coef(fit_baseline))
cat("  tau =", round(params_poly[1], 2), 
    ", sigma =", round(params_poly[2], 2),
    ", rho =", round(params_poly[3], 2), "\n")
cat("  Mean R² =", round(mean(fit_baseline$r_squared), 3), "\n\n")

cat("Spline baseline:\n")
params_spline <- colMeans(coef(fit_spline))
cat("  tau =", round(params_spline[1], 2),
    ", sigma =", round(params_spline[2], 2), 
    ", rho =", round(params_spline[3], 2), "\n")
cat("  Mean R² =", round(mean(fit_spline$r_squared), 3), "\n\n")

cat("Nuisance regressor baseline:\n")
params_nuisance <- colMeans(coef(fit_nuisance))
cat("  tau =", round(params_nuisance[1], 2),
    ", sigma =", round(params_nuisance[2], 2),
    ", rho =", round(params_nuisance[3], 2), "\n")
cat("  Mean R² =", round(mean(fit_nuisance$r_squared), 3), "\n")