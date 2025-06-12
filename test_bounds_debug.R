# Debug bounds enforcement
library(fmriparametric)

set.seed(123)

# Create test data
n_time <- 100
n_voxels <- 5
fmri_data <- matrix(rnorm(n_time * n_voxels), nrow = n_time, ncol = n_voxels)

# Create event design with some events
event_design <- matrix(0, nrow = n_time, ncol = 1)
event_design[seq(10, 90, by = 20), 1] <- 1

# Test with restrictive bounds
restrictive_bounds <- list(
  lower = c(2, 1, 0.1),    # tau >= 2, sigma >= 1, rho >= 0.1
  upper = c(8, 3, 0.5)     # tau <= 8, sigma <= 3, rho <= 0.5
)

# Fit with restrictive bounds
fit_restricted <- estimate_parametric_hrf(
  fmri_data = fmri_data,
  event_model = event_design,
  theta_bounds = restrictive_bounds,
  verbose = TRUE
)

# Check all parameters
params <- fit_restricted$estimated_parameters
print("Estimated parameters:")
print(params)

print("\nBounds:")
print(restrictive_bounds)

print("\nViolations:")
print(paste("tau < lower:", sum(params[, 1] < restrictive_bounds$lower[1])))
print(paste("tau > upper:", sum(params[, 1] > restrictive_bounds$upper[1])))
print(paste("sigma < lower:", sum(params[, 2] < restrictive_bounds$lower[2])))
print(paste("sigma > upper:", sum(params[, 2] > restrictive_bounds$upper[2])))
print(paste("rho < lower:", sum(params[, 3] < restrictive_bounds$lower[3])))
print(paste("rho > upper:", sum(params[, 3] > restrictive_bounds$upper[3])))