# Demonstration of robust parameter validation
# This demo shows how the package validates inputs and handles edge cases

library(fmriparametric)

cat("=== Robust Parameter Validation Demo ===\n\n")

# Create some test data
set.seed(123)
n_time <- 100
n_vox <- 10

# Generate simple fMRI data
fmri_data <- matrix(rnorm(n_time * n_vox), nrow = n_time)

# Create event model with a few events
events <- matrix(0, nrow = n_time, ncol = 1)
events[c(10, 30, 50, 70, 90), 1] <- 1

cat("1. Testing with valid inputs:\n")
tryCatch({
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = events,
    parametric_model = "lwu",
    verbose = FALSE
  )
  cat("   SUCCESS: Valid inputs accepted\n")
  print(summary(fit))
}, error = function(e) {
  cat("   ERROR:", e$message, "\n")
})

cat("\n2. Testing with invalid fmri_data:\n")
tryCatch({
  fit <- estimate_parametric_hrf(
    fmri_data = "not_a_matrix",
    event_model = events,
    verbose = FALSE
  )
}, error = function(e) {
  cat("   EXPECTED ERROR:", e$message, "\n")
})

cat("\n3. Testing with mismatched dimensions:\n")
tryCatch({
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = events[1:50, , drop = FALSE],  # Wrong length
    verbose = FALSE
  )
}, error = function(e) {
  cat("   EXPECTED ERROR:", e$message, "\n")
})

cat("\n4. Testing with custom bounds:\n")
tryCatch({
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = events,
    theta_bounds = list(
      lower = c(2, 0.5, 0),
      upper = c(8, 4, 1)
    ),
    verbose = FALSE
  )
  cat("   SUCCESS: Custom bounds accepted\n")
}, error = function(e) {
  cat("   ERROR:", e$message, "\n")
})

cat("\n=== Demo Complete ===\n")