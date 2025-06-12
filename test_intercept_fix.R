#!/usr/bin/env Rscript

# Test script to verify the intercept fix for R² = 0 issue

library(fmriparametric)

# Set debug mode
options(fmriparametric.debug = TRUE)

# Create test data with baseline
set.seed(42)
n_time <- 100
n_vox <- 2

# Create stimulus with 5 events
stim <- rep(0, n_time)
stim[c(10, 30, 50, 70, 90)] <- 1
S <- matrix(stim, ncol = 1)

# Create synthetic HRF response with baseline
hrf <- c(0, 0.2, 0.8, 1, 0.7, 0.3, 0, -0.1, -0.05, rep(0, 20))
expected_signal <- convolve(stim, rev(hrf), type = "open")[1:n_time]

# Create data with baseline of 100 and noise
baseline <- 100
Y <- matrix(baseline + 10 * expected_signal + rnorm(n_time * n_vox, sd = 2), 
            nrow = n_time, ncol = n_vox)

cat("Data properties:\n")
cat("Y mean:", mean(Y), "\n")
cat("Y range:", range(Y), "\n")
cat("Expected signal range:", range(expected_signal), "\n\n")

# Test with default intercept baseline
cat("Testing with baseline_model = 'intercept' (default):\n")
fit <- estimate_parametric_hrf(
  Y, S, 
  parametric_hrf = "lwu",
  verbose = FALSE,
  global_refinement = FALSE
)

cat("\nResults:\n")
cat("R-squared values:", fit$fit_quality$r_squared, "\n")
cat("Mean R-squared:", mean(fit$fit_quality$r_squared), "\n")

# Test without baseline (should fail)
cat("\n\nTesting with baseline_model = NULL:\n")
fit_no_baseline <- estimate_parametric_hrf(
  Y, S, 
  parametric_hrf = "lwu",
  baseline_model = NULL,
  verbose = FALSE,
  global_refinement = FALSE
)

cat("\nResults without baseline:\n")
cat("R-squared values:", fit_no_baseline$fit_quality$r_squared, "\n")
cat("Mean R-squared:", mean(fit_no_baseline$fit_quality$r_squared), "\n")

cat("\n\nSUMMARY:\n")
cat("With intercept: Mean R² =", mean(fit$fit_quality$r_squared), "\n")
cat("Without intercept: Mean R² =", mean(fit_no_baseline$fit_quality$r_squared), "\n")
cat("Improvement:", mean(fit$fit_quality$r_squared) - mean(fit_no_baseline$fit_quality$r_squared), "\n")