# Test the new default behavior
library(devtools)
load_all()

# Simple test data
set.seed(42)
Y <- matrix(rnorm(50 * 2), 50, 2)
events <- matrix(0, 50, 1)
events[c(10, 30), 1] <- 1

cat("=== Testing new default behavior ===\n\n")

# Test 1: Default should use refactored
cat("Test 1: Default (should use refactored)\n")
result_default <- estimate_parametric_hrf(
  fmri_data = Y,
  event_model = events,
  verbose = FALSE
)
cat("Success - R²:", round(result_default$r_squared, 3), "\n\n")

# Test 2: Explicit legacy should show warning
cat("Test 2: Explicit legacy (should show deprecation warning)\n")
result_legacy <- estimate_parametric_hrf(
  fmri_data = Y,
  event_model = events,
  verbose = FALSE,
  .implementation = "legacy"
)
cat("Success - R²:", round(result_legacy$r_squared, 3), "\n\n")

# Test 3: Compare option should fail
cat("Test 3: Compare option (should fail)\n")
tryCatch({
  estimate_parametric_hrf(
    fmri_data = Y,
    event_model = events,
    verbose = FALSE,
    .implementation = "compare"
  )
  cat("UNEXPECTED: Compare didn't fail!\n")
}, error = function(e) {
  cat("Expected error:", e$message, "\n\n")
})

# Test 4: Global option for legacy
cat("Test 4: Global option for legacy\n")
options(fmriparametric.use_legacy = TRUE)
result_global <- estimate_parametric_hrf(
  fmri_data = Y,
  event_model = events,
  verbose = TRUE
)
cat("Success - R²:", round(result_global$r_squared, 3), "\n\n")

# Clean up
options(fmriparametric.use_legacy = NULL)

cat("All tests completed successfully!\n")