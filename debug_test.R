# Debug the matrix issue
# =====================

# Load our function
source("R/estimate_parametric_hrf.R")

# Try a simple test case
set.seed(42)
n_time <- 50
n_vox <- 10

fmri_data <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
# Create a simple stimulus timing vector instead of design matrix
event_model <- matrix(c(rep(c(1,0,0,0,0), 10)), n_time, 1)

cat("Input dimensions:\n")
cat("  fmri_data:", dim(fmri_data), "\n")
cat("  event_model:", dim(event_model), "\n")

# Try to call the function and see where it fails
result <- tryCatch({
  estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    verbose = FALSE
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  traceback()
  return(NULL)
})

if (!is.null(result)) {
  cat("SUCCESS! Result has:\n")
  cat("  estimated_parameters:", class(result$estimated_parameters), dim(result$estimated_parameters), "\n")
} else {
  cat("FAILED - debugging needed\n")
}