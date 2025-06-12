# Debug refactored implementation
library(devtools)
load_all()

# Minimal test
set.seed(42)
Y <- matrix(rnorm(50 * 2), 50, 2)
events <- matrix(0, 50, 1)
events[c(10, 30), 1] <- 1

# Try with debug enabled
result <- tryCatch({
  estimate_parametric_hrf(
    fmri_data = Y,
    event_model = events,
    parametric_hrf = "lwu",
    global_refinement = FALSE,
    tiered_refinement = "none",
    compute_se = FALSE,
    verbose = TRUE,
    .implementation = "refactored"
  )
}, error = function(e) {
  cat("Error occurred:", e$message, "\n")
  traceback()
  NULL
})