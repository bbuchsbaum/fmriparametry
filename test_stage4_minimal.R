library(devtools)
load_all(quiet = TRUE)

# Test just the problematic part
set.seed(123)

# Minimal data
Y <- matrix(rnorm(50*2), 50, 2)
events <- matrix(0, 50, 1)
events[c(10, 30), 1] <- 1

# First test without any refinement
cat("Test 1: No refinement\n")
fit1 <- estimate_parametric_hrf(Y, events, 
                               tiered_refinement = "none",
                               global_refinement = FALSE,
                               verbose = FALSE)
cat("Success\n\n")

# Test with moderate only
cat("Test 2: Moderate refinement only\n")
fit2 <- estimate_parametric_hrf(Y, events, 
                               tiered_refinement = "moderate",
                               global_refinement = FALSE,
                               verbose = FALSE)
cat("Success\n\n")

# The problem case - but let's trace it
cat("Test 3: Aggressive refinement with tracing\n")

# Enable more debugging
options(error = function() {
  calls <- sys.calls()
  cat("\nCall stack:\n")
  for(i in length(calls):1) {
    cat(i, ": ", deparse(calls[[i]])[1], "\n", sep="")
  }
})

tryCatch({
  fit3 <- estimate_parametric_hrf(Y, events, 
                                 tiered_refinement = "aggressive",
                                 global_refinement = FALSE,
                                 verbose = TRUE)
  cat("Success\n")
}, error = function(e) {
  cat("\nError:", e$message, "\n")
  
  # Let's manually check what's happening
  cat("\nManual debugging:\n")
  
  # Prepare the data like Stage 0 would
  inputs <- .prepare_parametric_inputs(
    fmri_data = Y,
    event_model = events,
    hrf_eval_times = seq(0, 30, by = 0.1),
    hrf_span = 30
  )
  
  cat("Y_proj shape:", dim(inputs$Y_proj), "\n")
  cat("S_target shape:", dim(inputs$S_target), "\n")
  cat("S_target_proj shape:", dim(inputs$S_target_proj), "\n")
  
  # Check if S_target exists
  if (is.null(inputs$S_target)) {
    cat("WARNING: S_target is NULL!\n")
  }
})