library(devtools)
load_all(quiet = TRUE)

# Simple test to debug the issue
set.seed(123)
Y <- matrix(rnorm(100*5), 100, 5)
events <- matrix(0, 100, 1)
events[c(20, 50, 80), 1] <- 1

# Test without tiered refinement first
cat("Test 1: No refinement\n")
fit1 <- estimate_parametric_hrf(Y, events, tiered_refinement = "none", verbose = FALSE)
cat("Success - Mean R²:", mean(fit1$r_squared), "\n\n")

# Test with moderate refinement
cat("Test 2: Moderate refinement\n")
fit2 <- estimate_parametric_hrf(Y, events, tiered_refinement = "moderate", verbose = FALSE)
cat("Success - Mean R²:", mean(fit2$r_squared), "\n")
cat("Refinement info:", names(fit2$refinement_info), "\n\n")

# Test with aggressive refinement (the failing case)
cat("Test 3: Aggressive refinement\n")
tryCatch({
  fit3 <- estimate_parametric_hrf(Y, events, tiered_refinement = "aggressive", verbose = TRUE)
  cat("Success - Mean R²:", mean(fit3$r_squared), "\n")
}, error = function(e) {
  cat("Error occurred:", e$message, "\n")
  # Try to debug by looking at the gauss-newton input
  cat("\nDebugging: let's check what's happening in Stage 4...\n")
})

# Let's check if the issue is with the default theta_bounds
hrf_interface <- .get_hrf_interface("lwu")
default_bounds <- hrf_interface$default_bounds()
cat("\nDefault bounds:\n")
print(default_bounds)