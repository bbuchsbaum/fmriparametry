# Debug test for R-squared recalculation

library(testthat)
library(fmriparametric)

test_that("Debug R-squared recalculation", {
  skip_if_not_installed("fmrireg")
  
  # Minimal test case
  set.seed(123)
  n_time <- 50
  n_vox <- 2
  
  # Simple design matrix with 3 events
  design_matrix <- matrix(0, n_time, 1)
  design_matrix[c(10, 20, 30), 1] <- 1
  
  # Create data with clear signal
  Y <- matrix(rnorm(n_time * n_vox, 100, 5), n_time, n_vox)
  
  # Add strong signal that should be detectable
  for (i in c(10, 20, 30)) {
    Y[i:(i+5), ] <- Y[i:(i+5), ] + 50
  }
  
  # Run with refinement enabled
  result <- estimate_parametric_hrf(
    fmri_data = Y,
    event_model = design_matrix,
    parametric_hrf = "lwu",
    global_refinement = TRUE,
    tiered_refinement = "moderate",
    verbose = TRUE
  )
  
  # Print diagnostics
  cat("\n=== DEBUG INFO ===\n")
  cat("R-squared values:", result$r_squared, "\n")
  cat("Mean R-squared:", mean(result$r_squared), "\n")
  cat("Refinement info exists:", !is.null(result$refinement_info), "\n")
  
  if (!is.null(result$refinement_info)) {
    cat("Refinement stages:\n")
    print(names(result$refinement_info))
    
    if (!is.null(result$refinement_info$moderate)) {
      cat("Moderate refinement:\n")
      cat("  - Refined indices:", result$refinement_info$moderate$refined_indices, "\n")
      cat("  - Initial theta:\n")
      print(result$refinement_info$moderate$initial_theta)
      cat("  - Refined theta:\n") 
      print(result$refinement_info$moderate$refined_theta)
    }
  }
  
  # The R-squared should be > 0 with this strong signal
  expect_true(mean(result$r_squared) > 0.1)
})