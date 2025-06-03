#' Classify voxels into refinement tiers
#'
#' Classifies voxels into refinement tiers based on fit quality (R²) and 
#' uncertainty (standard errors) metrics. This determines which voxels need
#' additional refinement and what type of refinement to apply.
#'
#' @param r2_voxel Numeric vector of R-squared values for each voxel
#' @param se_theta_hat_voxel Matrix of standard errors (voxels x parameters) or NULL
#' @param refinement_opts List containing refinement thresholds and options
#'
#' @return List containing:
#'   - queue_labels: Character vector of queue assignments
#'   - queue_summary: Summary statistics for each queue
#'   - refinement_needed: Logical indicating if any refinement is needed
#' @keywords internal
.classify_refinement_queue <- function(
  r2_voxel,
  se_theta_hat_voxel = NULL,
  refinement_opts = list(
    apply_refinement = TRUE,
    r2_threshold_hard = 0.3,
    r2_threshold_moderate = 0.7,
    se_threshold_hard = 0.5,
    se_threshold_moderate = 0.3
  )
) {
  n_vox <- length(r2_voxel)
  
  # Initialize all voxels as "easy" (no refinement needed)
  queue_labels <- rep("easy", n_vox)
  
  if (!refinement_opts$apply_refinement) {
    return(list(
      queue_labels = queue_labels,
      queue_summary = table(queue_labels),
      refinement_needed = FALSE
    ))
  }
  
  # Extract thresholds
  r2_hard <- refinement_opts$r2_threshold_hard
  r2_moderate <- refinement_opts$r2_threshold_moderate
  se_hard <- refinement_opts$se_threshold_hard
  se_moderate <- refinement_opts$se_threshold_moderate
  
  # Calculate SE metrics if available
  se_metric <- NULL
  if (!is.null(se_theta_hat_voxel)) {
    # Use average absolute SE across parameters as the metric
    # se_theta_hat_voxel is expected to contain standard errors for each
    # parameter, so we take the mean across parameters per voxel
    se_metric <- rowMeans(se_theta_hat_voxel, na.rm = TRUE)
    
    # Cap extreme values
    se_metric[is.na(se_metric)] <- Inf
    se_metric[is.infinite(se_metric)] <- se_hard * 2
  }
  
  # Classification logic
  # Priority: hard > moderate > easy
  
  # First identify hard cases (poor R² OR high uncertainty)
  hard_r2 <- r2_voxel < r2_hard
  hard_se <- if (!is.null(se_metric)) se_metric > se_hard else rep(FALSE, n_vox)
  hard_cases <- hard_r2 | hard_se
  
  # Then identify moderate cases (moderate R² OR moderate uncertainty)
  # But not already classified as hard
  moderate_r2 <- r2_voxel >= r2_hard & r2_voxel < r2_moderate
  moderate_se <- if (!is.null(se_metric)) {
    se_metric > se_moderate & se_metric <= se_hard
  } else {
    rep(FALSE, n_vox)
  }
  moderate_cases <- (moderate_r2 | moderate_se) & !hard_cases
  
  # Assign queue labels
  queue_labels[hard_cases] <- "hard_GN"
  queue_labels[moderate_cases] <- "moderate_local_recenter"
  
  # Handle missing data
  # Voxels with NA R² are automatically hard cases
  na_voxels <- is.na(r2_voxel)
  if (any(na_voxels)) {
    queue_labels[na_voxels] <- "hard_GN"
  }
  
  # Create summary
  queue_summary <- table(queue_labels)
  queue_proportions <- prop.table(queue_summary)
  
  # Detailed summary by queue
  queue_details <- list()
  for (q in unique(queue_labels)) {
    idx <- which(queue_labels == q)
    queue_details[[q]] <- list(
      n = length(idx),
      proportion = length(idx) / n_vox,
      r2_mean = mean(r2_voxel[idx], na.rm = TRUE),
      r2_median = median(r2_voxel[idx], na.rm = TRUE),
      r2_range = range(r2_voxel[idx], na.rm = TRUE)
    )
    
    if (!is.null(se_metric)) {
      queue_details[[q]]$se_mean <- mean(se_metric[idx], na.rm = TRUE)
      queue_details[[q]]$se_median <- median(se_metric[idx], na.rm = TRUE)
    }
  }
  
  # Determine if refinement is needed
  refinement_needed <- any(queue_labels != "easy")
  
  # Return classification results
  list(
    queue_labels = queue_labels,
    queue_summary = queue_summary,
    queue_proportions = queue_proportions,
    queue_details = queue_details,
    refinement_needed = refinement_needed,
    classification_criteria = list(
      r2_thresholds = c(hard = r2_hard, moderate = r2_moderate),
      se_thresholds = if (!is.null(se_metric)) {
        c(hard = se_hard, moderate = se_moderate)
      } else NULL,
      se_available = !is.null(se_metric)
    )
  )
}

#' Print refinement queue summary
#'
#' Helper function to print a formatted summary of the refinement queue
#' classification results.
#'
#' @param queue_result Result from .classify_refinement_queue
#' @param verbose Logical whether to print
#' @keywords internal
.print_refinement_summary <- function(queue_result, verbose = TRUE) {
  if (!verbose) return(invisible(NULL))
  
  cat("\nRefinement Queue Classification:\n")
  cat("================================\n")
  
  # Overall summary
  cat("Total voxels:", sum(queue_result$queue_summary), "\n")
  cat("Refinement needed:", ifelse(queue_result$refinement_needed, "Yes", "No"), "\n\n")
  
  # Queue breakdown
  cat("Queue assignments:\n")
  for (q in names(queue_result$queue_summary)) {
    n <- queue_result$queue_summary[q]
    prop <- queue_result$queue_proportions[q]
    details <- queue_result$queue_details[[q]]
    
    cat(sprintf("  %-25s: %6d voxels (%5.1f%%) | Mean R² = %.3f\n",
                q, n, prop * 100, details$r2_mean))
  }
  
  # Classification criteria
  cat("\nClassification criteria:\n")
  crit <- queue_result$classification_criteria
  cat("  R² thresholds: hard <", crit$r2_thresholds["hard"], 
      ", moderate <", crit$r2_thresholds["moderate"], "\n")
  
  if (!is.null(crit$se_thresholds)) {
    cat("  SE thresholds: hard >", crit$se_thresholds["hard"],
        ", moderate >", crit$se_thresholds["moderate"], "\n")
  } else {
    cat("  SE thresholds: not used (SEs not available)\n")
  }
  
  cat("\n")
  invisible(NULL)
}