#' Simple voxel classification for refinement
#'
#' Classifies voxels into refinement categories based solely on R-squared values.
#' This replaces the complex queue classification system with a straightforward
#' approach that matches actual usage.
#'
#' @param r_squared Numeric vector of R-squared values
#' @param r2_threshold_easy Threshold above which voxels are considered "easy" (default: 0.7)
#' @param r2_threshold_hard Threshold below which voxels are considered "hard" (default: 0.3)
#'
#' @return Character vector with classifications: "easy", "moderate", or "hard"
#' @noRd
.classify_voxels_simple <- function(r_squared, 
                                   r2_threshold_easy = 0.7,
                                   r2_threshold_hard = 0.3) {
  # Validate inputs
  if (!is.numeric(r_squared)) {
    stop("r_squared must be numeric")
  }
  
  # Initialize all as moderate
  classes <- rep("moderate", length(r_squared))
  
  # Apply thresholds
  classes[r_squared > r2_threshold_easy] <- "easy"
  classes[r_squared < r2_threshold_hard] <- "hard"
  
  # Handle missing values
  classes[is.na(r_squared)] <- "hard"
  
  return(classes)
}

#' Get refinement summary
#'
#' Provides a summary of voxel classifications for refinement
#'
#' @param classifications Character vector of voxel classifications
#' @param r_squared Numeric vector of R-squared values
#'
#' @return List with summary statistics
#' @noRd
.get_refinement_summary <- function(classifications, r_squared) {
  # Count by class
  class_counts <- table(classifications)
  
  # R^2 statistics by class
  r2_by_class <- tapply(r_squared, classifications, function(x) {
    c(mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE))
  })
  
  list(
    counts = as.list(class_counts),
    proportions = as.list(prop.table(class_counts)),
    r2_stats = r2_by_class,
    total_voxels = length(classifications),
    refinement_needed = any(classifications != "easy")
  )
}