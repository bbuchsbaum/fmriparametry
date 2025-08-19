#' Local voxel-specific re-centering for moderate difficulty voxels
#'
#' Performs local re-centering using each voxel's current estimate as its own
#' expansion point. This is more computationally intensive than global or K-means
#' re-centering but can improve fits for voxels that don't fit well with any
#' global or cluster-based expansion point.
#'
#' @param theta_hat_voxel Matrix of current parameter estimates (voxels x parameters)
#' @param r2_voxel Numeric vector of current R-squared values
#' @param Y_proj Numeric matrix of projected BOLD data (timepoints x voxels)
#' @param S_target_proj Numeric matrix of projected stimulus design
#' @param scan_times Numeric vector of scan acquisition times
#' @param hrf_eval_times Numeric vector of HRF evaluation time points
#' @param hrf_interface List with HRF interface functions
#' @param theta_bounds List with elements `lower` and `upper`
#' @param queue_labels Character vector of refinement queue assignments
#' @param coeffs_voxel Matrix of current linear coefficients
#' @param lambda_ridge Ridge penalty for refitting
#' @param epsilon_beta Small value to avoid division by zero
#' @param verbose Logical whether to print progress
#'
#' @return List with updated estimates and refinement statistics
#' @noRd
.local_recentering_moderate <- function(
  theta_hat_voxel,
  r2_voxel,
  Y_proj,
  S_target_proj,
  scan_times,
  hrf_eval_times,
  hrf_interface,
  theta_bounds,
  queue_labels,
  coeffs_voxel,
  lambda_ridge = 0.01,
  epsilon_beta = 1e-6,
  verbose = FALSE
) {
  n_vox <- nrow(theta_hat_voxel)
  n_params <- ncol(theta_hat_voxel)
  n_time <- nrow(Y_proj)
  
  # Identify moderate voxels
  idx_moderate <- which(queue_labels == "moderate_local_recenter")
  n_moderate <- length(idx_moderate)
  
  if (n_moderate == 0) {
    if (verbose) cat("No moderate voxels to refine\n")
    return(list(
      theta_hat = theta_hat_voxel,
      r2 = r2_voxel,
      coeffs = coeffs_voxel,
      n_refined = 0,
      n_improved = 0
    ))
  }
  
  if (verbose) {
    cat("\nLocal re-centering for", n_moderate, "moderate voxels\n")
  }
  
  # Store original values for comparison
  theta_hat_orig <- theta_hat_voxel
  r2_orig <- r2_voxel
  coeffs_orig <- coeffs_voxel
  
  # Track improvements
  n_improved <- 0
  improvement_details <- numeric(n_moderate)
  
  # Process each moderate voxel individually
  for (i in seq_along(idx_moderate)) {
    v <- idx_moderate[i]
    
    # Use voxel's current estimate as expansion point
    theta_v <- theta_hat_voxel[v, ]
    
    # Skip if parameters are at bounds (likely stuck)
    if (any(theta_v == theta_bounds$lower) || any(theta_v == theta_bounds$upper)) {
      next
    }
    
    # Construct Taylor basis for this voxel
    X_taylor_v <- hrf_interface$taylor_basis(theta_v, hrf_eval_times)
    if (!is.matrix(X_taylor_v)) {
      X_taylor_v <- matrix(X_taylor_v, ncol = n_params + 1)
    }
    
    # Design matrix via convolution using centralized helper
    X_design_v <- .convolve_signal_with_kernels(S_target_proj[, 1], X_taylor_v, n_time)
    
    # Solve for this voxel
    tryCatch({
      # QR decomposition
      qr_decomp <- qr(X_design_v)
      Q <- qr.Q(qr_decomp)
      R <- qr.R(qr_decomp)
      R_inv <- solve(R + lambda_ridge * diag(ncol(R)))
      
      # Estimate coefficients
      y_v <- Y_proj[, v]
      coeffs_v_new <- R_inv %*% t(Q) %*% y_v
      
      # Extract parameters
      beta0_new <- coeffs_v_new[1]
      if (abs(beta0_new) < epsilon_beta) {
        # Skip if amplitude too small
        next
      }
      
      delta_theta <- coeffs_v_new[2:(n_params+1)] / beta0_new
      theta_v_new <- theta_v + delta_theta
      
      # Apply bounds
      theta_v_new <- pmax(theta_bounds$lower, pmin(theta_v_new, theta_bounds$upper))
      
      # Calculate new R^2
      fitted_v_new <- X_design_v %*% coeffs_v_new
      resid_v_new <- y_v - fitted_v_new
      ss_res_new <- sum(resid_v_new^2)
      ss_tot <- sum((y_v - mean(y_v))^2)
      r2_v_new <- if (ss_tot > 0) 1 - ss_res_new / ss_tot else 0
      
      # Update only if improved
      if (r2_v_new > r2_voxel[v]) {
        theta_hat_voxel[v, ] <- theta_v_new
        r2_voxel[v] <- r2_v_new
        coeffs_voxel[, v] <- coeffs_v_new
        n_improved <- n_improved + 1
        improvement_details[i] <- r2_v_new - r2_orig[v]
      }
      
    }, error = function(e) {
      # Skip voxel if numerical issues
      if (verbose) cat("  Warning: Local recentering failed for voxel", v, "\n")
    })
    
    # Progress reporting
    if (verbose && i %% 100 == 0) {
      cat("  Processed", i, "/", n_moderate, "voxels,",
          n_improved, "improved\n")
    }
  }
  
  if (verbose) {
    cat("  Local re-centering complete:\n")
    cat("    Refined:", n_moderate, "voxels\n")
    cat("    Improved:", n_improved, "voxels (", 
        round(100 * n_improved / n_moderate, 1), "%)\n")
    if (n_improved > 0) {
      cat("    Mean R^2 improvement:", 
          round(mean(improvement_details[improvement_details > 0]), 4), "\n")
    }
  }
  
  # Update queue labels for successfully refined voxels
  successfully_refined <- idx_moderate[which(r2_voxel[idx_moderate] > r2_orig[idx_moderate])]
  if (length(successfully_refined) > 0) {
    # Move successfully refined moderate voxels to easy queue
    queue_labels[successfully_refined] <- "easy"
  }
  
  # Return results
  list(
    theta_hat = theta_hat_voxel,
    r2 = r2_voxel,
    coeffs = coeffs_voxel,
    queue_labels = queue_labels,
    n_refined = n_moderate,
    n_improved = n_improved,
    prop_improved = n_improved / n_moderate,
    improvement_summary = list(
      mean_improvement = if (n_improved > 0) {
        mean(improvement_details[improvement_details > 0])
      } else 0,
      max_improvement = if (n_improved > 0) {
        max(improvement_details)
      } else 0
    )
  )
}