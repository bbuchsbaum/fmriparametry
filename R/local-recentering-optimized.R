#' Optimized local voxel-specific re-centering for moderate difficulty voxels
#'
#' Performs local re-centering using each voxel's current estimate as its own
#' expansion point. This version is optimized to process voxels in batches
#' when they share similar parameters, reducing computational overhead.
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
#' @param batch_tolerance Tolerance for grouping similar parameters (default: 1e-6)
#' @param verbose Logical whether to print progress
#'
#' @return List with updated estimates and refinement statistics
#' @noRd
.local_recentering_moderate_optimized <- function(
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
  batch_tolerance = 1e-6,
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
    cat("\nOptimized local re-centering for", n_moderate, "moderate voxels\n")
  }
  
  # Store original values for comparison
  theta_hat_orig <- theta_hat_voxel
  r2_orig <- r2_voxel
  coeffs_orig <- coeffs_voxel
  
  # Track improvements
  n_improved <- 0
  improvement_details <- numeric(n_moderate)
  
  # Extract parameters for moderate voxels
  theta_moderate <- theta_hat_voxel[idx_moderate, , drop = FALSE]
  
  # Group voxels with similar parameters (within tolerance)
  # This allows batch processing of voxels with nearly identical parameters
  group_ids <- integer(n_moderate)
  n_groups <- 0
  
  for (i in seq_len(n_moderate)) {
    if (group_ids[i] == 0) {
      n_groups <- n_groups + 1
      group_ids[i] <- n_groups
      
      # Find other voxels with similar parameters
      if (i < n_moderate) {
        for (j in (i + 1):n_moderate) {
          if (group_ids[j] == 0) {
            param_diff <- max(abs(theta_moderate[i, ] - theta_moderate[j, ]))
            if (param_diff < batch_tolerance) {
              group_ids[j] <- n_groups
            }
          }
        }
      }
    }
  }
  
  if (verbose) {
    cat("  Grouped into", n_groups, "batches for processing\n")
  }
  
  # Process each group
  for (g in seq_len(n_groups)) {
    group_voxels <- which(group_ids == g)
    group_idx <- idx_moderate[group_voxels]
    n_group <- length(group_idx)
    
    # Use the median parameters of the group as expansion point
    theta_group <- apply(theta_moderate[group_voxels, , drop = FALSE], 2, median)
    
    # Skip if parameters are at bounds (likely stuck)
    if (any(theta_group <= theta_bounds$lower + 1e-6) || 
        any(theta_group >= theta_bounds$upper - 1e-6)) {
      next
    }
    
    # Construct Taylor basis for this group
    X_taylor <- hrf_interface$taylor_basis(theta_group, hrf_eval_times)
    if (!is.matrix(X_taylor)) {
      X_taylor <- matrix(X_taylor, ncol = n_params + 1)
    }
    
    # Design matrix via convolution using centralized helper
    X_design <- .convolve_signal_with_kernels(S_target_proj[, 1], X_taylor, n_time)
    
    # Batch solve for all voxels in group using QR decomposition
    tryCatch({
      # QR decomposition
      qr_decomp <- qr(X_design)
      Q <- qr.Q(qr_decomp)
      R <- qr.R(qr_decomp)
      R_ridge <- R + lambda_ridge * diag(ncol(R))
      
      # Extract data for this group
      Y_group <- Y_proj[, group_idx, drop = FALSE]
      
      # Batch solve
      coeffs_new <- solve(R_ridge, t(Q) %*% Y_group)
      
      # Process results for each voxel in group
      for (k in seq_len(n_group)) {
        v <- group_idx[k]
        local_idx <- which(idx_moderate == v)
        
        # Extract coefficients for this voxel
        coeffs_v <- coeffs_new[, k]
        beta0_new <- coeffs_v[1]
        
        if (abs(beta0_new) < epsilon_beta) {
          # Skip if amplitude too small
          next
        }
        
        # Compute parameter update
        delta_theta <- coeffs_v[2:(n_params+1)] / beta0_new
        theta_v_new <- theta_group + delta_theta
        
        # Apply bounds
        theta_v_new <- pmax(theta_bounds$lower, pmin(theta_v_new, theta_bounds$upper))
        
        # Calculate new R^2
        fitted_v <- X_design %*% coeffs_v
        resid_v <- Y_group[, k] - fitted_v
        ss_res_new <- sum(resid_v^2)
        ss_tot <- sum((Y_group[, k] - mean(Y_group[, k]))^2)
        r2_v_new <- if (ss_tot > 0) 1 - ss_res_new / ss_tot else 0
        
        # Update only if improved
        if (r2_v_new > r2_voxel[v]) {
          theta_hat_voxel[v, ] <- theta_v_new
          r2_voxel[v] <- r2_v_new
          coeffs_voxel[, v] <- coeffs_v
          n_improved <- n_improved + 1
          improvement_details[local_idx] <- r2_v_new - r2_orig[v]
        }
      }
      
    }, error = function(e) {
      # Skip group if numerical issues
      if (verbose) cat("  Warning: Local recentering failed for group", g, "\n")
    })
    
    # Progress reporting
    if (verbose && g %% 10 == 0) {
      cat("  Processed", g, "/", n_groups, "groups,",
          n_improved, "voxels improved\n")
    }
  }
  
  if (verbose) {
    cat("  Optimized local re-centering complete:\n")
    cat("    Refined:", n_moderate, "voxels in", n_groups, "batches\n")
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
    n_groups = n_groups,
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