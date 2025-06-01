#' K-means based spatial re-centering for heterogeneous data
#'
#' Implements K-means clustering on parameter estimates to identify spatial 
#' clusters with distinct HRF characteristics, then re-centers the Taylor
#' expansion within each cluster for improved fits.
#'
#' @param theta_hat_voxel Matrix of current parameter estimates (voxels x parameters)
#' @param r2_voxel Numeric vector of R-squared values for each voxel
#' @param Y_proj Numeric matrix of projected BOLD data (timepoints x voxels)
#' @param S_target_proj Numeric matrix of projected stimulus design
#' @param scan_times Numeric vector of scan acquisition times
#' @param hrf_eval_times Numeric vector of HRF evaluation time points
#' @param hrf_interface List with HRF interface functions
#' @param theta_bounds List with elements `lower` and `upper`
#' @param coeffs_voxel Matrix of current linear coefficients
#' @param recenter_kmeans_passes Integer number of K-means iterations
#' @param kmeans_k Integer number of clusters
#' @param r2_threshold_kmeans Numeric R² threshold for selecting voxels for clustering
#' @param transform_params Logical whether to transform parameters for clustering
#' @param lambda_ridge Ridge penalty for refitting
#' @param epsilon_beta Small value to avoid division by zero
#' @param verbose Logical whether to print progress
#'
#' @return List with updated estimates and clustering information
#' @keywords internal
.kmeans_recentering <- function(
  theta_hat_voxel,
  r2_voxel,
  Y_proj,
  S_target_proj,
  scan_times,
  hrf_eval_times,
  hrf_interface,
  theta_bounds,
  coeffs_voxel,
  recenter_kmeans_passes = 2,
  kmeans_k = 5,
  r2_threshold_kmeans = 0.7,
  transform_params = TRUE,
  lambda_ridge = 0.01,
  epsilon_beta = 1e-6,
  verbose = FALSE
) {
  n_vox <- nrow(theta_hat_voxel)
  n_params <- ncol(theta_hat_voxel)
  n_time <- nrow(Y_proj)
  
  # Store original values
  theta_hat_orig <- theta_hat_voxel
  coeffs_orig <- coeffs_voxel
  r2_orig <- r2_voxel
  
  # Track clustering information
  cluster_info <- list(
    n_passes = 0,
    cluster_assignments = NULL,
    cluster_centers = list(),
    cluster_r2_improvements = list()
  )
  
  if (verbose) cat("\nStarting K-means re-centering with k =", kmeans_k, "\n")
  
  for (pass_km in seq_len(recenter_kmeans_passes)) {
    if (verbose) cat("K-means pass", pass_km, "\n")
    
    # Select good voxels for clustering
    idx_good <- which(r2_voxel >= r2_threshold_kmeans)
    
    if (length(idx_good) < kmeans_k * 10) {
      if (verbose) cat("  Too few good voxels (", length(idx_good), 
                       ") for meaningful clustering\n")
      break
    }
    
    # Prepare parameters for clustering
    theta_for_clustering <- theta_hat_voxel[idx_good, , drop = FALSE]
    
    if (transform_params) {
      # Transform parameters to improve clustering
      # For LWU model: tau is already in good scale, log-transform sigma, logit rho
      theta_transformed <- theta_for_clustering
      if (n_params >= 2) {
        # Log transform width parameter (ensure positive)
        theta_transformed[, 2] <- log(pmax(theta_transformed[, 2], 0.1))
      }
      if (n_params >= 3) {
        # Logit transform for bounded parameter (0-1.5 range)
        rho_scaled <- pmin(pmax(theta_transformed[, 3], 0.001), 1.499) / 1.5
        theta_transformed[, 3] <- log(rho_scaled / (1 - rho_scaled))
      }
    } else {
      theta_transformed <- theta_for_clustering
    }
    
    # Standardize for clustering
    theta_scaled <- scale(theta_transformed)
    
    # Perform K-means with multiple starts
    set.seed(123 + pass_km)  # For reproducibility
    km_result <- tryCatch({
      kmeans(theta_scaled, centers = kmeans_k, nstart = 20, iter.max = 50)
    }, error = function(e) {
      if (verbose) cat("  K-means failed:", e$message, "\n")
      NULL
    })
    
    if (is.null(km_result)) break
    
    # Get cluster centers in original parameter space
    centers_scaled <- km_result$centers
    centers_transformed <- sweep(centers_scaled, 2, 
                                 attr(theta_scaled, "scaled:scale"), "*")
    centers_transformed <- sweep(centers_transformed, 2, 
                                 attr(theta_scaled, "scaled:center"), "+")
    
    # Back-transform if necessary
    cluster_centers <- centers_transformed
    if (transform_params) {
      if (n_params >= 2) {
        cluster_centers[, 2] <- exp(cluster_centers[, 2])
      }
      if (n_params >= 3) {
        # Inverse logit
        cluster_centers[, 3] <- 1.5 * plogis(cluster_centers[, 3])
      }
    }
    
    # Ensure cluster centers are within bounds
    for (k in seq_len(kmeans_k)) {
      cluster_centers[k, ] <- pmax(theta_bounds$lower, 
                                   pmin(cluster_centers[k, ], theta_bounds$upper))
    }
    
    if (verbose) {
      cat("  Cluster sizes:", table(km_result$cluster), "\n")
      cat("  Cluster centers:\n")
      print(round(cluster_centers, 3))
    }
    
    # Create full cluster assignment vector
    cluster_assignment_full <- rep(NA, n_vox)
    cluster_assignment_full[idx_good] <- km_result$cluster
    
    # Process each cluster
    cluster_improvements <- numeric(kmeans_k)
    
    for (k in seq_len(kmeans_k)) {
      # Find voxels in this cluster
      idx_cluster <- which(cluster_assignment_full == k)
      
      if (length(idx_cluster) == 0) next
      
      if (verbose) cat("  Processing cluster", k, "with", 
                       length(idx_cluster), "voxels\n")
      
      # Use cluster center as expansion point
      theta_cluster <- cluster_centers[k, ]
      
      # Construct Taylor basis at cluster center
      X_taylor <- hrf_interface$taylor_basis(theta_cluster, hrf_eval_times)
      if (!is.matrix(X_taylor)) {
        X_taylor <- matrix(X_taylor, ncol = n_params + 1)
      }
      
      # Design matrix via convolution
      X_design <- matrix(0, nrow = n_time, ncol = ncol(X_taylor))
      for (j in seq_len(ncol(X_taylor))) {
        basis_col <- X_taylor[, j]
        conv_full <- stats::convolve(S_target_proj[, 1], rev(basis_col), type = "open")
        X_design[, j] <- conv_full[seq_len(n_time)]
      }
      
      # QR decomposition for this cluster
      qr_decomp <- qr(X_design)
      Q <- qr.Q(qr_decomp)
      R <- qr.R(qr_decomp)
      R_inv <- solve(R + lambda_ridge * diag(ncol(R)))
      
      # Refit voxels in this cluster
      Y_cluster <- Y_proj[, idx_cluster, drop = FALSE]
      coeffs_new <- R_inv %*% t(Q) %*% Y_cluster
      
      # Extract parameters
      beta0_new <- coeffs_new[1, ]
      beta0_safe <- ifelse(abs(beta0_new) < epsilon_beta, epsilon_beta, beta0_new)
      delta_theta <- coeffs_new[2:(n_params+1), , drop = FALSE] / 
                     matrix(rep(beta0_safe, each = n_params), nrow = n_params)
      
      theta_hat_new <- matrix(theta_cluster, nrow = length(idx_cluster), 
                              ncol = n_params, byrow = TRUE) + t(delta_theta)
      
      # Apply bounds
      theta_hat_new <- pmax(theta_bounds$lower, 
                            pmin(theta_hat_new, theta_bounds$upper))
      
      # Calculate new R² for cluster voxels
      fitted_new <- X_design %*% coeffs_new
      residuals_new <- Y_cluster - fitted_new
      SS_res_new <- colSums(residuals_new^2)
      
      Y_means <- colMeans(Y_cluster)
      SS_tot <- colSums((Y_cluster - matrix(Y_means, nrow = n_time, 
                                             ncol = length(idx_cluster), byrow = TRUE))^2)
      r2_new <- 1 - SS_res_new / SS_tot
      
      # Update only if improved
      improved <- r2_new > r2_voxel[idx_cluster]
      idx_update <- idx_cluster[improved]
      
      if (length(idx_update) > 0) {
        theta_hat_voxel[idx_update, ] <- theta_hat_new[improved, , drop = FALSE]
        coeffs_voxel[, idx_update] <- coeffs_new[, improved, drop = FALSE]
        r2_voxel[idx_update] <- r2_new[improved]
      }
      
      cluster_improvements[k] <- mean(r2_new[improved] - r2_orig[idx_cluster[improved]])
      
      if (verbose) {
        cat("    Updated", sum(improved), "of", length(idx_cluster), "voxels\n")
        cat("    Mean R² improvement:", round(cluster_improvements[k], 4), "\n")
      }
    }
    
    # Store clustering information
    cluster_info$n_passes <- pass_km
    cluster_info$cluster_assignments <- cluster_assignment_full
    cluster_info$cluster_centers[[pass_km]] <- cluster_centers
    cluster_info$cluster_r2_improvements[[pass_km]] <- cluster_improvements
    
    # Check for overall improvement
    overall_improvement <- mean(r2_voxel - r2_orig)
    if (verbose) {
      cat("  Overall R² improvement:", round(overall_improvement, 4), "\n")
    }
    
    # Early stopping if minimal improvement
    if (overall_improvement < 0.001) {
      if (verbose) cat("  Minimal improvement; stopping K-means re-centering\n")
      break
    }
  }
  
  # Return updated results
  list(
    theta_hat = theta_hat_voxel,
    coeffs = coeffs_voxel,
    r_squared = r2_voxel,
    cluster_info = cluster_info,
    improvement_summary = list(
      mean_r2_improvement = mean(r2_voxel - r2_orig),
      n_improved = sum(r2_voxel > r2_orig),
      prop_improved = mean(r2_voxel > r2_orig)
    )
  )
}