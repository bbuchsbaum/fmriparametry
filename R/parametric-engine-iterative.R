#' Internal parametric HRF fitting engine with iterative refinement
#'
#' Implements the core Taylor approximation with optional iterative global 
#' re-centering to improve parameter estimates. This enhanced version calculates
#' R-squared values and can perform multiple passes to refine the expansion point.
#'
#' @param Y_proj Numeric matrix of projected BOLD data (timepoints x voxels)
#' @param S_target_proj Numeric matrix of projected stimulus design (timepoints x regressors)
#' @param scan_times Numeric vector of scan acquisition times
#' @param hrf_eval_times Numeric vector of time points for HRF evaluation
#' @param hrf_interface List with at least element `taylor_basis` producing the
#'   Taylor basis matrix
#' @param theta_seed Numeric vector of starting parameters
#' @param theta_bounds List with elements `lower` and `upper`
#' @param lambda_ridge Numeric ridge penalty applied to the QR R matrix
#' @param epsilon_beta Small numeric to avoid division by zero when beta is near zero
#' @param recenter_global_passes Integer number of global re-centering iterations
#' @param recenter_epsilon Numeric convergence tolerance for re-centering
#' @param r2_threshold Numeric R-squared threshold for selecting good voxels
#' @param compute_residuals Logical whether to compute and return residuals
#' @param compute_se Logical whether to compute standard errors via Delta method
#' @param lambda_ridge_jacobian Numeric ridge penalty for SE calculation
#' @param recenter_kmeans_passes Integer number of K-means re-centering passes
#' @param kmeans_k Integer number of clusters for K-means
#' @param r2_threshold_kmeans Numeric R² threshold for K-means clustering
#' @param verbose Logical whether to print progress messages
#'
#' @return List with elements:
#'   - `theta_hat`: Matrix of parameter estimates (voxels x parameters)
#'   - `beta0`: Numeric vector of amplitudes
#'   - `r_squared`: Numeric vector of R-squared values
#'   - `residuals`: Optional matrix of residuals (timepoints x voxels)
#'   - `se_theta_hat`: Optional matrix of standard errors (voxels x parameters)
#'   - `convergence_info`: List with convergence details
#'   - `coeffs`: Matrix of linear coefficients (for SE calculation)
#' @keywords internal
.parametric_engine_iterative <- function(
  Y_proj,
  S_target_proj,
  scan_times,
  hrf_eval_times,
  hrf_interface,
  theta_seed,
  theta_bounds,
  lambda_ridge = 0.01,
  epsilon_beta = 1e-6,
  recenter_global_passes = 3,
  recenter_epsilon = 0.01,
  r2_threshold = 0.1,
  compute_residuals = TRUE,
  compute_se = FALSE,
  lambda_ridge_jacobian = 0.01,
  recenter_kmeans_passes = 0,
  kmeans_k = 5,
  r2_threshold_kmeans = 0.7,
  verbose = FALSE
) {
  n_time <- nrow(S_target_proj)
  n_vox <- ncol(Y_proj)
  n_params <- length(theta_seed)
  
  # Initialize tracking variables
  theta_current_global <- theta_seed
  r2_voxel <- rep(-Inf, n_vox)
  theta_trajectory <- list(theta_seed)
  
  # Storage for best estimates per voxel
  theta_hat_voxel <- matrix(theta_seed, nrow = n_vox, ncol = n_params, byrow = TRUE)
  beta0_voxel <- rep(0, n_vox)
  coeffs_voxel <- matrix(0, nrow = ncol(hrf_interface$taylor_basis(theta_seed, hrf_eval_times)), 
                         ncol = n_vox)
  
  # Calculate total sum of squares once
  Y_means <- colMeans(Y_proj)
  SS_tot <- colSums((Y_proj - matrix(Y_means, nrow = n_time, ncol = n_vox, byrow = TRUE))^2)
  
  if (verbose) cat("Starting iterative global re-centering...\n")
  
  # Main iteration loop
  for (pass_global in seq_len(recenter_global_passes)) {
    if (verbose) cat("Global pass", pass_global, "with theta =", round(theta_current_global, 3), "\n")
    
    # Perform single Taylor pass
    X_taylor <- hrf_interface$taylor_basis(theta_current_global, hrf_eval_times)
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
    
    # Global QR decomposition
    # Use fast C++ ridge solver
    coeffs_pass <- .ridge_linear_solve(X_design, Y_proj, lambda_ridge)
    beta0_pass <- coeffs_pass[1, ]
    beta0_safe <- ifelse(abs(beta0_pass) < epsilon_beta, epsilon_beta, beta0_pass)
    delta_theta <- coeffs_pass[2:(n_params+1), , drop = FALSE] / 
                   matrix(rep(beta0_safe, each = n_params), nrow = n_params)
    theta_hat_pass <- matrix(theta_current_global, nrow = n_vox, ncol = n_params, byrow = TRUE) + 
                      t(delta_theta)
    
    # Apply bounds
    theta_hat_pass <- pmax(theta_bounds$lower, pmin(theta_hat_pass, theta_bounds$upper))
    
    # Calculate R-squared for current pass
    fitted_values <- X_design %*% coeffs_pass
    residuals_pass <- Y_proj - fitted_values
    SS_res <- colSums(residuals_pass^2)
    r2_pass <- 1 - SS_res / SS_tot
    
    # Update best voxel estimates if R² improved
    improved <- r2_pass > r2_voxel
    if (any(improved)) {
      theta_hat_voxel[improved, ] <- theta_hat_pass[improved, , drop = FALSE]
      beta0_voxel[improved] <- beta0_pass[improved]
      r2_voxel[improved] <- r2_pass[improved]
      coeffs_voxel[, improved] <- coeffs_pass[, improved, drop = FALSE]
    }
    
    if (verbose) {
      cat("  R² range:", round(range(r2_pass, na.rm = TRUE), 3), 
          "Mean:", round(mean(r2_pass, na.rm = TRUE), 3), "\n")
      cat("  Improved voxels:", sum(improved), "/", n_vox, "\n")
    }
    
    # Re-center global θ₀ (if not final pass)
    if (pass_global < recenter_global_passes) {
      # Select good voxels
      idx_good <- which(r2_voxel >= r2_threshold)
      
      if (length(idx_good) >= 10) {
        # Compute robust median of good voxels
        theta_new <- apply(theta_hat_voxel[idx_good, , drop = FALSE], 2, median, na.rm = TRUE)
        
        # Apply bounds
        theta_new <- pmax(theta_bounds$lower, pmin(theta_new, theta_bounds$upper))
        
        # Check convergence
        if (max(abs(theta_new - theta_current_global)) < recenter_epsilon) {
          if (verbose) cat("Converged after", pass_global, "iterations\n")
          break
        }
        
        # Update global theta
        theta_current_global <- theta_new
        theta_trajectory[[pass_global + 1]] <- theta_new
      } else {
        warning("Too few good voxels (", length(idx_good), ") for re-centering; stopping early")
        break
      }
    }
  }
  
  # Compute final residuals if requested
  residuals_matrix <- NULL
  if (compute_residuals) {
    # Reconstruct final design matrix using final global theta
    X_taylor_final <- hrf_interface$taylor_basis(theta_current_global, hrf_eval_times)
    if (!is.matrix(X_taylor_final)) {
      X_taylor_final <- matrix(X_taylor_final, ncol = n_params + 1)
    }
    
    X_design_final <- matrix(0, nrow = n_time, ncol = ncol(X_taylor_final))
    for (j in seq_len(ncol(X_taylor_final))) {
      basis_col <- X_taylor_final[, j]
      conv_full <- stats::convolve(S_target_proj[, 1], rev(basis_col), type = "open")
      X_design_final[, j] <- conv_full[seq_len(n_time)]
    }
    
    # Calculate residuals for each voxel
    fitted_final <- X_design_final %*% coeffs_voxel
    residuals_matrix <- Y_proj - fitted_final
  }
  
  # Prepare convergence information
  convergence_info <- list(
    trajectory = theta_trajectory,
    n_iterations = length(theta_trajectory),
    final_global_theta = theta_current_global,
    converged = length(theta_trajectory) < recenter_global_passes
  )
  
  # Apply K-means re-centering if requested (disabled for stability)
  kmeans_info <- NULL
  if (recenter_kmeans_passes > 0 && kmeans_k > 1) {
    if (verbose) cat("\nK-means re-centering disabled in basic engine for stability\n")
    
    # Mark as not applied
    convergence_info$kmeans_applied <- FALSE
    convergence_info$kmeans_improvement <- "skipped_for_stability"
  }
  
  # Compute standard errors if requested
  se_theta_hat <- NULL
  if (compute_se) {
    if (verbose) cat("Computing standard errors via Delta method...\n")
    
    # Initialize SE matrix
    se_theta_hat <- matrix(NA, nrow = n_vox, ncol = n_params)
    
    # Get final design matrix (already computed if residuals were calculated)
    if (!compute_residuals) {
      X_taylor_final <- hrf_interface$taylor_basis(theta_current_global, hrf_eval_times)
      if (!is.matrix(X_taylor_final)) {
        X_taylor_final <- matrix(X_taylor_final, ncol = n_params + 1)
      }
      
      X_design_final <- matrix(0, nrow = n_time, ncol = ncol(X_taylor_final))
      for (j in seq_len(ncol(X_taylor_final))) {
        basis_col <- X_taylor_final[, j]
        conv_full <- stats::convolve(S_target_proj[, 1], rev(basis_col), type = "open")
        X_design_final[, j] <- conv_full[seq_len(n_time)]
      }
    }
    
    # Compute X'X once
    XtX <- crossprod(X_design_final)
    
    # For each voxel
    for (v in seq_len(n_vox)) {
      # Skip if poor fit
      if (r2_voxel[v] < 0) next
      
      # Get coefficients for this voxel
      coeffs_v <- coeffs_voxel[, v]
      
      # Calculate residuals if not already done
      if (is.null(residuals_matrix)) {
        fitted_v <- X_design_final %*% coeffs_v
        resid_v <- Y_proj[, v] - fitted_v
      } else {
        resid_v <- residuals_matrix[, v]
      }
      
      # Error variance
      sigma2_v <- sum(resid_v^2) / (n_time - (n_params + 1))
      
      # Covariance matrix of coefficients
      Sigma_coeffs_v <- sigma2_v * solve(XtX + lambda_ridge_jacobian * diag(ncol(XtX)))
      
      # Delta method for SE of theta
      # theta = theta_expansion + delta_theta
      # delta_theta_k = coeffs_{k+1} / coeffs_1
      
      # Jacobian of g(coeffs) = (coeffs_2/coeffs_1, ..., coeffs_{P+1}/coeffs_1)
      beta0_v <- coeffs_v[1]
      if (abs(beta0_v) < epsilon_beta) {
        # Skip SE calculation for near-zero amplitude
        next
      }
      
      # Construct Jacobian matrix
      J_g <- matrix(0, nrow = n_params, ncol = n_params + 1)
      
      # Partial derivatives
      for (k in seq_len(n_params)) {
        # d(delta_theta_k)/d(beta0) = -coeffs_{k+1} / beta0^2
        J_g[k, 1] <- -coeffs_v[k + 1] / (beta0_v^2)
        
        # d(delta_theta_k)/d(coeffs_{k+1}) = 1 / beta0
        J_g[k, k + 1] <- 1 / beta0_v
      }
      
      # Covariance of delta_theta via Delta method
      Sigma_delta_theta_v <- J_g %*% Sigma_coeffs_v %*% t(J_g)
      
      # Extract standard errors (square root of diagonal)
      se_theta_hat[v, ] <- sqrt(pmax(0, diag(Sigma_delta_theta_v)))
    }
  }
  
  # Return enhanced output
  list(
    theta_hat = theta_hat_voxel,
    beta0 = as.numeric(beta0_voxel),
    r_squared = r2_voxel,
    residuals = residuals_matrix,
    se_theta_hat = se_theta_hat,
    convergence_info = convergence_info,
    kmeans_info = kmeans_info,
    coeffs = coeffs_voxel,
    theta_expansion = theta_current_global  # Final expansion point
  )
}