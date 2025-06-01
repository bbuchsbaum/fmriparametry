#' Simple parametric HRF fitting engine (no dependencies)
#'
#' A completely self-contained implementation of the Taylor approximation
#' method that avoids all external dependencies and source() calls.
#'
#' @param Y_proj Numeric matrix of projected BOLD data (timepoints x voxels)
#' @param S_target_proj Numeric matrix of projected stimulus design (timepoints x regressors)
#' @param scan_times Numeric vector of scan acquisition times
#' @param hrf_eval_times Numeric vector of time points for HRF evaluation
#' @param hrf_interface List with HRF function interface
#' @param theta_seed Numeric vector of starting parameters
#' @param theta_bounds List with elements `lower` and `upper`
#' @param lambda_ridge Numeric ridge penalty (default: 0.01)
#' @param verbose Logical whether to print progress (default: FALSE)
#'
#' @return List with elements:
#'   - `theta_hat`: Matrix of parameter estimates (voxels x parameters)
#'   - `beta0`: Numeric vector of amplitudes
#'   - `r_squared`: Numeric vector of R-squared values
#'   - `residuals`: Matrix of residuals (timepoints x voxels)
#'   - `coeffs`: Matrix of linear coefficients
#' @keywords internal
.parametric_engine_simple <- function(
  Y_proj,
  S_target_proj,
  scan_times,
  hrf_eval_times,
  hrf_interface,
  theta_seed,
  theta_bounds,
  lambda_ridge = 0.01,
  verbose = FALSE
) {
  
  # Basic input validation
  n_time <- nrow(Y_proj)
  n_vox <- ncol(Y_proj)
  n_params <- length(theta_seed)
  
  if (verbose) cat("Simple engine: processing", n_vox, "voxels with", n_params, "parameters\n")
  
  # Get Taylor basis at seed parameters
  X_taylor <- hrf_interface$taylor_basis(theta_seed, hrf_eval_times)
  if (!is.matrix(X_taylor)) {
    X_taylor <- matrix(X_taylor, ncol = n_params + 1)
  }
  
  # Convolve with stimulus to create design matrix
  X_design <- matrix(0, nrow = n_time, ncol = ncol(X_taylor))
  for (j in seq_len(ncol(X_taylor))) {
    basis_col <- X_taylor[, j]
    for (k in seq_len(ncol(S_target_proj))) {
      stimulus_col <- S_target_proj[, k]
      convolved <- stats::convolve(basis_col, rev(stimulus_col), type = "open")
      X_design[, j] <- X_design[, j] + convolved[seq_len(n_time)]
    }
  }
  
  # QR decomposition with ridge regularization
  qr_decomp <- qr(X_design)
  Q <- qr.Q(qr_decomp)
  R <- qr.R(qr_decomp)
  
  # Add ridge penalty to diagonal
  R_ridge <- R + lambda_ridge * diag(ncol(R))
  
  # Solve for coefficients
  coeffs <- solve(R_ridge) %*% t(Q) %*% Y_proj
  
  # Extract amplitudes (first coefficient) and parameter changes
  beta0 <- coeffs[1, ]
  
  # Prevent division by very small amplitudes
  beta0_safe <- ifelse(abs(beta0) < 1e-6, sign(beta0) * 1e-6, beta0)
  
  # Parameter updates (scaled by amplitude)
  delta_theta <- coeffs[2:(n_params+1), , drop = FALSE] / 
                 matrix(rep(beta0_safe, each = n_params), nrow = n_params)
  
  # Updated parameters
  theta_hat <- matrix(theta_seed, nrow = n_vox, ncol = n_params, byrow = TRUE) + 
               t(delta_theta)
  
  # Apply bounds
  if (!is.null(theta_bounds)) {
    for (j in seq_len(n_params)) {
      theta_hat[, j] <- pmax(theta_bounds$lower[j], 
                            pmin(theta_hat[, j], theta_bounds$upper[j]))
    }
  }
  
  # Set parameter names
  if (!is.null(hrf_interface$parameter_names)) {
    colnames(theta_hat) <- hrf_interface$parameter_names
  }
  
  # Compute fitted values and R-squared
  fitted_values <- X_design %*% coeffs
  residuals <- Y_proj - fitted_values
  
  # Total sum of squares (for R-squared)
  Y_mean <- matrix(colMeans(Y_proj), nrow = n_time, ncol = n_vox, byrow = TRUE)
  SS_tot <- colSums((Y_proj - Y_mean)^2)
  SS_res <- colSums(residuals^2)
  
  # Avoid division by zero
  r_squared <- ifelse(SS_tot > 1e-10, 1 - SS_res / SS_tot, 0)
  r_squared <- pmax(0, pmin(1, r_squared))  # Clamp to [0,1]
  
  if (verbose) {
    cat("  R² range:", round(range(r_squared, na.rm = TRUE), 3), "\n")
  }
  
  # Return results
  list(
    theta_hat = theta_hat,
    beta0 = as.numeric(beta0),
    r_squared = as.numeric(r_squared),
    residuals = residuals,
    coeffs = coeffs
  )
}