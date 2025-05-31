#' Internal parametric HRF fitting engine (single pass)
#'
#' Implements the core single-pass Taylor approximation used to estimate
#' parametric HRF parameters for a set of voxels. This function is not exported
#' and is intended for internal use.
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
#'
#' @return List with elements `theta_hat` (matrix) and `beta0` (numeric vector)
#' @keywords internal
.parametric_engine <- function(
  Y_proj,
  S_target_proj,
  scan_times,
  hrf_eval_times,
  hrf_interface,
  theta_seed,
  theta_bounds,
  lambda_ridge = 0.01,
  epsilon_beta = 1e-6
) {
  # 1. Taylor basis
  X_taylor <- hrf_interface$taylor_basis(theta_seed, hrf_eval_times)
  if (!is.matrix(X_taylor)) {
    X_taylor <- matrix(X_taylor, ncol = 4)
  }

  n_time <- nrow(S_target_proj)
  n_vox <- ncol(Y_proj)

  # 2. Design matrix via convolution
  X_design <- matrix(0, nrow = n_time, ncol = ncol(X_taylor))
  for (j in seq_len(ncol(X_taylor))) {
    basis_col <- X_taylor[, j]
    conv_full <- stats::convolve(S_target_proj[, 1], rev(basis_col), type = "open")
    X_design[, j] <- conv_full[seq_len(n_time)]
  }

  # 3. Global QR decomposition
  qr_decomp <- qr(X_design)
  Q <- qr.Q(qr_decomp)
  R <- qr.R(qr_decomp)
  R_inv <- solve(R + lambda_ridge * diag(ncol(R)))

  # 4. Voxel-wise estimation
  coeffs <- R_inv %*% t(Q) %*% Y_proj
  beta0 <- coeffs[1, ]
  beta0_safe <- ifelse(abs(beta0) < epsilon_beta, epsilon_beta, beta0)
  delta_theta <- coeffs[2:4, , drop = FALSE] / matrix(rep(beta0_safe, each = 3), nrow = 3)
  theta_hat <- matrix(theta_seed, nrow = n_vox, ncol = length(theta_seed), byrow = TRUE) + t(delta_theta)

  # 5. Parameter bounds
  theta_hat <- pmax(theta_bounds$lower, pmin(theta_hat, theta_bounds$upper))

  list(theta_hat = theta_hat, beta0 = as.numeric(beta0))
}
