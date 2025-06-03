#' Internal parametric HRF fitting engine
#'
#' Core implementation of the Taylor approximation method for parametric HRF estimation.
#' This is the main workhorse function that performs voxel-wise parameter estimation.
#'
#' @param Y_proj Numeric matrix of projected BOLD data (timepoints x voxels)
#' @param S_target_proj Numeric matrix of projected stimulus design (timepoints x regressors)
#' @param scan_times Numeric vector of scan acquisition times
#' @param hrf_eval_times Numeric vector of time points for HRF evaluation
#' @param hrf_interface List with HRF function interface
#' @param theta_seed Numeric vector of starting parameters
#' @param theta_bounds List with elements `lower` and `upper`
#' @param lambda_ridge Numeric ridge penalty (default: 0.01)
#' @param epsilon_beta Small value to avoid division by zero when beta0 is extremely small (default: 1e-6)
#' @param verbose Logical whether to print progress (default: FALSE)
#'
#' @return List with elements:
#'   - `theta_hat`: Matrix of parameter estimates (voxels x parameters)
#'   - `beta0`: Numeric vector of amplitudes
#'   - `r_squared`: Numeric vector of R-squared values
#'   - `residuals`: Matrix of residuals (timepoints x voxels)
#'   - `coeffs`: Matrix of linear coefficients
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
  epsilon_beta = 1e-6,
  verbose = FALSE
) {
  # Validate inputs before computing dimensions
  .validate_input(Y_proj, "Y_proj", type = "matrix")
  .validate_input(S_target_proj, "S_target_proj", type = "matrix")
  if (nrow(S_target_proj) != nrow(Y_proj)) {
    stop(
      sprintf(
        "S_target_proj must have %d rows to match Y_proj", nrow(Y_proj)
      ),
      call. = FALSE
    )
  }

  n_time   <- nrow(Y_proj)
  n_vox    <- ncol(Y_proj)
  n_params <- length(theta_seed)

  if (nrow(S_target_proj) != n_time) {
    stop(
      sprintf(
        "S_target_proj has %d rows but expected %d to match Y_proj",
        nrow(S_target_proj), n_time
      ),
      call. = FALSE
    )
  }

  if (verbose) {
    cat("Parametric engine: ", n_vox, " voxels, ", n_params,
        " parameters\n", sep = "")
  }

  # 1. Taylor basis at expansion point
  X_basis <- hrf_interface$taylor_basis(theta_seed, hrf_eval_times)
  if (!is.matrix(X_basis)) {
    X_basis <- matrix(X_basis, ncol = n_params + 1)
  }

  # 2. Convolve stimulus with basis functions
  X_design <- .batch_convolution(S_target_proj, X_basis, n_time)
  
  # Ensure X_design is a matrix
  if (!is.matrix(X_design)) {
    stop(sprintf("X_design is not a matrix after batch convolution. Class: %s, is.null: %s", 
                 paste(class(X_design), collapse=","), is.null(X_design)))
  }

  # 3. Linear solution with ridge regularisation
  coeffs <- .ridge_linear_solve(X_design, Y_proj, lambda_ridge)

  # 4. Extract amplitude and parameter updates
  beta0 <- coeffs[1, ]
  # Avoid division by zero when beta0 is extremely small
  beta0_safe <- ifelse(
    abs(beta0) < epsilon_beta,
    ifelse(beta0 < 0, -epsilon_beta, epsilon_beta),
    beta0
  )
  delta_theta <- coeffs[2:(n_params + 1), , drop = FALSE] /
    matrix(rep(beta0_safe, each = n_params), nrow = n_params)

  theta_hat <- matrix(theta_seed, nrow = n_vox, ncol = n_params, byrow = TRUE) +
    t(delta_theta)

  # 5. Apply bounds if provided
  if (!is.null(theta_bounds)) {
    lower <- matrix(theta_bounds$lower,
                    nrow = n_vox, ncol = n_params, byrow = TRUE)
    upper <- matrix(theta_bounds$upper,
                    nrow = n_vox, ncol = n_params, byrow = TRUE)
    theta_hat <- pmin(upper, pmax(lower, theta_hat))
  }

  if (!is.null(hrf_interface$parameter_names)) {
    colnames(theta_hat) <- hrf_interface$parameter_names
  }

  # 6. Fit quality metrics
  fitted_values <- X_design %*% coeffs
  residuals <- Y_proj - fitted_values
  y_mean <- matrix(colMeans(Y_proj), nrow = n_time, ncol = n_vox, byrow = TRUE)
  ss_tot <- colSums((Y_proj - y_mean)^2)
  ss_res <- colSums(residuals^2)
  r_squared <- ifelse(ss_tot > 1e-10, 1 - ss_res / ss_tot, 0)
  r_squared <- pmax(0, pmin(1, r_squared))

  list(
    theta_hat = theta_hat,
    beta0 = if(is.matrix(beta0)) beta0 else matrix(beta0, ncol = 1),
    r_squared = as.numeric(r_squared),
    residuals = residuals,
    coeffs = coeffs
  )
}