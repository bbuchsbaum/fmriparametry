#' Internal parametric HRF fitting engine
#'
#' Core implementation of the Taylor approximation method for parametric HRF estimation.
#' This is the main workhorse function that performs voxel-wise parameter estimation.
#'
#' @param Y_proj Numeric matrix of projected BOLD data (timepoints x voxels)
#' @param S_target_proj Numeric matrix of projected stimulus design (timepoints x regressors)
#' @param hrf_eval_times Numeric vector of time points for HRF evaluation
#' @param hrf_interface List with HRF function interface
#' @param theta_seed Numeric vector of starting parameters
#' @param theta_bounds List with elements `lower` and `upper`
#' @param lambda_ridge Numeric ridge penalty (default: 0.01)
#' @param epsilon_beta Small value to avoid division by zero when beta0 is extremely small (default: 1e-6)
#' @param verbose Logical whether to print progress (default: FALSE)
#' @param baseline_model Character string specifying baseline model (default: "intercept")
#'
#' @note
#' This engine assumes evenly spaced scan times (constant TR).
#' If acquisition times vary, adjust the stimulus design or
#' interpolate the data before calling this function.
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
  hrf_eval_times,
  hrf_interface,
  theta_seed,
  theta_bounds,
  lambda_ridge = 0.01,
  epsilon_beta = 1e-6,
  verbose = FALSE,
  baseline_model = "intercept"
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

  # 2b. Add intercept/baseline column if requested
  has_intercept <- FALSE
  if (!is.null(baseline_model) && baseline_model == "intercept") {
    intercept_col <- matrix(1.0, nrow = n_time, ncol = 1)
    X_design <- cbind(intercept_col, X_design)
    has_intercept <- TRUE
  }

  # 3. Linear solution with ridge regularisation
  # Debug: Check types before calling C++
  if (getOption("fmriparametric.debug", FALSE)) {
    cat("\nDEBUG before ridge_linear_solve:\n")
    cat("X_design class:", class(X_design), "\n")
    cat("X_design typeof:", typeof(X_design), "\n")
    cat("X_design dim:", dim(X_design), "\n")
    cat("Y_proj class:", class(Y_proj), "\n")
    cat("Y_proj typeof:", typeof(Y_proj), "\n")
    cat("Y_proj dim:", dim(Y_proj), "\n")
  }
  
  # Ensure matrices are double-precision before C++ call
  storage.mode(X_design) <- "double"
  storage.mode(Y_proj) <- "double"
  
  coeffs <- .ridge_linear_solve(X_design, Y_proj, lambda_ridge)
  
  # DEBUG: Check coefficients
  if (getOption("fmriparametric.debug", FALSE)) {
    cat("\nDEBUG: Ridge regression results:\n")
    cat("Coeffs dim:", dim(coeffs), "\n")
    cat("Coeffs range:", range(coeffs), "\n")
    if (has_intercept) {
      cat("Intercept (row 1):", coeffs[1,], "\n")
      cat("Beta0 (row 2):", coeffs[2,], "\n")
    } else {
      cat("Beta0 (row 1):", coeffs[1,], "\n")
    }
  }

  # 4. Extract amplitude and parameter updates
  if (has_intercept) {
    # Skip the first row (intercept) when extracting HRF coefficients
    beta0 <- coeffs[2, ]  # HRF amplitude is now in row 2
    coeff_start_idx <- 3   # Derivative coefficients start at row 3
  } else {
    beta0 <- coeffs[1, ]   # HRF amplitude is in row 1
    coeff_start_idx <- 2   # Derivative coefficients start at row 2
  }
  
  # Avoid division by zero when beta0 is extremely small
  beta0_safe <- ifelse(
    abs(beta0) < epsilon_beta,
    ifelse(beta0 < 0, -epsilon_beta, epsilon_beta),
    beta0
  )
  delta_theta <- coeffs[coeff_start_idx:(coeff_start_idx + n_params - 1), , drop = FALSE] /
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
  
  # Calculate comprehensive fit metrics
  fit_metrics <- .calculate_fit_metrics(
    y_true = Y_proj,
    y_pred = fitted_values,
    n_predictors = n_params,
    has_intercept = has_intercept
  )
  
  r_squared <- fit_metrics$r_squared

  # DEBUG: Add detailed diagnostics
  if (any(r_squared == 0) || all(r_squared < 0.01)) {
    cat("\n=== PARAMETRIC ENGINE DEBUG (Low R² detected) ===\n")
    cat("Number of voxels with R² = 0:", sum(r_squared == 0), "\n")
    cat("Mean R²:", mean(r_squared), "\n")
    cat("Max R²:", max(r_squared), "\n")
    
    # Check a sample voxel
    voxel_vars <- apply(Y_proj, 2, var)
    v_idx <- which.max(voxel_vars)[1]  # Pick highest variance voxel
    if (!is.na(v_idx)) {
      cat("\nDiagnosing highest variance voxel", v_idx, ":\n")
      cat("Y data range:", range(Y_proj[, v_idx]), "\n")
      cat("Y variance:", var(Y_proj[, v_idx]), "\n")
      cat("Y mean:", mean(Y_proj[, v_idx]), "\n")
      cat("X_design dim:", dim(X_design), "\n")
      cat("Has intercept:", has_intercept, "\n")
      if (has_intercept) {
        cat("X_design[,1] (intercept) range:", range(X_design[, 1]), "\n")
        cat("X_design[,2] (HRF) range:", range(X_design[, 2]), "\n")
      } else {
        cat("X_design[,1] (HRF) range:", range(X_design[, 1]), "\n")
      }
      cat("Beta values:", coeffs[, v_idx], "\n")
      cat("Fitted values range:", range(fitted_values[, v_idx]), "\n")
      cat("Fitted values mean:", mean(fitted_values[, v_idx]), "\n")
      cat("Fitted values variance:", var(fitted_values[, v_idx]), "\n")
      cat("SS_res:", fit_metrics$rss[v_idx], "\n")
      cat("SS_tot:", fit_metrics$tss[v_idx], "\n")
      cat("R²:", r_squared[v_idx], "\n")
      
      # Check if design matrix is all zeros
      if (all(abs(X_design) < 1e-10)) {
        cat("\nWARNING: Design matrix X_design is all zeros!\n")
      }
      
      # Check convolution result
      cat("\nConvolution check:\n")
      cat("S_target_proj dim:", dim(S_target_proj), "\n")
      cat("S_target_proj[,1] range:", range(S_target_proj[,1]), "\n")
      cat("X_basis dim:", dim(X_basis), "\n")
      cat("X_basis[,1] (HRF) range:", range(X_basis[,1]), "\n")
    }
    cat("=================================\n\n")
  }

  list(
    theta_hat = theta_hat,
    beta0 = if(is.matrix(beta0)) beta0 else matrix(beta0, ncol = 1),
    r_squared = as.numeric(r_squared),
    residuals = residuals,
    coeffs = coeffs
  )
}