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
  verbose = FALSE
) {
  
  # Use simple, dependency-free engine for stability
  .parametric_engine_simple(
    Y_proj = Y_proj,
    S_target_proj = S_target_proj,
    scan_times = scan_times,
    hrf_eval_times = hrf_eval_times,
    hrf_interface = hrf_interface,
    theta_seed = theta_seed,
    theta_bounds = theta_bounds,
    lambda_ridge = lambda_ridge,
    verbose = verbose
  )
}