#' Single-voxel sanity check
#'
#' Fits the LWU parametric HRF model to a single time series with a
#' single-condition design. This helper exposes the minimal estimation
#' machinery for quick verification that the engine functions correctly.
#'
#' @param y Numeric vector of BOLD data (time points).
#' @param onsets Numeric vector of the same length with 1s at event onsets
#'   and 0 elsewhere.
#' @param parametric_hrf Character string specifying the HRF model to use.
#'   Currently only \code{"lwu"} is supported.
#' @param hrf_eval_times Numeric vector of times at which the HRF is
#'   evaluated. Defaults to \code{seq(0, 30, by = 0.5)}.
#' @param lambda_ridge Numeric ridge penalty used in the linear solve.
#' @param verbose Logical; print progress messages.
#'
#' @return List with elements \code{theta_hat}, \code{beta0} and
#'   \code{r_squared} for the single voxel.
#' @keywords internal
#' @export
single_voxel_sanity_check <- function(
  y,
  onsets,
  parametric_hrf = "lwu",
  hrf_eval_times = seq(0, 30, by = 0.5),
  lambda_ridge = 0.01,
  verbose = FALSE
) {
  stopifnot(is.numeric(y), is.numeric(onsets), length(y) == length(onsets))

  Y <- matrix(as.numeric(y), ncol = 1)
  S <- matrix(as.numeric(onsets), ncol = 1)

  # Get base HRF interface
  base_interface <- .get_hrf_interface(parametric_hrf)
  theta_seed <- base_interface$default_seed()
  theta_bounds <- base_interface$default_bounds()

  # Create properly wrapped HRF interface with bounds
  hrf_interface <- list(
    hrf_function = function(t, params_vector, ...) {
      base_interface$hrf_function(t, params_vector, bounds = theta_bounds, ...)
    },
    taylor_basis = function(params_vector0, t_hrf_eval, ...) {
      base_interface$taylor_basis(params_vector0, t_hrf_eval, bounds = theta_bounds, ...)
    },
    parameter_names = base_interface$parameter_names,
    default_seed = base_interface$default_seed,
    default_bounds = base_interface$default_bounds
  )

  .parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    hrf_eval_times = hrf_eval_times,
    hrf_interface = hrf_interface,
    theta_seed = theta_seed,
    theta_bounds = theta_bounds,
    lambda_ridge = lambda_ridge,
    verbose = verbose,
    baseline_model = "intercept"
  )
}
