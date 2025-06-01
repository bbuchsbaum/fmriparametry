#' Estimate parametric HRF parameters
#'
#' Main function for estimating voxel-wise parameters of parametric Hemodynamic 
#' Response Function (HRF) models from fMRI data. Currently implements the 
#' Lag-Width-Undershoot (LWU) model using a single-pass Taylor approximation method.
#'
#' @param fmri_data An fMRI dataset object from the \code{fmrireg} package (e.g., 
#'   \code{fmri_dataset} or \code{matrix_dataset}), or a numeric matrix with 
#'   time points in rows and voxels in columns.
#' @param event_model An event model object from \code{fmrireg::event_model()} 
#'   describing stimulus timing, or a numeric matrix with time points in rows.
#' @param parametric_hrf Character string specifying the parametric HRF model. 
#'   Currently only \code{"lwu"} (Lag-Width-Undershoot) is supported.
#' @param theta_seed Numeric vector of initial parameter values, or \code{NULL} 
#'   to use model defaults. For LWU: \code{c(tau, sigma, rho)}.
#' @param theta_bounds List with elements \code{lower} and \code{upper}, each 
#'   being numeric vectors of parameter bounds, or \code{NULL} for model defaults.
#' @param confound_formula Formula specifying confound regressors (not implemented 
#'   in Sprint 1).
#' @param baseline_model Character string specifying baseline model. Default is 
#'   \code{"intercept"}.
#' @param hrf_eval_times Numeric vector of time points (in seconds) at which to 
#'   evaluate the HRF, or \code{NULL} to use automatic spacing.
#' @param hrf_span Positive numeric giving the duration (in seconds) over which 
#'   to evaluate the HRF. Default is 30 seconds.
#' @param lambda_ridge Ridge regularization penalty for numerical stability. 
#'   Default is 0.01.
#' @param mask Optional logical/numeric mask for spatial subsetting of voxels.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#'
#' @return An object of class \code{parametric_hrf_fit} containing:
#'   \item{estimated_parameters}{Matrix of parameter estimates (voxels × parameters)}
#'   \item{amplitudes}{Numeric vector of amplitude estimates for each voxel}
#'   \item{parameter_names}{Character vector of parameter names}
#'   \item{hrf_model}{Character string identifying the HRF model used}
#'   \item{metadata}{List containing fitting details and settings}
#'
#' @details
#' The function implements a fast single-pass Taylor approximation method for 
#' estimating HRF parameters. The algorithm:
#' \enumerate{
#'   \item Constructs a Taylor expansion of the HRF around initial parameter values
#'   \item Solves a linear system to estimate parameter updates
#'   \item Applies parameter bounds to ensure physiological plausibility
#' }
#' 
#' For the LWU model, the three parameters are:
#' \itemize{
#'   \item \code{tau}: Time to peak (lag) in seconds
#'   \item \code{sigma}: Width of the response in seconds
#'   \item \code{rho}: Undershoot amplitude (0-1.5)
#' }
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(fmriparametric)
#' library(fmrireg)
#' 
#' # Basic example with matrix inputs
#' Y <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' events <- matrix(rbinom(100, 1, 0.1), ncol = 1)
#' 
#' fit <- estimate_parametric_hrf(
#'   fmri_data = Y,
#'   event_model = events,
#'   parametric_hrf = "lwu"
#' )
#' 
#' # Example with fmrireg objects
#' fmri_data <- fmrireg::read_fmri("functional.nii.gz")
#' events_df <- data.frame(
#'   onset = c(10, 30, 50, 70),
#'   duration = rep(2, 4),
#'   trial_type = "stimulus"
#' )
#' 
#' event_model <- fmrireg::event_model(
#'   onset ~ trial_type,
#'   data = events_df,
#'   sampling_rate = 0.5  # 2s TR
#' )
#' 
#' # Estimate with custom settings
#' fit <- estimate_parametric_hrf(
#'   fmri_data = fmri_data,
#'   event_model = event_model,
#'   parametric_hrf = "lwu",
#'   theta_seed = c(6, 2.5, 0.35),
#'   theta_bounds = list(
#'     lower = c(3, 1, 0),
#'     upper = c(9, 4, 1)
#'   ),
#'   lambda_ridge = 0.05,
#'   verbose = TRUE
#' )
#' 
#' # View results
#' print(fit)
#' summary(fit)
#' params <- coef(fit)
#' }
#' 
#' @seealso 
#' \code{\link{coef.parametric_hrf_fit}}, \code{\link{summary.parametric_hrf_fit}}, 
#' \code{\link{print.parametric_hrf_fit}}
#' 
#' @references
#' Glover, G. H. (1999). Deconvolution of impulse response in event-related BOLD fMRI. 
#' NeuroImage, 9(4), 416-429.
#' 
# NOT EXPORTED - Deprecated version
estimate_parametric_hrf_v1_DEPRECATED <- function(
  fmri_data,
  event_model,
  parametric_hrf = "lwu",
  theta_seed = NULL,
  theta_bounds = NULL,
  confound_formula = NULL,
  baseline_model = "intercept",
  hrf_eval_times = NULL,
  hrf_span = 30,
  lambda_ridge = 0.01,
  mask = NULL,
  verbose = FALSE
) {
  assertthat::assert_that(!missing(fmri_data), !missing(event_model))
  assertthat::assert_that(is.character(parametric_hrf), length(parametric_hrf) == 1)
  parametric_hrf <- tolower(parametric_hrf)

  if (identical(parametric_hrf, "lwu")) {
    hrf_interface <- list(
      hrf_function = .lwu_hrf_function,
      taylor_basis = .lwu_hrf_taylor_basis_function,
      parameter_names = .lwu_hrf_parameter_names(),
      default_seed = .lwu_hrf_default_seed(),
      default_bounds = .lwu_hrf_default_bounds()
    )
  } else {
    stop("Unsupported parametric_hrf: ", parametric_hrf, call. = FALSE)
  }

  if (is.null(theta_seed)) {
    theta_seed <- hrf_interface$default_seed
  } else {
    assertthat::assert_that(
      is.numeric(theta_seed),
      length(theta_seed) == length(hrf_interface$parameter_names),
      msg = "theta_seed must be numeric and match number of parameters"
    )
  }

  if (is.null(theta_bounds)) {
    theta_bounds <- hrf_interface$default_bounds
  } else {
    assertthat::assert_that(
      is.list(theta_bounds),
      all(c("lower", "upper") %in% names(theta_bounds)),
      is.numeric(theta_bounds$lower),
      is.numeric(theta_bounds$upper),
      length(theta_bounds$lower) == length(hrf_interface$parameter_names),
      length(theta_bounds$upper) == length(hrf_interface$parameter_names),
      msg = "theta_bounds must be a list with numeric 'lower' and 'upper' vectors"
    )
  }

  assertthat::assert_that(all(theta_bounds$lower < theta_bounds$upper),
                          msg = "theta_bounds 'lower' must be less than 'upper'")
  assertthat::assert_that(all(theta_seed >= theta_bounds$lower),
                          all(theta_seed <= theta_bounds$upper),
                          msg = "theta_seed must fall within theta_bounds")

  if (!is.null(hrf_eval_times)) {
    assertthat::assert_that(is.numeric(hrf_eval_times), all(hrf_eval_times >= 0))
  }
  assertthat::assert_that(is.numeric(hrf_span), length(hrf_span) == 1, hrf_span > 0)
  assertthat::assert_that(is.numeric(lambda_ridge), length(lambda_ridge) == 1, lambda_ridge >= 0)
  assertthat::assert_that(is.logical(verbose), length(verbose) == 1)

  if (verbose) message("Preparing inputs...")
  prepared <- .prepare_parametric_inputs(
    fmri_data = fmri_data,
    event_model = event_model,
    confound_formula = confound_formula,
    baseline_model = baseline_model,
    hrf_eval_times = hrf_eval_times,
    hrf_span = hrf_span,
    mask = mask
  )

  if (verbose) message("Running parametric engine...")
  fit_res <- .parametric_engine(
    Y_proj = prepared$Y_proj,
    S_target_proj = prepared$S_target_proj,
    scan_times = prepared$scan_times,
    hrf_eval_times = prepared$hrf_eval_times,
    hrf_interface = hrf_interface,
    theta_seed = theta_seed,
    theta_bounds = theta_bounds,
    lambda_ridge = lambda_ridge
  )

  if (verbose) message("Constructing output...")
  new_parametric_hrf_fit(
    estimated_parameters = fit_res$theta_hat,
    amplitudes = fit_res$beta0,
    parameter_names = hrf_interface$parameter_names,
    hrf_model = parametric_hrf,
    convergence = list(),
    metadata = list(
      call = match.call(),
      n_voxels = ncol(prepared$Y_raw),
      n_timepoints = nrow(prepared$Y_raw),
      theta_seed = theta_seed,
      theta_bounds = theta_bounds
    )
  )
}

