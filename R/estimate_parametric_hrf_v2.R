#' Estimate parametric HRF parameters (Sprint 2 version)
#'
#' Main function for estimating voxel-wise parameters of parametric Hemodynamic 
#' Response Function (HRF) models from fMRI data. This enhanced version implements
#' iterative global re-centering for improved parameter estimates.
#'
#' @param fmri_data An fMRI dataset object from the \code{fmrireg} package (e.g., 
#'   \code{fmri_dataset} or \code{matrix_dataset}), or a numeric matrix with 
#'   time points in rows and voxels in columns.
#' @param event_model An event model object from \code{fmrireg::event_model()} 
#'   describing stimulus timing, or a numeric matrix with time points in rows.
#' @param parametric_hrf Character string specifying the parametric HRF model. 
#'   Currently only \code{"lwu"} (Lag-Width-Undershoot) is supported.
#' @param theta_seed Numeric vector of initial parameter values, \code{"data_median"}
#'   for data-driven initialization, or \code{NULL} to use model defaults.
#' @param theta_bounds List with elements \code{lower} and \code{upper}, each 
#'   being numeric vectors of parameter bounds, or \code{NULL} for model defaults.
#' @param confound_formula Formula specifying confound regressors (not implemented).
#' @param baseline_model Character string specifying baseline model. Default is 
#'   \code{"intercept"}.
#' @param hrf_eval_times Numeric vector of time points (in seconds) at which to 
#'   evaluate the HRF, or \code{NULL} to use automatic spacing.
#' @param hrf_span Positive numeric giving the duration (in seconds) over which 
#'   to evaluate the HRF. Default is 30 seconds.
#' @param lambda_ridge Ridge regularization penalty for numerical stability. 
#'   Default is 0.01.
#' @param recenter_global_passes Integer number of global re-centering iterations.
#'   Default is 3. Set to 0 for single-pass (Sprint 1 behavior).
#' @param recenter_epsilon Numeric convergence tolerance for re-centering.
#'   Default is 0.01.
#' @param r2_threshold Numeric R-squared threshold for selecting good voxels
#'   during re-centering. Default is 0.1.
#' @param compute_se Logical whether to compute standard errors via Delta method.
#'   Default is TRUE.
#' @param mask Optional logical/numeric mask for spatial subsetting of voxels.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#'
#' @return An object of class \code{parametric_hrf_fit} containing:
#'   \item{estimated_parameters}{Matrix of parameter estimates (voxels × parameters)}
#'   \item{amplitudes}{Numeric vector of amplitude estimates for each voxel}
#'   \item{parameter_names}{Character vector of parameter names}
#'   \item{hrf_model}{Character string identifying the HRF model used}
#'   \item{r_squared}{Numeric vector of R-squared values for each voxel}
#'   \item{residuals}{Matrix of residuals (timepoints × voxels) if computed}
#'   \item{convergence_info}{List with convergence details}
#'   \item{metadata}{List containing fitting details and settings}
#'
#' @details
#' The function implements an iterative Taylor approximation method with global
#' re-centering for improved parameter estimation:
#' \enumerate{
#'   \item Optionally performs data-driven initialization
#'   \item Iteratively refines the global expansion point based on well-fit voxels
#'   \item Tracks R-squared values to identify improved estimates
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
#' # Basic example with iterative refinement
#' fit <- estimate_parametric_hrf_v2(
#'   fmri_data = Y,
#'   event_model = events,
#'   parametric_hrf = "lwu",
#'   recenter_global_passes = 3,
#'   verbose = TRUE
#' )
#' 
#' # Data-driven initialization
#' fit <- estimate_parametric_hrf_v2(
#'   fmri_data = fmri_data,
#'   event_model = event_model,
#'   theta_seed = "data_median",
#'   verbose = TRUE
#' )
#' }
#' 
#' @export
estimate_parametric_hrf_v2 <- function(
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
  recenter_global_passes = 3,
  recenter_epsilon = 0.01,
  r2_threshold = 0.1,
  compute_se = TRUE,
  mask = NULL,
  verbose = FALSE
) {
  # Input validation (same as v1)
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
  
  # Handle theta_bounds
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
  
  # Prepare inputs
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
  
  # Handle theta_seed including "data_median" option
  if (is.null(theta_seed)) {
    theta_seed_final <- hrf_interface$default_seed
  } else if (identical(theta_seed, "data_median")) {
    if (verbose) message("Computing data-driven seed...")
    
    # Preliminary pass with default seed
    prelim_res <- .parametric_engine_iterative(
      Y_proj = prepared$Y_proj,
      S_target_proj = prepared$S_target_proj,
      scan_times = prepared$scan_times,
      hrf_eval_times = prepared$hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_seed = hrf_interface$default_seed,
      theta_bounds = theta_bounds,
      lambda_ridge = lambda_ridge,
      recenter_global_passes = 1,  # Single pass for preliminary
      compute_residuals = FALSE,
      verbose = FALSE
    )
    
    # Select good voxels based on R²
    r2_threshold_prelim <- quantile(prelim_res$r_squared, 0.75, na.rm = TRUE)
    idx_good <- which(prelim_res$r_squared >= r2_threshold_prelim)
    
    if (length(idx_good) >= 10) {
      theta_seed_final <- apply(
        prelim_res$theta_hat[idx_good, , drop = FALSE],
        2,
        median,
        na.rm = TRUE
      )
      if (verbose) {
        message("Data-driven seed from ", length(idx_good), " good voxels: ", 
                paste(round(theta_seed_final, 3), collapse = ", "))
      }
    } else {
      warning("Too few good voxels for data-driven seed; using default")
      theta_seed_final <- hrf_interface$default_seed
    }
  } else {
    # Numeric seed provided
    assertthat::assert_that(
      is.numeric(theta_seed),
      length(theta_seed) == length(hrf_interface$parameter_names),
      msg = "theta_seed must be numeric and match number of parameters"
    )
    theta_seed_final <- theta_seed
  }
  
  # Validate seed is within bounds
  assertthat::assert_that(
    all(theta_seed_final >= theta_bounds$lower),
    all(theta_seed_final <= theta_bounds$upper),
    msg = "theta_seed must fall within theta_bounds"
  )
  
  # Additional parameter validation
  assertthat::assert_that(is.numeric(recenter_global_passes), 
                          length(recenter_global_passes) == 1,
                          recenter_global_passes >= 0)
  assertthat::assert_that(is.numeric(recenter_epsilon), 
                          length(recenter_epsilon) == 1,
                          recenter_epsilon > 0)
  assertthat::assert_that(is.numeric(r2_threshold), 
                          length(r2_threshold) == 1,
                          r2_threshold >= 0,
                          r2_threshold <= 1)
  assertthat::assert_that(is.logical(compute_se), length(compute_se) == 1)
  
  # Run main estimation
  if (verbose) message("Running parametric engine with iterative refinement...")
  fit_res <- .parametric_engine_iterative(
    Y_proj = prepared$Y_proj,
    S_target_proj = prepared$S_target_proj,
    scan_times = prepared$scan_times,
    hrf_eval_times = prepared$hrf_eval_times,
    hrf_interface = hrf_interface,
    theta_seed = theta_seed_final,
    theta_bounds = theta_bounds,
    lambda_ridge = lambda_ridge,
    recenter_global_passes = recenter_global_passes,
    recenter_epsilon = recenter_epsilon,
    r2_threshold = r2_threshold,
    compute_residuals = TRUE,
    compute_se = compute_se,
    verbose = verbose
  )
  
  # Construct enhanced output
  if (verbose) message("Constructing output...")
  new_parametric_hrf_fit(
    estimated_parameters = fit_res$theta_hat,
    amplitudes = fit_res$beta0,
    parameter_names = hrf_interface$parameter_names,
    hrf_model = parametric_hrf,
    r_squared = fit_res$r_squared,
    residuals = fit_res$residuals,
    parameter_ses = fit_res$se_theta_hat,
    convergence_info = fit_res$convergence_info,
    metadata = list(
      call = match.call(),
      n_voxels = ncol(prepared$Y_raw),
      n_timepoints = nrow(prepared$Y_raw),
      theta_seed = if(identical(theta_seed, "data_median")) "data_median" else theta_seed_final,
      theta_bounds = theta_bounds,
      recenter_global_passes = recenter_global_passes,
      coeffs = fit_res$coeffs,  # Store for SE calculation later
      theta_expansion = fit_res$theta_expansion
    )
  )
}