#' Estimate parametric HRF parameters
#'
#' Estimate hemodynamic response function parameters from fMRI time series using
#' a parametric HRF model. Optional refinement steps improve fits for challenging
#' voxels. The core engine performs a single Taylor approximation pass around
#' the provided parameter seed. Additional stages can iteratively update this
#' seed and apply specialized refinement algorithms.
#' 
#' @section Implementation Selection:
#' This function uses the Strangler Fig pattern to safely migrate from the
#' original monolithic implementation to a modular refactored version.
#' The refactored implementation is now the default. Legacy support is available
#' temporarily via `options(fmriparametric.use_legacy = TRUE)`.
#' 
#' @section Package Options:
#' Global iterative refinement (Stage 3) is controlled by the option
#' `fmriparametric.refine_global`. Set to `FALSE` to disable global re-centering
#' for all calls.
#'
#' @section Iteration and Convergence:
#' The core engine runs a single Taylor pass around `theta_seed`. If `global_refinement = TRUE`, 
#' the seed is updated for up to `global_passes` iterations until the maximum parameter change 
#' is below `convergence_epsilon`. K-means initialization can provide alternative seeds before 
#' this loop. Tiered refinement then applies local re-centering or Gauss-Newton optimization 
#' with its own iteration limits (e.g. `refinement_thresholds$gauss_newton_maxiter`).
#'
#' @param fmri_data An fMRI dataset object or numeric matrix (timepoints x voxels)
#' @param event_model Event timing design matrix or event model object
#' @param parametric_model Character string specifying HRF model (currently "lwu")
#' @param theta_seed Initial parameters or "data_driven" for automatic selection
#' @param theta_bounds List with 'lower' and 'upper' bounds for parameters
#' @param confound_formula Optional formula for nuisance regressors (deprecated, use baseline_model)
#' @param baseline_model Either an fmrireg::baseline_model object for complex confound regression,
#'   or NULL (default) for no confound regression. When provided as baseline_model object,
#'   supports drift modeling, motion parameters, and other nuisance regressors.
#' @param hrf_eval_times Time points for HRF evaluation
#' @param hrf_span Duration for HRF evaluation (default 30 seconds)
#' @param lambda_ridge Ridge penalty for stability (default 0.01)
#' @param mask Optional voxel selection mask
#' @param global_refinement Logical: perform iterative global refinement? (Sprint 2)
#' @param global_passes Number of global refinement iterations (default 3)
#' @param convergence_epsilon Convergence criterion for global refinement (default 0.01)
#' @param kmeans_refinement Logical: perform K-means local refinement? (Sprint 3)
#' @param kmeans_k Number of clusters for K-means (default 5)
#' @param kmeans_passes Number of K-means refinement passes (default 2)
#' @param tiered_refinement Character: refinement strategy ("none", "moderate", "aggressive")
#' @param refinement_thresholds List of R^2 and SE thresholds for tiered refinement
#' @param parallel Logical: use parallel processing?
#' @param n_cores Number of cores for parallel processing (NULL = auto-detect)
#' @param compute_se Logical: compute standard errors?
#' @param safety_mode Character: "maximum", "balanced", or "performance".
#'   Determines the level of input validation. "maximum" triggers
#'   comprehensive checks, "balanced" performs standard checks, and
#'   "performance" uses minimal validation.
#' @param progress Logical: show progress bar?
#' @param verbose Logical: print detailed messages?
#'
#' @return Object of class 'parametric_hrf_fit' containing:
#'   - estimated_parameters: Matrix of HRF parameters (voxels x parameters)
#'   - amplitudes: Response amplitudes for each voxel
#'   - standard_errors: Parameter standard errors (if computed)
#'   - fit_quality: List with R-squared and other metrics
#'   - convergence_info: Detailed convergence information
#'   - refinement_info: Information about refinement strategies applied
#'   - metadata: Complete analysis metadata
#'
#' @examples
#' # Basic example
#' set.seed(1)
#' fmri_data <- matrix(rnorm(40), nrow = 20, ncol = 2)
#' events <- matrix(0, nrow = 20, ncol = 1)
#' events[c(5, 15), 1] <- 1
#' fit <- estimate_parametric_hrf(fmri_data, events, parametric_model = "lwu",
#'                                verbose = FALSE)
#' summary(fit)
#' 
#' # Using data-driven initialization
#' fit_dd <- estimate_parametric_hrf(fmri_data, events, 
#'                                   theta_seed = "data_driven",
#'                                   verbose = FALSE)
#' 
#' # With refinement options
#' fit_refined <- estimate_parametric_hrf(fmri_data, events,
#'                                        global_refinement = TRUE,
#'                                        tiered_refinement = "moderate",
#'                                        verbose = FALSE)
#' 
#' \dontrun{
#' # Using baseline_model for complex confound regression
#' library(fmrireg)
#' sframe <- sampling_frame(blocklens = c(150, 150), TR = 2)
#' bmodel <- baseline_model(basis = "bs", degree = 5, sframe = sframe)
#' 
#' fit_baseline <- estimate_parametric_hrf(fmri_data, events,
#'                                         baseline_model = bmodel,
#'                                         verbose = FALSE)
#' }
#'
#' @export
estimate_parametric_hrf <- function(
  fmri_data,
  event_model,
  parametric_model = "lwu",
  # Basic parameters
  theta_seed = NULL,
  theta_bounds = NULL,
  confound_formula = NULL,
  baseline_model = NULL,
  hrf_eval_times = NULL,
  hrf_span = 30,
  lambda_ridge = 0.01,
  mask = NULL,
  # Global refinement
  global_refinement = TRUE,
  global_passes = 3,
  convergence_epsilon = 0.01,
  # K-means refinement
  kmeans_refinement = FALSE,
  kmeans_k = 5,
  kmeans_passes = 2,
  # Tiered refinement
  tiered_refinement = c("none", "moderate", "aggressive"),
  refinement_thresholds = list(
    r2_easy = 0.7,
    r2_hard = 0.3,
    se_low = 0.3,
    se_high = 0.7,
    gauss_newton_maxiter = 10
  ),
  # Parallel processing
  parallel = FALSE,
  n_cores = NULL,
  # Output options
  compute_se = TRUE,
  # Safety and diagnostics
  safety_mode = c("balanced", "maximum", "performance"),
  progress = TRUE,
  verbose = TRUE
) {
  
  
  # Ensure matrices are double-precision for C++ backend
  if (is.matrix(fmri_data) && typeof(fmri_data) != "double") {
    storage.mode(fmri_data) <- "double"
  }
  if (is.matrix(event_model) && typeof(event_model) != "double") {
    storage.mode(event_model) <- "double"
  }
  
  # Prepare arguments list for reuse
  args <- list(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_model = parametric_model,
    theta_seed = theta_seed,
    theta_bounds = theta_bounds,
    confound_formula = confound_formula,
    baseline_model = baseline_model,
    hrf_eval_times = hrf_eval_times,
    hrf_span = hrf_span,
    lambda_ridge = lambda_ridge,
    mask = mask,
    global_refinement = global_refinement,
    global_passes = global_passes,
    convergence_epsilon = convergence_epsilon,
    kmeans_refinement = kmeans_refinement,
    kmeans_k = kmeans_k,
    kmeans_passes = kmeans_passes,
    tiered_refinement = tiered_refinement,
    refinement_thresholds = refinement_thresholds,
    parallel = parallel,
    n_cores = n_cores,
    compute_se = compute_se,
    safety_mode = safety_mode,
    progress = progress,
    verbose = verbose
  )
  
  
  # DEFAULT: Use refactored implementation
  do.call(.run_hrf_estimation_engine, args)
}