#' fmriparametric: Extensible Parametric HRF Estimation for fMRI Data
#'
#' The fmriparametric package provides robust and efficient tools for estimating 
#' parameters of parametric Hemodynamic Response Function (HRF) models from fMRI 
#' data. Using an iterative linear Taylor approximation method, it enables 
#' voxel-wise estimation of interpretable HRF parameters (e.g., lag, width, 
#' undershoot amplitude) with uncertainty quantification.
#'
#' @section Main Functions:
#' The package provides the following main functions:
#' \itemize{
#'   \item \code{\link{estimate_parametric_hrf}}: Main function for parametric HRF estimation
#'   \item \code{\link{coef.parametric_hrf_fit}}: Extract parameter estimates
#'   \item \code{\link{summary.parametric_hrf_fit}}: Summarize fit results
#'   \item \code{\link{print.parametric_hrf_fit}}: Print fit summary
#' }
#'
#' @section HRF Models:
#' Currently supported parametric HRF models:
#' \itemize{
#'   \item \strong{LWU (Lag-Width-Undershoot)}: 3-parameter model with lag (τ), 
#'         width (σ), and undershoot (ρ) parameters
#' }
#'
#' @section Package Design:
#' The package is designed with extensibility in mind:
#' \itemize{
#'   \item Modular architecture for adding new HRF models
#'   \item Efficient vectorized operations for whole-brain analysis
#'   \item Integration with the fmrireg ecosystem
#'   \item Robust numerical methods with ridge regularization
#' }
#'
#' @name fmriparametric-package
#' @aliases fmriparametric
#'
#' @import methods
#' @import stats
#' @importFrom Matrix Matrix
#' @importFrom assertthat assert_that is.number is.flag
#' @importFrom stats median quantile convolve cor fft kmeans mad mvfft nextn plogis prcomp reshape runif sd time var
#' @importFrom grDevices rainbow rgb
#' @importFrom graphics abline hist legend lines mtext par
#' @importFrom utils head memory.size sessionInfo tail
#' @importFrom digest digest
#' @useDynLib fmriparametric, .registration = TRUE
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' library(fmriparametric)
#' library(fmrireg)
#' 
#' # Create simple test data
#' fmri_data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' event_data <- matrix(rbinom(100, 1, 0.1), ncol = 1)
#' 
#' # Estimate parametric HRF
#' fit <- estimate_parametric_hrf(
#'   fmri_data = fmri_data,
#'   event_model = event_data,
#'   parametric_hrf = "lwu"
#' )
#' 
#' # View results
#' print(fit)
#' summary(fit)
#' params <- coef(fit)
#' }
"_PACKAGE"

# Suppress R CMD check NOTEs about global variables
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # Variables from various functions
    ".parametric_engine_direct", "Count", "Queue", "amplitudes",
    "ff", "get_fmri_dimensions", "hrf", "load_fmri_chunk",
    "parameter", "queue", "r2", "r2_initial", "theta_current",
    "value", "voxel", "voxel_idx",
    # Variables from ggplot2 usage
    "voxel_id", "time", "response", "residual"
  ))
}

# Initialize package options
.on_load_options <- function(libname, pkgname) {
  op <- options()
  op.fmrip <- list(
    fmriparametric.refine_global = TRUE,
    fmriparametric.verbose = TRUE
  )
  toset <- !(names(op.fmrip) %in% names(op))
  if (any(toset)) options(op.fmrip[toset])
  invisible(NULL)
}