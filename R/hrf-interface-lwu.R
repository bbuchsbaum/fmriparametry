#' Internal LWU HRF function wrapper
#'
#' @param t Numeric vector of time points
#' @param params_vector Numeric vector of LWU parameters (tau, sigma, rho)
#' @param bounds List with 'lower' and 'upper' numeric vectors
#' @param ... Additional arguments passed to `fmrireg::hrf_lwu`
#' @return Numeric vector of HRF values
#' @keywords internal
.lwu_hrf_function <- function(t, params_vector, bounds, ...) {
  assertthat::assert_that(length(params_vector) == 3)
  assertthat::assert_that(!is.null(bounds), msg = "Bounds must be provided to .lwu_hrf_function")
  
  # Enforce parameter bounds to prevent fmrireg errors
  params_vector <- pmax(bounds$lower, pmin(params_vector, bounds$upper))
  
  # Ensure sigma is strictly > 0.05 as required by fmrireg
  params_vector[2] <- max(params_vector[2], 0.051)
  
  fmrihrf::hrf_lwu(t = t,
                   tau = params_vector[1],
                   sigma = params_vector[2],
                   rho = params_vector[3],
                   normalize = "none",
                   ...)
}

#' Internal LWU HRF Taylor basis wrapper
#'
#' @param params_vector0 Numeric vector of LWU parameters at expansion point
#' @param t_hrf_eval Numeric vector of time points for basis evaluation
#' @param bounds List with 'lower' and 'upper' numeric vectors
#' @param ... Additional arguments passed to `fmrireg::hrf_basis_lwu`
#' @return Matrix with length(t_hrf_eval) rows and 4 columns
#' @keywords internal
.lwu_hrf_taylor_basis_function <- function(params_vector0, t_hrf_eval, bounds, ...) {
  assertthat::assert_that(length(params_vector0) == 3)
  assertthat::assert_that(!is.null(bounds), msg = "Bounds must be provided to .lwu_hrf_taylor_basis_function")
  
  # Enforce parameter bounds to prevent fmrireg errors
  params_vector0 <- pmax(bounds$lower, pmin(params_vector0, bounds$upper))
  # Avoid exact boundary values which can break numerical derivatives
  # Use larger epsilon for rho since numDeriv uses larger step sizes
  eps <- c(0.01, 0.01, 0.01)  # Ensure enough space for numerical derivatives
  params_vector0 <- pmax(bounds$lower + eps,
                         pmin(params_vector0, bounds$upper - eps))
  
  # Ensure params_vector0 has names for fmrireg
  names(params_vector0) <- c("tau", "sigma", "rho")
  
  basis <- fmrihrf::hrf_basis_lwu(theta0 = params_vector0,
                                  t = t_hrf_eval,
                                  normalize_primary = "none",
                                  ...)
  if (!is.matrix(basis)) {
    basis <- matrix(basis, ncol = 4)
  }
  # Ensure double precision for C++ compatibility
  storage.mode(basis) <- "double"
  basis
}

#' Internal LWU HRF parameter names
#' @return Character vector of parameter names
#' @keywords internal
.lwu_hrf_parameter_names <- function() {
  c("tau", "sigma", "rho")
}

#' Internal LWU HRF default seed
#' @return Numeric vector of default parameter seed
#' @keywords internal
.lwu_hrf_default_seed <- function() {
  c(6, 2.5, 0.35)
}

#' Internal LWU HRF default bounds
#' @return List with elements `lower` and `upper`
#' @keywords internal
.lwu_hrf_default_bounds <- function() {
  list(
    lower = c(0, 0.051, 0),  # sigma must be strictly > 0.05 for fmrireg
    upper = c(20, 10, 1.5)
  )
}
