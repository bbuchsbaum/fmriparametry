#' Internal LWU HRF function wrapper
#'
#' @param t Numeric vector of time points
#' @param params_vector Numeric vector of LWU parameters (tau, sigma, rho)
#' @param ... Additional arguments passed to `fmrireg::hrf_lwu`
#' @return Numeric vector of HRF values
#' @keywords internal
.lwu_hrf_function <- function(t, params_vector, ...) {
  assertthat::assert_that(length(params_vector) == 3)
  
  # Enforce LWU parameter bounds to prevent fmrireg errors
  bounds <- .lwu_hrf_default_bounds()
  params_vector <- pmax(bounds$lower, pmin(params_vector, bounds$upper))
  
  fmrireg::hrf_lwu(t = t,
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
#' @param ... Additional arguments passed to `fmrireg::hrf_basis_lwu`
#' @return Matrix with length(t_hrf_eval) rows and 4 columns
#' @keywords internal
.lwu_hrf_taylor_basis_function <- function(params_vector0, t_hrf_eval, ...) {
  assertthat::assert_that(length(params_vector0) == 3)
  
  # Enforce LWU parameter bounds to prevent fmrireg errors
  bounds <- .lwu_hrf_default_bounds()
  params_vector0 <- pmax(bounds$lower, pmin(params_vector0, bounds$upper))
  # Avoid exact boundary values which can break numerical derivatives
  eps <- 1e-6
  params_vector0 <- pmax(bounds$lower + eps,
                         pmin(params_vector0, bounds$upper - eps))
  
  basis <- fmrireg::hrf_basis_lwu(theta0 = params_vector0,
                                  t = t_hrf_eval,
                                  normalize_primary = "none",
                                  ...)
  if (!is.matrix(basis)) {
    basis <- matrix(basis, ncol = 4)
  }
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
    lower = c(0, 0.05, 0),
    upper = c(20, 10, 1.5)
  )
}
