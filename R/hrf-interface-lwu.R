#' Internal LWU HRF function wrapper
#'
#' Fast LWU HRF formula without external assertion overhead.
#'
#' @keywords internal
.lwu_hrf_formula <- function(t, tau, sigma, rho, normalize = "none") {
  if (!is.numeric(t)) {
    stop("t must be numeric", call. = FALSE)
  }
  if (!is.numeric(tau) || length(tau) != 1L) {
    stop("tau must be a single numeric value", call. = FALSE)
  }
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
    stop("sigma must be a positive numeric scalar", call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1L) {
    stop("rho must be a single numeric value", call. = FALSE)
  }
  if (!normalize %in% c("none", "height", "area")) {
    stop("normalize must be one of 'none', 'height', or 'area'", call. = FALSE)
  }

  if (normalize == "area") {
    warning(
      "`normalize = \"area\"` is not yet fully implemented for hrf_lwu and will behave like `normalize = \"none\"`.",
      call. = FALSE
    )
  }

  # .Call overhead can dominate for short vectors; keep a tiny-vector fast path.
  if (length(t) <= 128L && normalize %in% c("none", "height")) {
    response <- exp(-((t - tau)^2) / (2 * sigma^2)) -
      rho * exp(-((t - (tau + 2 * sigma))^2) / (2 * (1.6 * sigma)^2))

    if (normalize == "height") {
      max_abs_val <- max(abs(response), na.rm = TRUE)
      if (max_abs_val > 1e-10) {
        response <- response / max_abs_val
      }
    }

    return(response)
  }

  lwu_hrf_formula_cpp(
    t = t,
    tau = as.numeric(tau),
    sigma = as.numeric(sigma),
    rho = as.numeric(rho),
    normalize = normalize
  )
}

#' Internal LWU HRF function wrapper
#'
#' @param t Numeric vector of time points
#' @param params_vector Numeric vector of LWU parameters (tau, sigma, rho)
#' @param bounds List with 'lower' and 'upper' numeric vectors
#' @param ... Additional arguments passed to `fmrireg::hrf_lwu`
#' @return Numeric vector of HRF values
#' @keywords internal
.lwu_hrf_function <- function(t, params_vector, bounds, ...) {
  if (!is.numeric(params_vector) || length(params_vector) != 3) {
    stop("length(params_vector) not equal to 3", call. = FALSE)
  }
  if (is.null(bounds) || is.null(bounds$lower) || is.null(bounds$upper)) {
    stop("Bounds must be provided to .lwu_hrf_function", call. = FALSE)
  }
  
  # Enforce parameter bounds to prevent fmrireg errors
  params_vector <- .apply_bounds(params_vector, bounds)
  
  # Ensure sigma is strictly > 0.05 as required by fmrireg
  params_vector[2] <- max(params_vector[2], 0.051)

  dots <- list(...)
  normalize <- dots$normalize %||% "none"
  if (!normalize %in% c("none", "height", "area")) {
    stop("normalize must be one of 'none', 'height', or 'area'", call. = FALSE)
  }

  .lwu_hrf_formula(
    t = t,
    tau = params_vector[1],
    sigma = params_vector[2],
    rho = params_vector[3],
    normalize = normalize
  )
}

# Fast finite-difference LWU Taylor basis used in the core engine path.
.lwu_hrf_taylor_basis_fast <- function(params_vector0, t_hrf_eval, bounds, ...) {
  if (!is.numeric(params_vector0) || length(params_vector0) != 3) {
    stop("length(params_vector0) not equal to 3", call. = FALSE)
  }
  if (is.null(bounds) || is.null(bounds$lower) || is.null(bounds$upper)) {
    stop("Bounds must be provided to .lwu_hrf_taylor_basis_fast", call. = FALSE)
  }

  basis <- lwu_taylor_basis_fd_cpp(
    params_vector0 = as.numeric(params_vector0),
    t_hrf_eval = as.numeric(t_hrf_eval),
    lower = as.numeric(bounds$lower),
    upper = as.numeric(bounds$upper),
    rel_step = 1e-4,
    min_step = 1e-5,
    bound_eps = 1e-8
  )
  storage.mode(basis) <- "double"
  basis
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
  if (!is.numeric(params_vector0) || length(params_vector0) != 3) {
    stop("length(params_vector0) not equal to 3", call. = FALSE)
  }
  if (is.null(bounds) || is.null(bounds$lower) || is.null(bounds$upper)) {
    stop("Bounds must be provided to .lwu_hrf_taylor_basis_function", call. = FALSE)
  }

  params_vector0 <- as.numeric(params_vector0)
  in_bounds <- all(params_vector0 >= bounds$lower & params_vector0 <= bounds$upper)

  if (in_bounds) {
    # In-bounds path: use fast C++ finite-difference basis (used by engine and
    # matched by Rcpp kernel tests).
    return(.lwu_hrf_taylor_basis_fast(params_vector0, t_hrf_eval, bounds, ...))
  }

  # Out-of-bounds compatibility path: preserve legacy clamp semantics used by
  # edge-case wrapper tests.
  theta0 <- .apply_bounds(
    params = params_vector0,
    bounds = bounds,
    epsilon = c(0.01, 0.01, 0.01)
  )
  names(theta0) <- c("tau", "sigma", "rho")

  basis <- fmrihrf::hrf_basis_lwu(
    theta0 = theta0,
    t = t_hrf_eval,
    normalize_primary = "none",
    ...
  )
  if (!is.matrix(basis)) {
    basis <- matrix(basis, ncol = 4)
  }
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
