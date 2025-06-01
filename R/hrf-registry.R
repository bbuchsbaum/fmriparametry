# HRF model plugin registry

# Environment storing available HRF interfaces
.hrf_registry <- new.env(parent = emptyenv())

# Register default LWU interface on load
.hrf_registry$lwu <- list(
  hrf_function = .lwu_hrf_function,
  taylor_basis = .lwu_hrf_taylor_basis_function,
  parameter_names = .lwu_hrf_parameter_names(),
  default_seed = .lwu_hrf_default_seed,
  default_bounds = .lwu_hrf_default_bounds
)

#' Register a custom HRF model
#'
#' Adds a user supplied HRF interface to the registry so it can be used
#' by `estimate_parametric_hrf` via the `parametric_hrf` argument.
#'
#' @param name Character name of the HRF model
#' @param interface List with elements `hrf_function`, `taylor_basis`,
#'   `parameter_names`, `default_seed` and `default_bounds`.
#' @export
register_hrf_model <- function(name, interface) {
  stopifnot(is.character(name), length(name) == 1)
  stopifnot(is.list(interface))
  .hrf_registry[[name]] <- interface
  invisible(TRUE)
}

# Helper to retrieve a registered interface
.get_hrf_interface <- function(name) {
  if (!exists(name, envir = .hrf_registry, inherits = FALSE)) {
    stop(sprintf("HRF model '%s' is not registered", name), call. = FALSE)
  }
  get(name, envir = .hrf_registry, inherits = FALSE)
}
