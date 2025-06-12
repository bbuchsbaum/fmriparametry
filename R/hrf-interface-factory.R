#' Create an HRF Interface
#'
#' Factory function to generate the hrf_interface list for a given HRF model.
#' This approach uses closures to bake model-specific information, like parameter
#' bounds, directly into the interface functions.
#'
#' @param model Character string, e.g., "lwu".
#' @param user_bounds A list with `lower` and `upper` numeric vectors. If NULL,
#'   model-specific defaults are used.
#' @return A list containing hrf_function, taylor_basis, parameter_names, etc.
#' @keywords internal
.create_hrf_interface <- function(model, user_bounds = NULL) {
  if (model == "lwu") {
    # 1. Determine the definitive bounds to use
    default_bounds <- .lwu_hrf_default_bounds()
    
    # If user provided bounds, validate and merge with defaults
    if (!is.null(user_bounds)) {
      # Validate user bounds
      .validate_theta_bounds(user_bounds, n_params = 3)
      
      # Use user bounds but ensure sigma >= 0.051 for numerical stability
      final_bounds <- user_bounds
      if (final_bounds$lower[2] < 0.051) {
        warning("Adjusting sigma lower bound from ", final_bounds$lower[2], 
                " to 0.051 for numerical stability", call. = FALSE)
        final_bounds$lower[2] <- 0.051
      }
    } else {
      final_bounds <- default_bounds
    }
    
    # 2. Return the interface list with functions that capture `final_bounds`
    list(
      hrf_function = function(t, params_vector, ...) {
        # This function now uses the `final_bounds` from its parent environment
        .lwu_hrf_function(t, params_vector, bounds = final_bounds, ...)
      },
      taylor_basis = function(params_vector0, t_hrf_eval, ...) {
        # This function also uses the `final_bounds`
        .lwu_hrf_taylor_basis_function(params_vector0, t_hrf_eval, bounds = final_bounds, ...)
      },
      parameter_names = .lwu_hrf_parameter_names(),
      default_seed = .lwu_hrf_default_seed(),
      default_bounds = default_bounds, # Keep defaults for reference
      active_bounds = final_bounds      # Store the bounds being used
    )
  } else {
    stop(paste("Unsupported HRF model:", model), call. = FALSE)
  }
}