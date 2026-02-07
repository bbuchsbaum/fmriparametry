#' Create unified convergence configuration
#'
#' Creates a standardized convergence configuration to ensure consistent
#' convergence criteria across all refinement algorithms.
#'
#' @param param_tol Maximum parameter change tolerance (default: 1e-4)
#' @param objective_tol Relative objective function change tolerance (default: 1e-4)
#' @param max_iter Maximum iterations (default: 5)
#' @param r2_tol Minimum R-squared improvement to continue (default: 1e-5)
#'
#' @return List with convergence configuration
#' @keywords internal
.create_convergence_config <- function(
  param_tol = 1e-4,
  objective_tol = 1e-4,
  max_iter = 5,
  r2_tol = 1e-5
) {
  list(
    param_tol = param_tol,
    objective_tol = objective_tol,
    max_iter = max_iter,
    r2_tol = r2_tol
  )
}

#' Check convergence based on unified criteria
#'
#' @param current Current parameter values
#' @param previous Previous parameter values
#' @param obj_current Current objective value (optional)
#' @param obj_previous Previous objective value (optional)
#' @param config Convergence configuration from .create_convergence_config()
#'
#' @return List with converged (logical) and reason (character)
#' @keywords internal
.check_convergence <- function(current, previous, obj_current = NULL, 
                              obj_previous = NULL, config) {
  # Parameter change
  if (is.matrix(current) && is.matrix(previous)) {
    param_change <- max(abs(current - previous))
  } else {
    param_change <- max(abs(current - previous))
  }
  
  if (param_change < config$param_tol) {
    return(list(converged = TRUE, reason = "parameter_tolerance"))
  }
  
  # Objective change (if provided)
  if (!is.null(obj_current) && !is.null(obj_previous)) {
    obj_change <- abs(obj_current - obj_previous) / (abs(obj_previous) + 1e-10)
    if (obj_change < config$objective_tol) {
      return(list(converged = TRUE, reason = "objective_tolerance"))
    }
  }
  
  return(list(converged = FALSE, reason = "not_converged"))
}