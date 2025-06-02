#' Engineering Standards for fmriparametric
#'
#' This file defines the engineering standards and utilities that ensure
#' code quality throughout the package.

# Validation Standards --------------------------------------------------------

#' Validate function inputs with standardized error messages
#'
#' @param x Object to validate
#' @param arg_name Character name of argument for error messages
#' @param type Expected type(s)
#' @param dims Expected dimensions (NULL to skip)
#' @param constraints List of constraints (range, properties, etc.)
#' @keywords internal
.validate_input <- function(x, arg_name, type = NULL, dims = NULL, 
                           constraints = NULL, null_ok = FALSE) {
  # NULL handling
  if (is.null(x)) {
    if (null_ok) return(invisible(TRUE))
    stop(sprintf("Argument '%s' cannot be NULL", arg_name), call. = FALSE)
  }
  
  # Type validation
  if (!is.null(type)) {
    if (!inherits(x, type)) {
      stop(sprintf(
        "Argument '%s' must be of type %s, got %s",
        arg_name, 
        paste(type, collapse = " or "),
        paste(class(x), collapse = ", ")
      ), call. = FALSE)
    }
  }
  
  # Dimension validation
  if (!is.null(dims)) {
    actual_dims <- if (is.matrix(x) || is.array(x)) dim(x) else length(x)
    
    if (length(dims) != length(actual_dims) || !all(dims == actual_dims | dims == -1)) {
      stop(sprintf(
        "Argument '%s' must have dimensions %s, got %s",
        arg_name,
        paste(dims, collapse = " x "),
        paste(actual_dims, collapse = " x ")
      ), call. = FALSE)
    }
  }
  
  # Constraint validation
  if (!is.null(constraints)) {
    .validate_constraints(x, arg_name, constraints)
  }
  
  invisible(TRUE)
}

#' Validate numerical constraints
#' @keywords internal
.validate_constraints <- function(x, arg_name, constraints) {
  # Range constraints
  if (!is.null(constraints$range)) {
    if (any(x < constraints$range[1] | x > constraints$range[2], na.rm = TRUE)) {
      stop(sprintf(
        "Argument '%s' contains values outside range [%g, %g]",
        arg_name, constraints$range[1], constraints$range[2]
      ), call. = FALSE)
    }
  }
  
  # Finite values
  if (isTRUE(constraints$finite)) {
    if (any(!is.finite(x))) {
      stop(sprintf(
        "Argument '%s' must contain only finite values",
        arg_name
      ), call. = FALSE)
    }
  }
  
  # Positive values
  if (isTRUE(constraints$positive)) {
    if (any(x <= 0, na.rm = TRUE)) {
      stop(sprintf(
        "Argument '%s' must contain only positive values",
        arg_name
      ), call. = FALSE)
    }
  }
}

# Numerical Robustness --------------------------------------------------------


# Engineering Standards Options -----------------------------------------------

#' Set engineering standards options
#' @export
set_engineering_options <- function(
  verbose = FALSE,
  validate = TRUE,
  profile = FALSE,
  precision = "double",
  debug = FALSE
) {
  options(
    fmriparametric.verbose = verbose,
    fmriparametric.validate = validate,
    fmriparametric.profile = profile,
    fmriparametric.precision = precision,
    fmriparametric.debug = debug
  )
  
  if (debug) {
    message("Engineering standards: DEBUG mode enabled")
    message("  - All validations active")
    message("  - Profiling enabled")
    message("  - Verbose output")
    options(
      fmriparametric.verbose = TRUE,
      fmriparametric.validate = TRUE,
      fmriparametric.profile = TRUE
    )
  }
  
  invisible(NULL)
}