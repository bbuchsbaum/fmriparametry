# ============================================================================
# ENGINEERING STANDARDS MODULE
# ============================================================================
#
# Purpose: Central repository for code quality standards and utilities
#
# This module establishes consistent engineering practices across the 
# fmriparametric package to ensure reliability, maintainability, and
# robustness of the codebase.
#
# Key principles:
# 1. Defensive programming - validate inputs early and comprehensively
# 2. Fail fast - detect problems as soon as possible
# 3. Clear error messages - help users understand and fix issues
# 4. Numerical robustness - handle edge cases and numerical issues gracefully
# 5. Consistent interfaces - standardize function signatures and returns
#
# Standards implemented:
# - Input validation with type checking and constraint enforcement
# - Standardized error messaging through diagnostics module
# - Numerical stability checks and safe operations
# - Output quality assertions
#
# Usage:
# All internal functions should use these utilities rather than
# implementing their own validation or error handling. This ensures
# consistency and makes it easier to update standards across the package.
#
# ============================================================================

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
#' @param null_ok Logical whether NULL is acceptable
#' @noRd
.validate_input <- function(x, arg_name, type = NULL, dims = NULL, 
                           constraints = NULL, null_ok = FALSE) {
  # Build spec for unified validation
  spec <- list(
    name = arg_name,
    null_ok = null_ok
  )
  
  # Handle type - convert to simpler format for .validate_data
  if (!is.null(type)) {
    # Take first type if multiple given
    spec$type <- type[1]
  }
  
  # Handle dimensions
  if (!is.null(dims) && length(dims) == 2) {
    spec$dims <- list(nrow = dims[1], ncol = dims[2])
  }
  
  # Merge constraints
  if (!is.null(constraints)) {
    spec <- c(spec, constraints)
  }
  
  # Use unified validation
  result <- .validate_data(x, spec, safety_mode = "balanced", caller = "validate_input")
  
  invisible(TRUE)
}


# Numerical Robustness --------------------------------------------------------
