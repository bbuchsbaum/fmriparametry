#' DEPRECATED VERSION NOTICE
#'
#' @description
#' Warning about deprecated function versions
#'
#' @details  
#' The following files contain DEPRECATED versions of estimate_parametric_hrf:
#' \itemize{
#'   \item estimate_parametric_hrf_v1_deprecated.R (original Sprint 1 version)
#'   \item estimate_parametric_hrf_v2.R (Sprint 2 features - NOW MERGED)
#'   \item estimate_parametric_hrf_v3.R (Sprint 3 features - NOW MERGED)
#'   \item estimate_parametric_hrf_rock_solid.R (safety features - NOW MERGED)
#' }
#'
#' ALL features have been consolidated into the ONE TRUE version in:
#' estimate_parametric_hrf.R
#'
#' Having multiple functions with the same name is UNIMPECCABLE.
#' We have resolved this engineering malpractice.
#'
#' DO NOT USE THE DEPRECATED VERSIONS.
#' They remain only for historical reference.
#'
#' @keywords internal

.deprecated_version_warning <- function() {
  warning(
    "You are using a DEPRECATED version of estimate_parametric_hrf. ",
    "Please use the consolidated version in estimate_parametric_hrf.R",
    call. = FALSE
  )
}

#' Sprint 2 wrapper
#'
#' Provides backward compatibility with the old Sprint 2 API. It simply
#' calls [estimate_parametric_hrf()] with the supplied arguments and
#' issues a deprecation warning.
#'
#' @inheritParams estimate_parametric_hrf
#' @return A `parametric_hrf_fit` object.
#' @export
estimate_parametric_hrf_v2 <- function(...) {
  .deprecated_version_warning()
  estimate_parametric_hrf(...)
}

#' Sprint 3 wrapper
#'
#' Backwards compatible wrapper for the Sprint 3 implementation.
#' Delegates to [estimate_parametric_hrf()] after issuing a deprecation
#' warning.
#'
#' @inheritParams estimate_parametric_hrf
#' @return A `parametric_hrf_fit` object.
#' @export
estimate_parametric_hrf_v3 <- function(...) {
  .deprecated_version_warning()
  estimate_parametric_hrf(...)
}

#' Rock-solid wrapper
#'
#' Legacy entry point that mapped to the robust "rock solid" estimator.
#' The functionality is now incorporated into [estimate_parametric_hrf()].
#'
#' @inheritParams estimate_parametric_hrf
#' @return A `parametric_hrf_fit` object.
#' @export
estimate_parametric_hrf_rock_solid <- function(...) {
  .deprecated_version_warning()
  estimate_parametric_hrf(...)
}