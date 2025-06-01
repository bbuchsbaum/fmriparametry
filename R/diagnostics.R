#' Diagnostics and messaging utilities
#'
#' These helper functions centralize error, warning, and information
#' messages as well as optional progress reporting. Verbosity can be
#' controlled via the option `fmriparametric.verbose`.
#'
#' @keywords internal
NULL

#' Emit an informational message respecting package verbosity
#' @keywords internal
.diag_inform <- function(msg, ...) {
  if (isTRUE(getOption("fmriparametric.verbose", TRUE))) {
    rlang::inform(sprintf(msg, ...))
  }
}

#' Emit a warning message using rlang
#' @keywords internal
.diag_warn <- function(msg, ...) {
  rlang::warn(sprintf(msg, ...))
}

#' Emit an error message using rlang
#' @keywords internal
.diag_abort <- function(msg, ...) {
  rlang::abort(sprintf(msg, ...))
}

#' Create a progressr progressor if available
#' @keywords internal
.diag_progressor <- function(steps) {
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::progressor(steps)
  } else {
    function(...) NULL
  }
}
