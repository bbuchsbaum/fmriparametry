#' Deprecated version notice
#'
#' The following files contain obsolete versions of `estimate_parametric_hrf`:
#' \itemize{
#'   \item estimate_parametric_hrf_v1_deprecated.R
#'   \item estimate_parametric_hrf_v2.R
#'   \item estimate_parametric_hrf_v3.R
#'   \item estimate_parametric_hrf_rock_solid.R
#' }
#'
#' All features have been merged into `estimate_parametric_hrf.R`. These files
#' remain for historical reference only.
#'
#' @keywords internal

.deprecated_version_warning <- function() {
  warning(
    "You are using a DEPRECATED version of estimate_parametric_hrf. ",
    "Please use the consolidated version in estimate_parametric_hrf.R",
    call. = FALSE
  )
}