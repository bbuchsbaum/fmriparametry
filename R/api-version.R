#' Public API version
#'
#' This function returns the semantic version of the package, which also
#' serves as the version for the exported API. It provides a reliable way
#' for compatibility checks by reading the version directly from the
#' package's DESCRIPTION file.
#'
#' @return A `character` string representing the package version (e.g., "0.2.0").
#' @export
#' @examples
#' if (fmriparametric_api_version() < "0.3.0") {
#'   warning("This script was designed for a newer version of fmriparametry.")
#' }
fmriparametric_api_version <- function() {
  utils::packageDescription("fmriparametry")$Version
}

