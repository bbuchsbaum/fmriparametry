# Load HRF estimation implementations
# This file ensures legacy implementation is available when needed

# Load legacy implementation when needed
.load_legacy_implementation <- function() {
  if (!exists(".estimate_hrf_legacy", envir = parent.frame())) {
    # system.file() now looks for "legacy/estimate-hrf-legacy.R"
    # in the top-level installed package directory.
    legacy_file <- system.file("legacy", "estimate-hrf-legacy.R", 
                               package = "fmriparametric")
    
    if (legacy_file == "") {
      stop("Legacy implementation not found. Ensure it is located at 'inst/legacy/estimate-hrf-legacy.R'.", 
           call. = FALSE)
    }
    source(legacy_file, local = FALSE)
  }
}