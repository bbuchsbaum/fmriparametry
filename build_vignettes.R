#!/usr/bin/env Rscript

# Build vignettes for fmriparametric package
# Run this script from the package root directory

cat("Building fmriparametric vignettes...\n")

# Check if we're in the right directory
if (!file.exists("DESCRIPTION")) {
  stop("Please run this script from the fmriparametric package root directory")
}

# Check for required packages
required_pkgs <- c("knitr", "rmarkdown", "devtools")
missing_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[,"Package"]]

if (length(missing_pkgs) > 0) {
  cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs)
}

# Load devtools
library(devtools)

# Build vignettes
cat("\nBuilding vignettes...\n")
tryCatch({
  build_vignettes()
  cat("\nVignettes built successfully!\n")
  cat("\nTo view the vignettes, use:\n")
  cat("  browseVignettes('fmriparametric')\n")
  cat("\nOr view individual vignettes:\n")
  cat("  vignette('quick-start', package = 'fmriparametric')\n")
  cat("  vignette('parametric-hrf-estimation', package = 'fmriparametric')\n")
  cat("  vignette('technical-details', package = 'fmriparametric')\n")
}, error = function(e) {
  cat("\nError building vignettes:\n")
  cat(e$message, "\n")
  cat("\nTry building manually with:\n")
  cat("  devtools::build_vignettes()\n")
})
