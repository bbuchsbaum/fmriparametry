# Global Code Audit

This audit reviews the repository for redundancy, lack of modularity, and potential cruft that could lead to bugs.

## Redundant or Extraneous Files

- `..Rcheck/00check.log` is an R CMD check log that should not be version controlled.
- Two files (`R/performance_enhancements.R` and `R/performance_optimizations.R`) both define `.qr_cache` and `.cached_qr_solve`, creating duplicate functionality.
- `tests/testthat/test-placeholder.R` contains a trivial placeholder test.

## Lack of Modularity

- `R/estimate_parametric_hrf.R` is over 1000 lines long, handling all stages of the estimation workflow in one function. Splitting this into smaller helpers would improve maintainability.
- Performance-related helpers are scattered across multiple files (`performance_enhancements.R`, `performance_optimizations.R`, `smart_performance_dispatcher.R`), some of which overlap in functionality.

## TODOs / Incomplete Sections

- `R/estimate_parametric_hrf.R` contains TODO comments for computing final R-squared and RMSE values.
- `R/test_compatibility_layer.R` includes a TODO for residual computation.

## Miscellaneous Cruft

- Multiple standalone R scripts in the repository root (`engineering_excellence_demo.R`, `test_basic.R`, etc.) are for manual testing or benchmarking. While excluded from the R package build via `.Rbuildignore`, keeping them in the repo could confuse users.
- Placeholder comments remain in `tests/testthat/test-integration.R` for future features.

## Recommendations

1. Remove `..Rcheck` from version control and add it to `.gitignore`.
2. Consolidate duplicate implementations of `.cached_qr_solve` into a single module to avoid inconsistencies.
3. Break down `estimate_parametric_hrf()` into smaller functions (data prep, fitting, refinement, inference).
4. Address remaining TODOs to ensure result metrics are accurately computed.
5. Consider pruning or moving example scripts to a dedicated `examples/` directory to keep the top level clean.
6. Remove placeholder tests or flesh them out with real assertions.

