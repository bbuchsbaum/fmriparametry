# Refactoring Plan: fmriparametry

## Executive Summary

This document outlines a comprehensive refactoring plan for the fmriparametry R package, addressing significant technical debt including code duplication, architectural issues, and bugs. The plan follows a phased approach to ensure safety and maintain functionality throughout the process.

## 1. High-Level Goals

- **Improve maintainability**: Consolidate duplicate code and establish clear module boundaries
- **Fix known bugs**: Particularly the R-squared calculation issue
- **Simplify architecture**: Refactor the monolithic `estimate_parametric_hrf()` function
- **Establish consistent patterns**: Standardize error handling, validation, and naming conventions
- **Create robust testing**: Build comprehensive test suite before refactoring
- **Remove dead code**: Eliminate unused functions and compatibility layers

## 2. Phase 0: Preparation & Tooling

### Setup and Analysis
- [ ] Create dedicated `refactor` branch for all work
- [ ] Configure development environment with proper tooling
  - [ ] Set up `lintr` configuration for code style checking
  - [ ] Configure `styler` for automatic formatting
  - [ ] Ensure `testthat` is properly configured
- [ ] Run dependency analysis to identify dead code
  - [ ] Use custom R script to trace function calls from `estimate_parametric_hrf`
  - [ ] Document suspected dead code in `analysis/dead_code_analysis.md`
  - [ ] Mark (don't delete) potentially dead functions with deprecation warnings

### Test Infrastructure
- [ ] Create standardized test datasets in `tests/testthat/testdata/`:
  - [ ] `minimal_2voxel_known_hrf.rds`: 2 voxels, 10 timepoints, known parameters
  - [ ] `realistic_100voxel_varying_snr.rds`: 100 voxels, 500 timepoints, varying SNR
  - [ ] `edge_cases/`: Directory with various edge case datasets
    - [ ] All-zero time series
    - [ ] Time series with NA/NaN values
    - [ ] Single voxel data
    - [ ] Singular design matrices
- [ ] Create `scripts/generate_test_data.R` for reproducible test data generation

### Documentation
- [ ] Document current public API surface in `docs/current_api.md`
- [ ] Establish naming conventions (see Section 7)

## 3. Phase 1: Build Safety Net & Fix Critical Bugs

### Characterization Tests
- [ ] Write comprehensive tests for current behavior of `estimate_parametric_hrf()`
  - [ ] Happy path with known simple dataset
  - [ ] Edge cases (zeros, NAs, single voxel)
  - [ ] Parameter variations (with/without temporal derivatives)
  - [ ] Different HRF models and basis functions
- [ ] Mark tests that validate buggy behavior with TODOs

### Bug Fixes
- [ ] Fix R-squared recalculation bug
  - [ ] Add test case demonstrating the bug
  - [ ] Implement fix to recalculate R-squared after refinement
  - [ ] Update tests to assert correct behavior
- [ ] Document any other critical bugs discovered during testing

### API Documentation
- [ ] Complete roxygen2 documentation for all public functions
- [ ] Generate initial pkgdown site to review documentation

## 4. Phase 2: Consolidate Utilities ("Make the Change Easy")

### Create New Structure
- [ ] Create consolidated utility files:
  - [ ] `R/utils-validation.R`: Input validation and checks
  - [ ] `R/utils-numerical.R`: Safe solvers, matrix operations
  - [ ] `R/utils-memory.R`: Chunked processing, memory management
  - [ ] `R/utils-parallel.R`: Parallel backend management
  - [ ] `R/utils-performance.R`: Timing and profiling utilities

### Migration Process
For each utility function:
1. Identify function in current location
2. Write focused unit test in `tests/testthat/test-utils-*.R`
3. Move function to appropriate `utils-*.R` file
4. Replace original with call to new location
5. Run full test suite
6. Commit

### Specific Consolidations
- [ ] Merge all versions of `.chunked_processing()` into canonical implementation
- [ ] Consolidate QR solving functions (`.cached_qr_solve()`, `.smart_qr_solve()`, etc.)
- [ ] Unify parallel processing helpers
- [ ] Standardize validation functions

## 5. Phase 3: Refactor Main Function

### Configuration Object
- [ ] Create `hrf_config` S3 class to encapsulate parameters
  - [ ] Constructor: `hrf_config()` with validation
  - [ ] Group related parameters (refinement_opts, parallel_opts, output_opts)
  - [ ] Print method for clear parameter display

### Extract Helper Functions
- [ ] Create `R/estimate-helpers.R` for orchestration helpers
- [ ] Extract major blocks into functions:
  - [ ] `.prepare_estimation()`: Initial setup and validation
  - [ ] `.initialize_parameters()`: Seed generation and initialization
  - [ ] `.run_refinement_pipeline()`: Main refinement logic
  - [ ] `.finalize_results()`: Post-processing and output packaging

### Refactor Refinement Logic
- [ ] Move refinement functions to dedicated files:
  - [ ] `R/refine-kmeans.R`: K-means initialization
  - [ ] `R/refine-local.R`: Local recentering
  - [ ] `R/refine-gauss-newton.R`: Already exists, may need updates

### Simplify Main Function
- [ ] Reduce `estimate_parametric_hrf()` to high-level orchestrator
- [ ] Target: < 50 lines focusing on workflow coordination

## 6. Phase 4: Polish and Finalize

### Dead Code Removal
- [ ] Remove all functions marked as deprecated in Phase 0
- [ ] Delete compatibility layers (e.g., `test_compatibility_layer.R`)
- [ ] Remove duplicate implementations identified in Phase 2
- [ ] Use `covr::package_coverage()` to find untested code

### Standardize Conventions
- [ ] Run `styler::style_pkg()` for consistent formatting
- [ ] Run `lintr::lint_package()` and fix issues
- [ ] Standardize error handling:
  - [ ] Use `rlang::abort()` for errors
  - [ ] Use `rlang::warn()` for warnings
  - [ ] Use `rlang::inform()` for messages
- [ ] Review and improve all function/variable names

### Performance Optimization
- [ ] Replace magic numbers with named constants or options
- [ ] Implement proper FFT crossover logic based on complexity
- [ ] Add package options for performance tuning

### Final Documentation
- [ ] Update all roxygen2 documentation
- [ ] Write vignettes for major use cases
- [ ] Generate final pkgdown site
- [ ] Update README with current examples

## 7. Naming Conventions

### Files
- **R files**: `snake_case.R` (e.g., `estimate_parametric_hrf.R`)
- **Test files**: `test-snake-case.R` (e.g., `test-utils-validation.R`)

### Functions
- **Public functions**: `verb_noun()` (e.g., `estimate_parametric_hrf()`, `validate_inputs()`)
- **Internal functions**: `.verb_noun()` (e.g., `.prepare_parametric_inputs()`)
- **S3 methods**: `generic.class()` (e.g., `print.parametric_hrf_fit()`)

### Variables
- **Regular variables**: `snake_case` (e.g., `design_matrix`, `n_voxels`)
- **Constants**: `UPPER_SNAKE_CASE` (e.g., `DEFAULT_LAMBDA`, `MAX_ITERATIONS`)

### Classes
- **S3 classes**: `snake_case` (e.g., `parametric_hrf_fit`, `hrf_config`)

## 8. Decision Log

### Test Data Strategy
**Decision**: Create three categories of standardized test data (minimal, realistic, edge cases) to ensure reproducible testing.
**Rationale**: Random or poorly understood test data makes it impossible to distinguish real regressions from noise.

### Parallel Code Paths
**Decision**: Default to one canonical implementation, use Strategy Pattern only where clear trade-offs exist.
**Process**:
1. Investigate why each parallel path exists
2. If one is simply older/worse, eliminate it
3. If legitimate trade-off exists (speed vs accuracy), implement Strategy Pattern
4. Document decision for each case

### Error Handling
**Decision**: Standardize on rlang for all condition handling.
**Rationale**: Provides richer error information and better debugging experience.

### Configuration Management
**Decision**: Use S3 configuration object instead of 20+ individual parameters.
**Rationale**: Simplifies function signatures and makes parameter grouping explicit.

## 9. Progress Tracking

### Pull Requests
- [ ] PR #1: Phase 0 - Setup and initial analysis
- [ ] PR #2: Phase 1 - Tests and R-squared bug fix
- [ ] PR #3: Phase 2 - Utility consolidation
- [ ] PR #4: Phase 3 - Main function refactoring
- [ ] PR #5: Phase 4 - Final cleanup

### Milestones
- [ ] All tests passing with fixed R-squared calculation
- [ ] All utilities consolidated with no duplication
- [ ] Main function < 50 lines
- [ ] 100% documentation coverage
- [ ] All dead code removed

## 10. Risk Mitigation

### Backward Compatibility
- Maintain existing public API throughout refactoring
- Use deprecation warnings for any necessary breaking changes
- Provide migration guide if API changes are required

### Testing Strategy
- Never commit code that breaks existing tests
- Add new tests before refactoring each component
- Use continuous integration to catch regressions

### Rollback Plan
- All work on `refactor` branch
- Can always return to `main` if issues arise
- Tag stable checkpoints throughout process

---

*Last Updated: [Date]*
*Next Review: [Date]*