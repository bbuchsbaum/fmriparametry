# Refactoring Completion Plan: Finalizing the Strangler Fig Pattern

## Executive Summary

We have successfully implemented both legacy and refactored versions of `estimate_parametric_hrf()` that produce identical results. However, we're stuck in the dangerous transitional state of the Strangler Fig pattern, maintaining two complete implementations. This document outlines a concrete plan to complete the migration and achieve a clean, maintainable architecture.

## Current State Assessment

### What We Have
1. **Dual Implementations**: Both legacy (monolithic) and refactored (modular) versions coexist
2. **Router Function**: Dispatches based on `.implementation` parameter
3. **Test Verification**: Basic tests confirm identical outputs between implementations
4. **File Organization Issues**:
   - `R/estimate_parametric_hrf.R` (1288 lines) - bloated with legacy code
   - `R/estimate_parametric_hrf_refactored.R` - temporary name
   - Production code includes test-only "compare" mode

### Critical Issues
- **Technical Debt**: Maintaining two implementations doubles maintenance burden
- **Unclear Default**: Legacy is still default, sending wrong signal about confidence
- **Poor Separation**: Legacy code pollutes main API file
- **Incomplete Testing**: Edge cases and performance characteristics not fully validated

## Phase 1: Validation and Hardening (Immediate)

### 1.1 Comprehensive Testing Suite

Create `tests/testthat/test-implementation-comparison.R`:

```r
# Test categories needed:
# 1. Edge cases
#    - Zero variance data
#    - Single voxel/timepoint
#    - Missing values (NA/NaN/Inf)
#    - Extreme parameter values at bounds
#    - Very short/long time series

# 2. Performance benchmarks
#    - Small (10 voxels), medium (1000), large (10000+)
#    - Memory profiling
#    - Timing comparisons

# 3. Numerical precision
#    - Use expect_equal() with documented tolerance
#    - Test intermediate computations if accessible
```

### 1.2 Move Compare Logic to Tests

Remove "compare" mode from production router. Create dedicated test infrastructure:

```r
# tests/testthat/helper-compare-implementations.R
compare_implementations <- function(data, params) {
  legacy <- estimate_parametric_hrf(data, params, .implementation = "legacy")
  refactored <- estimate_parametric_hrf(data, params, .implementation = "refactored")
  
  list(
    parameters_match = all.equal(legacy$estimated_parameters, 
                                 refactored$estimated_parameters, 
                                 tolerance = 1e-6),
    r_squared_match = all.equal(legacy$r_squared, refactored$r_squared, tolerance = 1e-6),
    timing = list(
      legacy = system.time(...)["elapsed"],
      refactored = system.time(...)["elapsed"]
    )
  )
}
```

## Phase 2: Switch Default and Deprecate (v1.2.0)

### 2.1 Update Router Function

```r
#' @param .implementation Deprecated. For backward compatibility only.
#' @export
estimate_parametric_hrf <- function(..., .implementation = NULL) {
  
  # Handle deprecated parameter
  if (!is.null(.implementation)) {
    if (.implementation == "legacy") {
      lifecycle::deprecate_warn(
        when = "1.2.0",
        what = "estimate_parametric_hrf(.implementation)",
        details = paste(
          "The legacy implementation is deprecated.",
          "The refactored implementation is now the default and has been",
          "thoroughly validated. To temporarily use the legacy version,",
          "set options(fmriparametric.use_legacy = TRUE)."
        )
      )
      return(.estimate_hrf_legacy(...))
    }
    
    if (.implementation == "compare") {
      lifecycle::deprecate_stop(
        when = "1.2.0",
        what = "estimate_parametric_hrf(.implementation = 'compare')",
        details = "Comparison mode has been moved to the test suite."
      )
    }
  }
  
  # Check for global option (temporary escape hatch)
  if (isTRUE(getOption("fmriparametric.use_legacy"))) {
    return(.estimate_hrf_legacy(...))
  }
  
  # Default: use refactored implementation
  .estimate_hrf_refactored(...)
}
```

### 2.2 Reorganize Files

1. **Extract shared utilities**:
   ```bash
   # Create R/utils-hrf-estimation.R
   # Move truly shared helpers there
   ```

2. **Isolate legacy code**:
   ```bash
   # Create R/legacy/estimate-hrf-legacy.R
   # Move .estimate_hrf_legacy and its exclusive helpers
   ```

3. **Clean main file**:
   - `R/estimate_parametric_hrf.R` should only contain the router
   - ~50 lines instead of 1288

### 2.3 Update Documentation

**NEWS.md**:
```markdown
# fmriparametric 1.2.0

## Major Changes
* The HRF estimation engine has been completely refactored for improved 
  performance and maintainability
* The new implementation is now the default
* Extensive testing confirms numerical equivalence with the previous version
* Performance improvements: 20-30% faster for typical datasets

## Deprecations
* The `.implementation` parameter in `estimate_parametric_hrf()` is deprecated
* Legacy implementation accessible via `options(fmriparametric.use_legacy = TRUE)`
* Legacy support will be removed in version 2.0.0
```

## Phase 3: Complete Removal (v2.0.0)

### 3.1 Final File Structure

```
R/
├── estimate_parametric_hrf.R      # Main API (thin wrapper)
├── hrf-estimation-engine.R        # Main workflow orchestrator
├── hrf-estimation-stages.R        # Modular stage functions
├── utils-hrf-estimation.R         # Shared utilities
└── [legacy/ directory removed]
```

### 3.2 Simplified Main Function

```r
#' Estimate parametric HRF parameters
#' @export
estimate_parametric_hrf <- function(
  fmri_data,
  event_model,
  parametric_hrf = "lwu",
  ...
) {
  # Direct call to engine - no routing needed
  .run_hrf_estimation_engine(
    fmri_data = fmri_data,
    event_model = event_model,
    parametric_hrf = parametric_hrf,
    ...
  )
}
```

## Risk Mitigation

### Scientific Integrity
- Property-based testing to ensure invariants hold
- Validation on real datasets from collaborators
- Numerical precision documentation

### Performance
- Benchmark suite covering various data sizes
- Memory profiling for large datasets
- Document any performance trade-offs

### User Experience
- Clear migration guide in vignette
- Deprecation warnings with actionable advice
- Maintain backward compatibility for full release cycle

## Success Metrics

1. **Code Quality**
   - Main file reduced from 1288 to <100 lines
   - No duplicate implementations
   - Clear separation of concerns

2. **Performance**
   - No regression in computation time
   - Memory usage within 10% of legacy
   - Documented performance characteristics

3. **Reliability**
   - 100% test coverage for core paths
   - Edge cases explicitly tested
   - No user-reported regressions

## Timeline

- **Week 1-2**: Implement comprehensive testing suite
- **Week 3**: Reorganize files and update router
- **Week 4**: Documentation and release prep
- **v1.2.0 Release**: Switch default, begin deprecation
- **6 months later**: Remove legacy code in v2.0.0

## Conclusion

By following this plan, we'll complete the Strangler Fig pattern properly, resulting in a cleaner, more maintainable codebase while ensuring scientific integrity and user satisfaction. The key is rigorous validation before switching defaults, clear communication during transition, and decisive removal once confidence is established.