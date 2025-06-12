# Edge Case Test Summary for fmriparametric

This document summarizes the comprehensive edge case tests created for the fmriparametric package, focusing on numerical stability and algorithmic robustness.

## Test Files Created

### 1. `test-gauss-newton-edge-cases.R`
Tests for the Gauss-Newton refinement algorithm and its helper functions.

**Key Edge Cases Covered:**
- `.calculate_objective_gn`:
  - Invalid theta parameters (wrong length, non-numeric)
  - Zero stimulus (singular system)
  - Near-singular HRF values
  - Normal operation validation
  
- `.get_jacobian_and_residuals`:
  - Zero stimulus returning NULL
  - Jacobian matrix validation
  - Parameters at bounds
  - Finite value checks

- `.gauss_newton_refinement`:
  - Empty refinement (no hard voxels)
  - Invalid voxel indices
  - Singular systems
  - Parameters hitting bounds
  - Line search failures
  - Numerical precision issues (small signals, collinear derivatives)

### 2. `test-parametric-engine-edge-cases.R`
Tests for the core parametric fitting engine.

**Key Edge Cases Covered:**
- Zero variance data (constant voxels)
- Near-zero variance
- Zero stimulus (no events)
- Extreme parameter values at bounds
- High-dimensional data (many voxels)
- Numerical precision limits:
  - Very large magnitude data (1e6)
  - Very small magnitude data (1e-6)
  - Mixed scales across voxels
- Sparse event designs (single event, far-apart events)
- Different baseline models (with/without intercept)

### 3. `test-numerical-stability-edge-cases.R`
Tests for numerical helper functions and core algorithms.

**Key Edge Cases Covered:**
- `.fast_batch_convolution`:
  - Empty signal
  - Single time point
  - Zero signal
  - Very long signals (FFT path)
  - Extreme values
  
- `.ridge_linear_solve`:
  - Singular matrices
  - Overdetermined systems
  - Underdetermined systems
  - Near-zero values
  - Identity matrix
  
- `.batch_convolution`:
  - Empty kernels
  - All-zero signals
  - Dimension mismatches
  
- HRF interface functions:
  - Parameters at exact bounds
  - Boundary adjustments
  - Short time vectors
  
- `.calculate_fit_metrics`:
  - Perfect fit (R² = 1)
  - Constant data with/without intercept
  - More predictors than observations
  - Pre-computed TSS
  - Negative R² clamping

### 4. `test-integration-edge-cases.R`
Integration tests for the main `estimate_parametric_hrf` function.

**Key Edge Cases Covered:**
- Pathological data patterns:
  - Perfectly periodic signals
  - Step functions
  - Extreme outliers
- Extreme parameter scenarios:
  - Very narrow bounds
  - Seed outside bounds
- Data structure edge cases:
  - Single time point (error expected)
  - More event columns than voxels
  - Very long time series
- Numerical precision limits:
  - Data near machine precision
  - Mixed scales
  - Constant columns with variable data
- Refinement integration:
  - Mix of easy and hard voxels
  - Queue label validation
- S3 method edge cases:
  - Minimal fit objects
  - Missing components

### 5. `test-validation-edge-cases.R`
Tests for input validation functions.

**Key Edge Cases Covered:**
- `.validate_fmri_data`:
  - NULL, empty, non-numeric inputs
  - Insufficient time points/voxels
  - High NA content (>50%)
  - Infinite values
  - Constant and zero voxels
  - Data frame conversion
  
- `.validate_event_model`:
  - NULL input
  - Dimension mismatches
  - No events, high/low event density
  - Object type handling
  
- `.validate_theta_bounds`:
  - Invalid structure
  - Wrong dimensions
  - Non-numeric bounds
  - Lower >= upper violations
  - Non-physiological warnings
  
- `.validate_numeric_param`:
  - NULL handling
  - Type validation
  - NA/Inf values
  - Range checking
  
- `.rock_solid_validate_inputs`:
  - Comprehensive validation
  - Default value application
  - Invalid model names
  - Parameter range violations

## Key Testing Patterns

1. **Boundary Testing**: All functions tested at parameter bounds and limits
2. **Singular/Degenerate Cases**: Zero variance, zero stimulus, constant data
3. **Numerical Precision**: Machine epsilon, mixed scales, extreme magnitudes
4. **Dimension Edge Cases**: Single voxel, single time point, high dimensions
5. **Error Recovery**: Graceful handling of invalid inputs
6. **Warning Conditions**: Non-fatal but concerning patterns

## Coverage Focus Areas

Based on algorithmic importance analysis:
- **Gauss-Newton refinement**: Critical for hard voxels
- **Ridge regression stability**: Core to all fits
- **Convolution accuracy**: Foundation of design matrix
- **R² calculation**: Key quality metric
- **Bounds enforcement**: Parameter validity

## Test Execution

All tests can be run with:
```r
testthat::test_file("tests/testthat/test-gauss-newton-edge-cases.R")
testthat::test_file("tests/testthat/test-parametric-engine-edge-cases.R")
testthat::test_file("tests/testthat/test-numerical-stability-edge-cases.R")
testthat::test_file("tests/testthat/test-integration-edge-cases.R")
testthat::test_file("tests/testthat/test-validation-edge-cases.R")
```

Or all at once:
```r
testthat::test_local(filter = "edge-cases")
```