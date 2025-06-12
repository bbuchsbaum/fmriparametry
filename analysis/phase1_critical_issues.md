# Phase 1 Audit: Critical Issues Found

## Overview

Phase 1 audit of the critical path files has revealed several significant issues that impact correctness, performance, and maintainability.

## Critical Issues (Must Fix)

### 1. **Mathematical Error in Gauss-Newton Jacobian** ⚠️
- **File**: `gauss-newton-refinement.R`, lines 304-305
- **Issue**: Missing factor of 2 in derivative calculation
- **Impact**: Incorrect search direction, slower convergence, suboptimal results
- **Fix**: 
  ```r
  # Correct calculation
  dbeta_dtheta_k <- (crossprod(dx_dtheta_k, y) - 2 * beta * crossprod(x_hrf, dx_dtheta_k)) / denom
  ```

### 2. **Inconsistent R² Calculation**
- **Files**: Multiple locations with different implementations
- **Issue**: Post-refinement R² doesn't include intercept adjustment
- **Impact**: Incorrect model fit assessment after refinement
- **Fix**: Create single canonical R² function and ensure intercept is included

### 3. **Monolithic Function (1098 lines)**
- **File**: `estimate_parametric_hrf.R`
- **Issue**: Single function handles validation, initialization, estimation, refinement, and results
- **Impact**: Unmaintainable, untestable, difficult to debug
- **Fix**: Refactor into separate stage functions

## High Priority Performance Issues

### 4. **Non-vectorized Post-Refinement R² Calculation**
- **File**: `estimate_parametric_hrf.R`, lines 607-632
- **Issue**: Nested loops over voxels and stimulus columns
- **Impact**: Major performance bottleneck for large datasets
- **Fix**: Use batch convolution and matrix operations

### 5. **Inefficient Batch Convolution**
- **File**: `batch_convolution.R`
- **Issue**: R loop over signals instead of vectorization
- **Impact**: Unnecessary overhead, memory churn
- **Fix**: Sum signals first, then single convolution

### 6. **Sequential Voxel Processing in Refinement**
- **Files**: `gauss-newton-refinement.R`, `local-recentering.R`
- **Issue**: For loops over voxels
- **Impact**: Cannot utilize multiple cores
- **Fix**: Parallelize using foreach or future.apply

## Numerical Stability Concerns

### 7. **Ridge Solver Edge Case**
- **File**: `ridge_linear_solver.R`
- **Issue**: No check for singular matrices when lambda=0
- **Impact**: Potential crash or garbage results with collinear design
- **Recommendation**: Document behavior or add pre-check

### 8. **Matrix Dimension Dropping**
- **Multiple files**
- **Issue**: Reactive checks for matrix dimensions throughout code
- **Impact**: Fragile code, unexpected vector conversions
- **Fix**: Use `drop = FALSE` consistently

## Next Steps

1. **Immediate**: Fix the Jacobian calculation error (Issue #1)
2. **Short term**: Create unified R² function and fix post-refinement calculation
3. **Medium term**: Refactor estimate_parametric_hrf.R into modular functions
4. **Performance**: Vectorize computations and add parallelization

These issues represent fundamental problems that affect the package's correctness and usability. The Jacobian error in particular must be fixed immediately as it compromises the mathematical validity of the Gauss-Newton refinement.