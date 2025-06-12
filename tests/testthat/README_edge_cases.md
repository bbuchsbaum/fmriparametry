# Edge Case Tests for fmriparametric

This directory contains comprehensive edge case tests for the critical algorithmic components of the fmriparametric package.

## Test Files Created

### 1. test-gauss-newton-edge-cases.R
Tests for the Gauss-Newton nonlinear optimization algorithm:
- Singular systems (zero stimulus, near-zero predictor magnitude)
- Degenerate inputs (zero variance data)
- Parameter boundary enforcement
- Numerical precision with extreme scales
- Convergence behavior

### 2. test-parametric-engine-edge-cases.R
Tests for the core Taylor approximation engine:
- Zero stimulus handling
- Tiny amplitude division safety (epsilon_beta)
- Constant Y (zero variance) data
- Extreme parameter seeds at boundaries
- Different baseline models (with/without intercept)
- Sparse event designs

### 3. test-convolution-edge-cases.R
Tests for convolution operations:
- Empty/single-point signals
- FFT vs C++ implementation equivalence
- Energy preservation
- Boundary condition handling
- Multi-kernel and multi-signal operations

### 4. test-ridge-solver-edge-cases.R
Tests for ridge regression solver:
- Singular/rank-deficient matrices
- Extreme data scales (tiny and huge)
- Comparison with standard ridge formula
- Zero matrices handling
- Highly correlated predictors
- Single predictor edge cases

## Key Edge Cases Covered

1. **Numerical Stability**
   - Division by near-zero values
   - Extreme data scales (1e-8 to 1e6)
   - Machine precision limits

2. **Degenerate Data**
   - Zero variance signals
   - Empty/zero matrices
   - Rank-deficient designs

3. **Boundary Conditions**
   - Parameters at exact bounds
   - Parameters outside valid ranges
   - Single time point data

4. **Algorithm Equivalence**
   - FFT vs time-domain convolution
   - Different code paths produce same results

5. **Sparse Designs**
   - Single brief events
   - Very long inter-stimulus intervals
   - Minimal signal content

## Running the Tests

```r
# Run all edge case tests
testthat::test_file("tests/testthat/test-gauss-newton-edge-cases.R")
testthat::test_file("tests/testthat/test-parametric-engine-edge-cases.R")
testthat::test_file("tests/testthat/test-convolution-edge-cases.R")
testthat::test_file("tests/testthat/test-ridge-solver-edge-cases.R")

# Or run all tests including edge cases
devtools::test()
```

## Test Coverage Focus

Based on algorithmic importance ranking:
1. **parametric-engine.R** - Core estimation algorithm
2. **gauss-newton-refinement.R** - Complex optimization
3. **batch/fast_batch_convolution.R** - Numerical kernels
4. **ridge_linear_solver.R** - Linear algebra operations

These tests ensure robustness against pathological inputs that could occur in real fMRI data analysis.