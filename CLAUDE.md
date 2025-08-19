# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`fmriparametric` is an R package for parametric Hemodynamic Response Function (HRF) estimation from fMRI data using iterative linear Taylor approximation. The package focuses on the Lag-Width-Undershoot (LWU) model with 3 parameters: lag (τ), width (σ), and undershoot (ρ).

## Common Development Commands

### Testing
```bash
# Run all tests
R -e "testthat::test_local()"

# Run specific test file
R -e "testthat::test_file('tests/testthat/test-estimate-parametric-hrf.R')"

# Run single test
R -e "testthat::test_file('tests/testthat/test-estimate-parametric-hrf.R', reporter = 'summary')" --args "test_name"

# Run tests with coverage
R -e "covr::package_coverage()"
```

### Package Development
```bash
# Build and check package
R CMD build .
R CMD check fmriparametric_*.tar.gz

# Document and install
R -e "devtools::document()"
R -e "devtools::install()"

# Full check with devtools
R -e "devtools::check()"

# Build vignettes
R -e "devtools::build_vignettes()"
```

### Installation
```r
# Install dependencies from GitHub
devtools::install_github("bbuchsbaum/fmrireg")
devtools::install_github("bbuchsbaum/fmrihrf") 
devtools::install_github("bbuchsbaum/fmridataset")

# Install this package
devtools::install_github("bbuchsbaum/fmriparametry")
```

## Architecture

### Core Workflow

The package implements a staged HRF estimation workflow:

1. **Stage 0**: Input validation and data preparation (`.stage0_validate_and_prepare`)
2. **Stage 1**: Core estimation with Taylor approximation (`.stage1_initial_estimation`)
3. **Stage 2**: K-means initialization (optional) (`.stage2_kmeans_initialization`)
4. **Stage 3**: Global refinement iterations (`.stage3_global_refinement`)
5. **Stage 4**: Tiered refinement for difficult voxels (`.stage4_tiered_refinement`)

### Key Components

**Main Entry Points:**
- `estimate_parametric_hrf()` - Primary user-facing function with Strangler Fig pattern for legacy support
- `.run_hrf_estimation_engine()` - Modular refactored implementation (default)

**HRF Interfaces:**
- `.create_hrf_interface()` - Factory for HRF model interfaces
- `.lwu_hrf_function()` - LWU HRF evaluation
- `.lwu_hrf_taylor_basis_function()` - Taylor basis with derivatives

**Core Algorithms:**
- `.bayesian_engine()` - Main fitting engine with Taylor approximation
- `.gauss_newton_refinement()` - Iterative refinement for difficult voxels  
- `.local_recentering_moderate()` - Moderate refinement strategy
- `.ridge_linear_solve()` - Ridge-regularized linear solver

**Data Preparation:**
- Input validation with tiered approach (minimal/standard/comprehensive)
- Confound regression via `baseline_model` parameter
- Design matrix construction from event models

**Output:**
- `parametric_hrf_fit` S3 class with methods: print, summary, coef, fitted, residuals, predict, plot
- Diagnostic reports: timing, memory usage, convergence tracking

### Parameter Bounds

LWU model default bounds:
- Lag (τ): [0.5, 10] seconds
- Width (σ): [0.5, 5] seconds  
- Undershoot (ρ): [0, 1.5]

Bounds are enforced throughout optimization with stability adjustments when parameters approach boundaries.

### Performance Optimizations

- Vectorized operations using Matrix package for sparse matrices
- Optional parallel processing for voxel-wise operations
- Fast batch convolution via FFT (`.fast_batch_convolution`)
- Caching of validation results and HRF evaluations
- Convergence-based early stopping

### Testing Infrastructure

Tests organized by functionality:
- Edge cases documented in `tests/testthat/README_edge_cases.md`
- Comprehensive test coverage for bounds enforcement, convergence, numerical stability
- Integration tests with real and simulated fMRI data
- Characterization tests for algorithm behavior

## Important Implementation Notes

1. **Package Name**: Directory is `fmriparametry` but package name is `fmriparametric`

2. **Refinement Strategies**: 
   - "none": Basic Taylor approximation only
   - "moderate": Local recentering for poor fits
   - "aggressive": Gauss-Newton for challenging voxels

3. **Safety Modes**:
   - "performance": Minimal validation
   - "balanced": Standard validation (default)
   - "maximum": Comprehensive validation with diagnostics

4. **Ridge Regularization**: Default λ = 0.01, adjustable for numerical stability

5. **Convergence Criteria**: 
   - Global refinement: max parameter change < epsilon (default 0.01)
   - Gauss-Newton: relative change in objective < 1e-6 or max iterations

6. **Dependencies**: Requires GitHub packages fmrireg, fmrihrf, fmridataset from bbuchsbaum