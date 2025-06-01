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
```

### Linting and Type Checking
No specific linting or type checking commands are currently configured. Consider adding lintr or other R linting tools.

## Architecture

### Core Components

1. **Main Entry Point**: `estimate_parametric_hrf()` in `R/estimate_parametric_hrf.R`
   - User-facing function for HRF parameter estimation
   - Handles initialization, iteration control, and output packaging

2. **Data Preparation**: `.prepare_parametric_inputs()` in `R/prepare-parametric-inputs.R`
   - Extracts BOLD data from fmri_dataset or matrix_dataset
   - Constructs design matrices from event models
   - Handles confound regression
   - Creates timing vectors for HRF evaluation

3. **Fitting Engine**: `.parametric_engine()` in `R/parametric-engine.R`
   - Implements iterative Taylor approximation algorithm
   - Single-pass estimation with ridge regularization
   - Enforces parameter bounds
   - Returns parameter estimates with standard errors

4. **HRF Interface**: Functions in `R/hrf-interface-lwu.R`
   - `.hrf_lwu_wrapper()`: Evaluates LWU HRF at given parameters
   - `.hrf_taylor_basis_lwu_wrapper()`: Constructs Taylor basis (HRF + derivatives)
   - `.get_lwu_defaults()`: Provides default seeds and bounds

5. **Output Class**: `parametric_hrf_fit` S3 class
   - Defined in `R/parametric-hrf-fit-class.R`
   - Methods in `R/parametric-hrf-fit-methods.R`: print, coef, summary

### Key Design Patterns

- **Extensibility**: Architecture supports adding new HRF models beyond LWU
- **Separation of Concerns**: Clear separation between data prep, fitting engine, and model-specific code
- **Defensive Programming**: Extensive input validation with assertthat
- **Performance**: Vectorized operations using Matrix package for sparse matrices

### External Dependencies

- **fmrireg**: Core fMRI functionality (HRF models, event modeling, data structures)
- **Matrix**: Sparse matrix operations for efficiency
- **RcppEigen**: Fast linear algebra operations
- **assertthat**: Input validation

### Testing Strategy

Tests are organized by component:
- `test-estimate-parametric-hrf.R`: Main function integration tests
- `test-prepare-inputs.R`: Data preparation logic
- `test-parametric-engine.R`: Core fitting algorithm
- `test-lwu-interface.R`: HRF model interface
- `test-parametric-hrf-fit-class.R`: Output class construction
- `test-parametric-hrf-fit-methods.R`: S3 methods

## Important Notes

1. **Package Name Mismatch**: Directory is `fmriparametry` but package name in DESCRIPTION is `fmriparametric`

2. **Development Status**: Currently implements single-pass Taylor approximation. Advanced features from proposal (adaptive seeding, K-means recentering, tiered refinement) are not yet implemented.

3. **Parameter Bounds**: LWU model enforces σ > 0.05 and 0 ≤ ρ ≤ 1.5

4. **Ridge Regularization**: Default lambda = 0.01 for numerical stability