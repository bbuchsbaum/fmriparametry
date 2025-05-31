# Sprint 1: Core LWU Fitting Implementation

## Sprint Overview

**Duration:** 2 weeks  
**Start Date:** [TBD]  
**End Date:** [TBD]

### Sprint Goal
Implement a minimal, functional version of `estimate_parametric_hrf` that can fit the Lag-Width-Undershoot (LWU) HRF model to voxel data using a single Taylor expansion pass. This sprint focuses on establishing the core infrastructure and basic functionality without advanced optimization features.

### Sprint Objectives
1. Set up the `fmriparametric` package structure
2. Implement the LWU HRF interface wrapper functions
3. Create data preparation utilities
4. Develop single-pass Taylor approximation engine
5. Define basic output structures and methods
6. Establish comprehensive unit testing

### Key Simplifications for Sprint 1
- **No iterative refinement:** Single Taylor expansion pass only
- **No K-means recentering:** Fixed initial seed
- **No tiered refinement:** Basic fitting only
- **Simplified error handling:** Defer standard errors to Sprint 2
- **Basic input options:** Numeric seeds and bounds only
- **Single-threaded:** No parallelization (`cores = 1`)

### Success Criteria
- [ ] Package builds and passes R CMD check
- [ ] Basic LWU parameter estimation works on test data
- [ ] Unit test coverage > 80% for implemented functions
- [ ] Core functions documented with roxygen2

---

## Sprint 1 Tickets

### Epic 1: Package Setup and HRF Interface (Est: 2 days)

#### Ticket SP1-101: Package Skeleton Setup
**Priority:** Critical  
**Estimate:** 4 hours  
**Dependencies:** None

**Description:**  
Create the foundational R package structure for `fmriparametric`.

**Acceptance Criteria:**
- [ ] Package directory structure created (R/, tests/, man/, data-raw/)
- [ ] DESCRIPTION file with proper metadata and dependencies:
  - `fmrireg (>= 0.2.0)` in Imports
  - `assertthat`, `stats`, `Matrix` in Imports
  - `numDeriv` in Suggests
- [ ] NAMESPACE file initialized
- [ ] Basic .Rbuildignore and .gitignore files
- [ ] Package loads successfully with `devtools::load_all()`

#### Ticket SP1-102: LWU HRF Interface Wrapper Functions
**Priority:** Critical  
**Estimate:** 6 hours  
**Dependencies:** SP1-101  

**Description:**  
Implement internal wrapper functions that provide the generic parametric HRF interface for the LWU model.

**Implementation Details:**
```r
# File: R/hrf-interface-lwu.R
.lwu_hrf_function <- function(t, params_vector, ...) 
.lwu_hrf_taylor_basis_function <- function(params_vector0, t_hrf_eval, ...)
.lwu_hrf_parameter_names <- function()
.lwu_hrf_default_seed <- function()
.lwu_hrf_default_bounds <- function()
```

**Acceptance Criteria:**
- [ ] `.lwu_hrf_function()` correctly calls `fmrireg::hrf_lwu` with proper parameter mapping
- [ ] `.lwu_hrf_taylor_basis_function()` returns matrix with correct dimensions (length(t) × 4)
- [ ] Parameter names are `c("tau", "sigma", "rho")`
- [ ] Default seed is `c(6, 2.5, 0.35)`
- [ ] Default bounds are physiologically reasonable
- [ ] All functions have unit tests with > 90% coverage

### Epic 2: Main Function and Data Preparation (Est: 3 days)

#### Ticket SP1-103: Main Function Skeleton & Input Validation
**Priority:** Critical  
**Estimate:** 8 hours  
**Dependencies:** SP1-102  

**Description:**  
Implement the main user-facing function with input parsing and validation.

**Function Signature:**
```r
estimate_parametric_hrf <- function(
  fmri_data,
  event_model,
  parametric_hrf = "lwu",
  theta_seed = NULL,
  theta_bounds = NULL,
  confound_formula = NULL,
  baseline_model = "intercept",
  hrf_eval_times = NULL,
  hrf_span = 30,
  lambda_ridge = 0.01,
  mask = NULL,
  verbose = FALSE
)
```

**Acceptance Criteria:**
- [ ] Function validates all inputs using `assertthat`
- [ ] Resolves `parametric_hrf = "lwu"` to interface functions
- [ ] Handles NULL theta_seed by using defaults
- [ ] Validates theta_bounds are within physiological ranges
- [ ] Clear error messages for invalid inputs
- [ ] Unit tests for all validation paths

#### Ticket SP1-104: Data Preparation Implementation
**Priority:** Critical  
**Estimate:** 12 hours  
**Dependencies:** SP1-103  

**Description:**  
Implement `.prepare_parametric_inputs()` to prepare BOLD data and event designs for fitting.

**Function Signature:**
```r
.prepare_parametric_inputs <- function(
  fmri_data,
  event_model,
  confound_formula = NULL,
  baseline_model = "intercept",
  hrf_eval_times = NULL,
  hrf_span = 30,
  mask = NULL
)
```

**Implementation Steps:**
1. Extract BOLD time series matrix `Y_raw` and `scan_times`
2. Generate stimulus regressor `S_target` from event model
3. Construct confound matrix `Z` if specified
4. Project out confounds using QR decomposition
5. Generate HRF evaluation time points

**Acceptance Criteria:**
- [ ] Correctly extracts data from `fmri_dataset` objects
- [ ] Handles both sparse and dense event representations
- [ ] Confound projection preserves data structure
- [ ] Returns list with all required components
- [ ] Handles edge cases (no events, all zero voxels)
- [ ] Memory efficient for large datasets

### Epic 3: Core Fitting Engine (Est: 3 days)

#### Ticket SP1-105: Parametric Engine - Single Pass Implementation
**Priority:** Critical  
**Estimate:** 16 hours  
**Dependencies:** SP1-104  

**Description:**  
Implement the core Taylor approximation engine for single-pass parameter estimation.

**Function Signature:**
```r
.parametric_engine <- function(
  Y_proj,
  S_target_proj,
  scan_times,
  hrf_eval_times,
  hrf_interface,
  theta_seed,
  theta_bounds,
  lambda_ridge = 0.01,
  epsilon_beta = 1e-6
)
```

**Algorithm Steps:**

1. **Taylor Basis Construction**
   ```r
   X_taylor <- hrf_interface$taylor_basis(theta_seed, hrf_eval_times)
   # Dimensions: length(hrf_eval_times) × 4
   ```

2. **Design Matrix Construction**
   - Convolve S_target with each Taylor basis column
   - Use efficient convolution via fmrireg utilities
   - Result: X_design (N_timepoints × 4)

3. **Global QR Decomposition**
   ```r
   qr_decomp <- qr(X_design)
   Q <- qr.Q(qr_decomp)
   R <- qr.R(qr_decomp)
   R_inv <- solve(R + lambda_ridge * diag(ncol(R)))
   ```

4. **Voxel-wise Estimation**
   ```r
   # Vectorized computation
   coeffs <- R_inv %*% t(Q) %*% Y_proj
   beta0 <- coeffs[1, ]
   delta_theta <- coeffs[2:4, ] / beta0
   theta_hat <- theta_seed + t(delta_theta)
   ```

5. **Parameter Bounds**
   ```r
   theta_hat <- pmax(theta_bounds$lower, pmin(theta_hat, theta_bounds$upper))
   ```

**Acceptance Criteria:**
- [ ] Correctly implements Taylor approximation mathematics
- [ ] Handles rank-deficient design matrices gracefully
- [ ] Vectorized operations for efficiency
- [ ] Proper handling of near-zero amplitudes
- [ ] Parameter bounds correctly applied
- [ ] Returns structured output with parameters and amplitudes

### Epic 4: Output Structure and Methods (Est: 1.5 days)

#### Ticket SP1-106: S3 Class Definition
**Priority:** High  
**Estimate:** 4 hours  
**Dependencies:** SP1-105  

**Description:**  
Define the `parametric_hrf_fit` S3 class structure.

**Class Structure:**
```r
structure(
  list(
    estimated_parameters = matrix(),  # V × 3 for LWU
    amplitudes = numeric(),          # V × 1
    parameter_names = character(),   # c("tau", "sigma", "rho")
    hrf_model = "lwu",
    convergence = list(),           # Placeholder for Sprint 2
    metadata = list(
      call = match.call(),
      n_voxels = integer(),
      n_timepoints = integer(),
      theta_seed = numeric(),
      theta_bounds = list()
    )
  ),
  class = "parametric_hrf_fit"
)
```

**Acceptance Criteria:**
- [ ] Constructor function `new_parametric_hrf_fit()` implemented
- [ ] Validates all required fields
- [ ] Consistent structure documented
- [ ] Helper functions for field access

#### Ticket SP1-107: Integration in Main Function
**Priority:** High  
**Estimate:** 4 hours  
**Dependencies:** SP1-103, SP1-104, SP1-105, SP1-106  

**Description:**  
Integrate all components in the main `estimate_parametric_hrf()` function.

**Implementation Flow:**
```r
estimate_parametric_hrf <- function(...) {
  # 1. Input validation (SP1-103)
  # 2. Prepare data (SP1-104)
  prepared <- .prepare_parametric_inputs(...)
  
  # 3. Run engine (SP1-105)
  fit_results <- .parametric_engine(...)
  
  # 4. Create output (SP1-106)
  new_parametric_hrf_fit(
    estimated_parameters = fit_results$theta_hat,
    amplitudes = fit_results$beta0,
    ...
  )
}
```

**Acceptance Criteria:**
- [ ] Main function integrates all components
- [ ] Error handling at each step
- [ ] Informative progress messages if verbose=TRUE
- [ ] Returns valid parametric_hrf_fit object

#### Ticket SP1-108: Basic S3 Methods
**Priority:** Medium  
**Estimate:** 4 hours  
**Dependencies:** SP1-106  

**Description:**  
Implement basic S3 methods for the parametric_hrf_fit class.

**Methods to Implement:**
```r
# Print method
print.parametric_hrf_fit <- function(x, ...) {
  cat("Parametric HRF Fit\n")
  cat("Model:", x$hrf_model, "\n")
  cat("Voxels:", nrow(x$estimated_parameters), "\n")
  cat("\nParameter Summary:\n")
  print(summary(x$estimated_parameters))
}

# Basic coef method
coef.parametric_hrf_fit <- function(object, ...) {
  object$estimated_parameters
}

# Basic summary method  
summary.parametric_hrf_fit <- function(object, ...) {
  # Return summary statistics
}
```

**Acceptance Criteria:**
- [ ] Print method shows key information clearly
- [ ] Coef method returns parameter matrix
- [ ] Summary provides useful statistics
- [ ] All methods handle edge cases

### Epic 5: Testing and Documentation (Est: 2.5 days)

#### Ticket SP1-109: Unit Tests - Data Preparation
**Priority:** High  
**Estimate:** 6 hours  
**Dependencies:** SP1-104  

**Description:**  
Comprehensive unit tests for data preparation functions.

**Test Cases:**
```r
test_that(".prepare_parametric_inputs handles basic case", {
  # Setup test data
  # Test extraction, projection, timing
})

test_that("confound projection works correctly", {
  # Test with and without confounds
  # Verify orthogonalization
})

test_that("event design construction is correct", {
  # Test sparse and dense events
  # Verify summing of factorial terms
})
```

**Acceptance Criteria:**
- [ ] Tests cover all major code paths
- [ ] Edge cases tested (empty data, no events)
- [ ] Performance benchmarks included
- [ ] Tests pass consistently
#### Ticket SP1-110: Unit Tests - Parametric Engine
**Priority:** High  
**Estimate:** 8 hours  
**Dependencies:** SP1-105  

**Description:**  
Unit tests for the core fitting engine with synthetic data.

**Test Scenarios:**

1. **Perfect Recovery Test**
   ```r
   # Generate data with known LWU parameters
   true_params <- c(tau = 6, sigma = 2.5, rho = 0.35)
   # Seed at true values → should recover exactly
   ```

2. **Convergence Test**
   ```r
   # Start from offset seed
   # Verify movement toward true parameters
   ```

3. **Bounds Test**
   ```r
   # Test parameter clamping
   # Verify bounds are respected
   ```

4. **Numerical Stability**
   ```r
   # Test with poorly conditioned data
   # Verify no crashes or NaN values
   ```

**Acceptance Criteria:**
- [ ] Recovery error < 1% for perfect case
- [ ] Convergence improves parameter estimates
- [ ] Bounds always respected
- [ ] Handles edge cases gracefully
#### Ticket SP1-111: Integration Tests
**Priority:** High  
**Estimate:** 6 hours  
**Dependencies:** SP1-107  

**Description:**  
End-to-end integration tests for the complete workflow.

**Test Cases:**

1. **Basic Workflow Test**
   ```r
   test_that("full estimation workflow works", {
     # Create simple test dataset
     data <- fmrireg::matrix_dataset(...)
     events <- fmrireg::event_model(...)
     
     # Run estimation
     fit <- estimate_parametric_hrf(data, events)
     
     # Verify output structure
     expect_s3_class(fit, "parametric_hrf_fit")
     expect_equal(ncol(fit$estimated_parameters), 3)
   })
   ```

2. **Realistic Data Test**
   ```r
   # Use more complex simulated data
   # Multiple voxels with varying SNR
   # Verify reasonable parameter ranges
   ```

**Acceptance Criteria:**
- [ ] Complete workflow executes without errors
- [ ] Output structure is correct
- [ ] Parameter estimates are physiologically plausible
- [ ] Performance is acceptable (< 10s for 1000 voxels)

#### Ticket SP1-112: Documentation
**Priority:** Medium  
**Estimate:** 4 hours  
**Dependencies:** All implementation tickets  

**Description:**  
Create comprehensive documentation for all exported functions.

**Documentation Requirements:**
- [ ] Roxygen2 documentation for all exported functions
- [ ] Internal function documentation
- [ ] README.md with installation and basic usage
- [ ] NEWS.md file started
- [ ] Basic vignette outline

---

## Sprint 1 Deliverables

### Code Deliverables
1. **Package Structure**
   - [x] `fmriparametric` R package skeleton
   - [x] Proper dependency management
   - [x] Build passes R CMD check

2. **Core Functionality**
   - [x] LWU HRF interface wrapper functions
   - [x] Data preparation pipeline
   - [x] Single-pass Taylor approximation engine
   - [x] Main `estimate_parametric_hrf()` function

3. **Output & Methods**
   - [x] `parametric_hrf_fit` S3 class
   - [x] Basic print, coef, and summary methods

4. **Testing**
   - [x] Unit tests with >80% coverage
   - [x] Integration tests for full workflow
   - [x] Performance benchmarks

5. **Documentation**
   - [x] Function documentation via roxygen2
   - [x] Basic package documentation
   - [x] README with examples

### Technical Achievements
- Single-pass LWU parameter estimation functional
- Efficient vectorized implementation
- Robust handling of edge cases
- Clean, extensible architecture

### Known Limitations (To Address in Sprint 2)
- No iterative refinement
- No standard error calculation
- Single-threaded execution only
- Limited to numeric seed specification
- No convergence diagnostics

### Sprint Retrospective Topics
1. What went well?
2. What challenges were encountered?
3. Are we on track for Sprint 2 goals?
4. Any architectural decisions to revisit?

---

## Transition to Sprint 2

Sprint 2 will build upon this foundation by adding:
- Iterative global recentering
- K-means recentering for robustness  
- Tiered refinement for difficult voxels
- Complete standard error calculation
- R² and residual computation
- Parallel processing support
- Enhanced S3 methods (plot, predict)
- Performance optimizations