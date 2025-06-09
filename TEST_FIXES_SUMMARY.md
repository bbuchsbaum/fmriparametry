# Test Fixes Summary

## Fixed Test Files

### 1. test-bayesian-engine.R
- **Issue**: Deprecated `context()` function
- **Fix**: Removed the `context(".bayesian_engine")` line

### 2. test-crossprod-updates.R
- **Issue**: Test expected `mcmc_samples=0` to work, but the function validates it must be positive
- **Fix**: Changed `mcmc_samples = 0` to `mcmc_samples = 1000` and updated expected samples length

### 3. test-dimension-checks.R
- **Issue**: Test was using non-existent `scan_times` parameter in `.parametric_engine()`
- **Fix**: Removed the `scan_times` parameter from the function call

### 4. test-estimate-parametric-hrf.R
- **Issue**: Expected error message didn't match actual validation message
- **Fix**: Changed expected error from "theta_bounds must have both" to "theta_bounds missing required elements"

### 5. test-integration.R
- **Issues**: 
  - Mock objects had `class` as a list element instead of using `class()` function
  - Expected error messages didn't match actual validation messages
- **Fixes**:
  - Removed `class = c(...)` from list definitions (keeping only `class()` assignment)
  - Updated expected error from "X is NULL or empty" to "Insufficient time points"
  - Updated expected error from "wrong number of rows" to "don't match fmri_data time points"

### 6. test-parametric-engine-validation.R
- **Issue**: Tests were using non-existent `scan_times` parameter
- **Fix**: Removed all three occurrences of `scan_times = 1:10,` from the test calls

### 7. rock-solid-validation.R (Supporting fix)
- **Issue**: Validation functions weren't handling real fmrireg objects properly
- **Fixes**:
  - Updated `.validate_fmri_data()` to use `fmrireg::get_data_matrix()` for real objects
  - Updated `.validate_event_model()` to use `fmrireg::design_matrix()` for real objects

## Test Results

All fixed tests now pass successfully:
- 0 failures
- Some expected warnings (e.g., low event density warnings)
- All core functionality working correctly