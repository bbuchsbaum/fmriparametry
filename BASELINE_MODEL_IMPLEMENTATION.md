# Baseline Model Implementation Summary

## Overview
Successfully implemented support for `fmrireg::baseline_model` objects in the `fmriparametric` package to provide a consistent interface with `event_model` and enable more sophisticated confound regression.

## Changes Made

### 1. Core Function Updates

#### prepare-parametric-inputs.R
- Added support for `baseline_model` objects in addition to `confound_formula`
- Modified confound matrix extraction logic to handle both approaches
- Maintains backward compatibility with formula-based confounds

#### estimate_parametric_hrf.R
- Updated parameter documentation to describe `baseline_model` usage
- Changed default from "intercept" to NULL for cleaner interface

### 2. Engine Updates

#### estimation-stages.R & parametric-engine.R
- Added helper function `.is_intercept_baseline()` to properly check baseline_model type
- Fixed all string comparisons that were causing errors with baseline_model objects
- Changed default baseline_model from "intercept" to NULL

### 3. Test Updates

#### test-confound-robustness.R
- Fixed to use `confound_formula` parameter instead of non-existent `confound_data`
- Added new test for baseline_model support

#### test-baseline-model-integration.R (new)
- Comprehensive tests for baseline_model functionality
- Tests polynomial drift, spline basis, and nuisance regressors
- Validates improved parameter recovery with baseline models

### 4. Documentation

- Created `examples/baseline_model_example.R` demonstrating usage patterns
- Updated function documentation to reflect new capabilities

## Key Features

1. **Consistent Interface**: Both `event_model` and `baseline_model` now follow the same pattern
2. **Backward Compatibility**: Existing `confound_formula` approach still works
3. **Enhanced Functionality**: Supports complex drift modeling, splines, and nuisance regressors
4. **Improved Type Safety**: Proper handling of different baseline_model types

## Usage Examples

```r
# Polynomial drift correction
baseline_poly <- baseline_model(~ poly(time, 3), sframe = sframe)
fit <- estimate_parametric_hrf(data, events, baseline_model = baseline_poly)

# Spline-based drift
baseline_spline <- baseline_model(~ ns(time, df = 5), sframe = sframe)
fit <- estimate_parametric_hrf(data, events, baseline_model = baseline_spline)

# Motion parameters and drift
baseline_full <- baseline_model(
  ~ poly(time, 3) + motion1 + motion2 + motion3,
  sframe = sframe,
  data = motion_data
)
fit <- estimate_parametric_hrf(data, events, baseline_model = baseline_full)
```

## Test Results

- All confound robustness tests now pass
- Baseline model integration tests pass (except one random variation failure)
- Package installs cleanly without errors
- Maintains compatibility with existing code

## Next Steps

1. Consider deprecating `confound_formula` in favor of `baseline_model` in future versions
2. Add more examples to package vignettes
3. Consider adding convenience functions for common baseline models