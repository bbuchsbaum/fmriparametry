# R-Squared Bug Fix Summary

## Issue Identified

The parametric HRF estimation was producing R² = 0 for all voxels due to a missing intercept term in the design matrix.

### Root Cause

1. fMRI data typically has a large baseline offset (e.g., mean ~100-1000)
2. The HRF Taylor basis only contained the HRF and its derivatives (4 columns)
3. No intercept/baseline term was included in the regression
4. This caused:
   - Fitted values to have a very different mean than the data
   - SS_residual >> SS_total
   - Negative R² values that were clamped to 0

## Solution Implemented

Working with Gemini, we added proper baseline handling to the parametric engine:

1. **Added `baseline_model` parameter** to `.parametric_engine()` in `R/parametric-engine.R`
   - Default value: `"intercept"`
   - When set, prepends a column of 1s to the design matrix

2. **Updated coefficient extraction logic**:
   - With intercept: coeffs[1,] = intercept, coeffs[2,] = HRF amplitude
   - Without intercept: coeffs[1,] = HRF amplitude (original behavior)

3. **Propagated the parameter** through the call chain:
   - `estimate_parametric_hrf()` passes `baseline_model` to the engine
   - Helper functions updated to include the parameter

## Test Results

Before fix:
- All R² values = 0
- Fitted values mean ~10 vs data mean ~100

After fix:
- R² values now positive (0.33 - 0.51 in test)
- Model correctly captures baseline
- Proper variance explained calculation

## Remaining Issues

1. The refinement stages have a separate bug causing dimension errors
2. The `baseline_model` parameter in `prepare_parametric_inputs()` is still marked as "unused"
3. Consider whether baseline should be projected out vs modeled

## Code Changes

Key changes in `parametric-engine.R`:
```r
# Add intercept if requested
has_intercept <- FALSE
if (!is.null(baseline_model) && baseline_model == "intercept") {
  intercept_col <- matrix(1, nrow = n_time, ncol = 1)
  X_design <- cbind(intercept_col, X_design)
  has_intercept <- TRUE
}
```

This fix resolves the core R² = 0 issue, allowing the package to produce meaningful goodness-of-fit metrics.