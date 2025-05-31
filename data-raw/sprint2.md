
# Sprint 2: Iterative Refinement, Diagnostics, and Performance

## Sprint Overview

**Duration:** 2 weeks  
**Start Date:** [TBD]  
**End Date:** [TBD]  
**Prerequisites:** Successful completion of Sprint 1

### Sprint Goal
Enhance `estimate_parametric_hrf` with iterative global re-centering of the Taylor expansion point θ₀, comprehensive diagnostics (R², residuals), complete standard error estimation via Delta method, and performance optimizations. Transform the single-pass implementation into a robust, production-ready tool.

### Sprint Objectives
1. Implement iterative global re-centering of Taylor expansion point
2. Add data-driven parameter initialization ("data_median")
3. Calculate R², residuals, and model diagnostics
4. Implement complete standard error estimation
5. Add parallel processing support
6. Expand S3 methods (summary, coef, fitted, residuals, plot)
7. Comprehensive testing and documentation

### Building on Sprint 1
- Assumes working `.prepare_parametric_inputs` helper
- Assumes functional single-pass logic in `.parametric_engine`
- Extends basic `parametric_hrf_fit` class with diagnostics

### Success Criteria
- [ ] Iterative refinement improves parameter estimates by >30%
- [ ] R² calculation validates model fit quality
- [ ] Standard errors match theoretical expectations
- [ ] All S3 methods produce publication-ready output
- [ ] Test coverage remains >85%

---

## Sprint 2 Tickets

### Epic 1: Iterative Global Re-Centering (Est: 3 days)

#### Ticket SP2-201: Implement Global θ₀ Re-Centering Loop
**Priority:** Critical  
**Estimate:** 16 hours  
**Dependencies:** SP1-105 (parametric engine)

**Description:**  
Modify `.parametric_engine` to support iterative refinement through global re-centering of the Taylor expansion point.

**Implementation Details:**
Wrap the single-pass Taylor approximation logic (from Sprint 1) in an iterative loop controlled by `recenter_global_passes` and `recenter_epsilon`.

**Algorithm Steps:**

1. **Initialize tracking variables:**
   - `theta_current_global ← theta_seed`
   - `R2_voxel ← vector(-Inf, n_voxels)`
   - `theta_trajectory ← list()`

2. **For each pass (pass_global = 1 to recenter_global_passes):**

   a. **Perform single Taylor pass** using `theta_current_global` as expansion point:
      - Execute steps from SP1-105 (Taylor basis, design matrix, QR solve)
      - Obtain `hat_theta_v_pass` and `hat_beta0_v_pass` for all voxels

   b. **Calculate R² for current pass:**
      ```r
      for (v in 1:n_voxels) {
        Fitted_v_pass <- X_design %*% current_pass_coeffs[, v]
        Residuals_v_pass <- Y_proj[, v] - Fitted_v_pass
        SS_res <- sum(Residuals_v_pass^2)
        SS_tot <- sum((Y_proj[, v] - mean(Y_proj[, v]))^2)
        R2_v_pass <- 1 - SS_res / SS_tot
      }
      ```

   c. **Update best voxel estimates** if R² improved:
      ```r
      if (R2_v_pass > R2_voxel[v]) {
        theta_hat_voxel[, v] <- hat_theta_v_pass[, v]
        beta0_voxel[v] <- hat_beta0_v_pass[v]
        R2_voxel[v] <- R2_v_pass
        coeffs_voxel[, v] <- current_pass_coeffs[, v]
      }
      ```

   d. **Re-center global θ₀** (if not final pass):
      - Select good voxels: `idx_good <- which(R2_voxel >= r2_threshold)`
      - If good voxels exist:
        - Compute robust median: `theta_new <- apply(theta_hat[, idx_good], 1, median)`
        - Apply bounds: `theta_new <- clamp(theta_new, bounds)`
        - Check convergence: if `max(abs(theta_new - theta_current)) < epsilon`, break
        - Update: `theta_current_global <- theta_new`
      - Else: warn and break

3. **Store convergence information:**
   - Trajectory of theta values
   - Number of iterations completed
   - Final global theta

**Acceptance Criteria:**
- [ ] Re-centering improves average R² by >20%
- [ ] Convergence detection works correctly
- [ ] Handles edge cases (no good voxels)
- [ ] Trajectory stored for diagnostics
- [ ] Performance overhead <2x single pass

#### Ticket SP2-202: Data-Driven Initialization ("data_median")
**Priority:** High  
**Estimate:** 8 hours  
**Dependencies:** SP2-201

**Description:**  
Implement data-driven parameter seeding when user specifies `theta_seed = "data_median"`.

**Algorithm:**

1. **Preliminary pass with default seed:**
   ```r
   if (theta_seed == "data_median") {
     # Use hardcoded default as initial seed
     default_theta0 <- .lwu_hrf_default_seed()  # e.g., c(6, 2.5, 0.35)
     
     # Single pass to get initial estimates
     prelim_results <- .parametric_engine(
       ...,
       theta_seed = default_theta0,
       recenter_global_passes = 0
     )
   }
   ```

2. **Select good voxels based on R²:**
   ```r
   # Use upper quartile as threshold
   r2_threshold <- quantile(prelim_results$R2_voxel, 0.75, na.rm = TRUE)
   idx_good <- which(prelim_results$R2_voxel >= r2_threshold)
   ```

3. **Compute robust median of good voxels:**
   ```r
   if (length(idx_good) >= 10) {
     theta_data_median <- apply(
       prelim_results$theta_hat[idx_good, , drop = FALSE],
       2,
       median,
       na.rm = TRUE
     )
   } else {
     warning("Too few good voxels; using default seed")
     theta_data_median <- default_theta0
   }
   ```

4. **Main estimation with data-driven seed:**
   ```r
   # Now run full estimation with data-driven seed
   main_results <- .parametric_engine(
     ...,
     theta_seed = theta_data_median,
     recenter_global_passes = user_specified_passes
   )
   ```

**Acceptance Criteria:**
- [ ] Correctly identifies good voxels from preliminary pass
- [ ] Robust median calculation handles outliers
- [ ] Fallback to default when insufficient good voxels
- [ ] Data-driven seed improves final results vs default
- [ ] Clear user messaging about seed selection

### Epic 2: Diagnostics and Standard Errors (Est: 2.5 days)

#### Ticket SP2-203: Compute Final Residuals
**Priority:** Critical  
**Estimate:** 6 hours  
**Dependencies:** SP2-201

**Description:**  
Calculate final residuals after all re-centering iterations complete.

**Implementation:**

After the re-centering loop, for each voxel v:

1. **Retrieve best estimates:**
   - Parameters: `theta_v = theta_hat_voxel[, v]`
   - Coefficients: `coeffs_v = coeffs_voxel[, v]`

2. **Reconstruct design matrix:**
   - For Sprint 2, use final global expansion point
   - `X_design_final = construct_design_matrix(theta_current_global_final)`
   - Note: Sprint 3 will handle voxel-specific expansion points

3. **Calculate fitted values and residuals:**
   ```r
   Fitted_v_final <- X_design_final %*% coeffs_v
   residuals_matrix[, v] <- Y_proj[, v] - Fitted_v_final
   ```

**Output Enhancement:**
- Add `residuals_matrix` (N × V) to `.parametric_engine` output
- Store for use in S3 methods and diagnostics

**Acceptance Criteria:**
- [ ] Residuals sum to approximately zero per voxel
- [ ] No systematic patterns in residuals
- [ ] Residuals uncorrelated with fitted values
- [ ] Handles missing data appropriately
- [ ] Memory efficient for large datasets

*   **Ticket SP2-204: Implement Standard Error Calculation in `.parametric_engine`**
    *   This implements step 3.1.2.f from the full proposal.
    *   **Inputs:** Final `theta_hat_voxel`, `coeffs_voxel`, `Y_proj`, `S_target_proj`, `scan_times_bold`, `hrf_eval_times`, `lambda_ridge_jacobian`, `P`.
    *   **For each voxel `v` (can be parallelized):**
        1.  Final estimates: \(\hat\theta_v = \text{theta\_hat_voxel}[,v]\). The \(\hat{\text{coeffs}}_v = \text{coeffs\_voxel}[,v]\) are the linear coefficients obtained when expanding around some \(\theta_{exp,v}\) (for Sprint 2, this is `theta_current_global_final`).
        2.  Reconstruct the design matrix \(X_{design}(\theta_{exp,v})\) that led to \(\hat{\text{coeffs}}_v\).
        3.  Residuals: \(e_v = Y_{proj,v} - X_{design}(\theta_{exp,v}) \hat{\text{coeffs}}_v\).
        4.  Error variance: \(\hat\sigma^2_{error,v} = \sum e_v^2 / (N - (P+1))\).
        5.  Covariance of linear coefficients: \(\Sigma_{\text{coeffs},v} = \hat\sigma^2_{error,v} (X_{design}(\theta_{exp,v})^T X_{design}(\theta_{exp,v}) + \text{diag}(rep(\lambda_{ridge\_jacobian},P+1)) )^{-1}\).
        6.  SEs for \(\hat\theta_v\) via Delta Method:
            *   The parameters are \(\hat\theta_v = \theta_{exp,v} + \Delta\hat\theta_v\), where \(\Delta\hat\theta_{v,k} = \hat{\text{coeffs}}_{v,k+1} / \hat{\text{coeffs}}_{v,1}\) for \(k=1 \ldots P\).
            *   Since \(\theta_{exp,v}\) is treated as fixed for this SE calculation step, \(SE(\hat\theta_{v,k}) = SE(\Delta\hat\theta_{v,k})\).
            *   The transformation is \(g(\text{coeffs}) = (\text{coeffs}_2/\text{coeffs}_1, \ldots, \text{coeffs}_{P+1}/\text{coeffs}_1)^T\).
            *   Compute the Jacobian of \(g\), \(J_g \in \mathbb{R}^{P \times (P+1)}\), evaluated at \(\hat{\text{coeffs}}_v\).
                For \(k=1 \ldots P\):
                \( \partial (\Delta\hat\theta_{v,k}) / \partial (\hat{\text{coeffs}}_{v,1}) = -\hat{\text{coeffs}}_{v,k+1} / (\hat{\text{coeffs}}_{v,1})^2 \)
                \( \partial (\Delta\hat\theta_{v,k}) / \partial (\hat{\text{coeffs}}_{v,k+1}) = 1 / \hat{\text{coeffs}}_{v,1} \)
                Other partials are 0.
            *   Covariance of \(\Delta\hat\theta_v\): \(\Sigma_{\Delta\theta,v} = J_g \Sigma_{\text{coeffs},v} J_g^T\).
            *   `se_theta_hat_voxel[,v] = sqrt(diag(Sigma_delta_theta_v))`.
    *   Output: Add `se_theta_hat_map = t(se_theta_hat_voxel)` to the list returned by `.parametric_engine`. Enable/disable via `compute_se` argument.

**III. Output Object Enhancements and S3 Methods:**

*   **Ticket SP2-205: Enhance `parametric_hrf_fit` S3 Class**
    *   Add slots for `r_squared`, `residuals`, `parameter_ses`, `convergence_info` to the class definition/constructor.
    *   `estimate_parametric_hrf` populates these new fields from the output of the enhanced `.parametric_engine`.

*   **Ticket SP2-206: Implement `summary.parametric_hrf_fit` Method**
    *   Show distributions (min, median, mean, max) of estimated \(\hat\tau, \hat\sigma, \hat\rho\), \(\hat\beta_0\), \(R^2\).
    *   Report number of global re-centering passes performed, final \(\Delta\theta_0\).

*   **Ticket SP2-207: Implement `coef.parametric_hrf_fit` Method**
    *   Allow `type = "parameters"` (default) to return the V x P matrix of \(\hat\theta_v\).
    *   Allow `type = "amplitude"` to return the V x 1 vector of \(\hat\beta_{0,v}\).
    *   Column names for parameters should come from `hrf_parameter_names()`.

*   **Ticket SP2-208: Implement `fitted.parametric_hrf_fit` and `residuals.parametric_hrf_fit` Methods**
    *   `fitted()`: Reconstructs fitted values \(X_{design}(\hat\theta_v) \hat{\text{coeffs}}_v\) for each voxel. This requires storing or recomputing \(X_{design}(\hat\theta_v)\). *Simplification for Sprint 2: If `residuals_matrix` is stored, `fitted = Y_proj - residuals_matrix`.*
    *   `residuals()`: Returns the stored `residuals_matrix`.

*   **Ticket SP2-209: Implement `plot.parametric_hrf_fit` Method**
    *   For a selected voxel (or mean/median voxel):
        *   Retrieve \(\hat\theta_v\) and \(\hat\beta_{0,v}\).
        *   Generate time points `t_plot` (e.g., `hrf_eval_times`).
        *   Calculate estimated HRF shape: `hrf_shape_est = .lwu_hrf_function(t_plot, params_vector = theta_hat_voxel[,v])`.
        *   Calculate full scaled HRF: `hrf_scaled_est = beta0_voxel[v] * hrf_shape_est`.
        *   Plot `hrf_scaled_est` vs `t_plot`.
        *   Optionally overlay fitted data and actual data for that voxel.

**IV. Unit Tests for Sprint 2:**

*   **Ticket SP2-210: Unit Tests for Global Re-Centering**
    *   Verify that `theta_current_global` converges or stops after `recenter_global_passes`.
    *   Test `theta_seed = "data_median"` behavior.
    *   Check that `R2_voxel` generally improves or stabilizes across re-centering passes.
*   **Ticket SP2-211: Unit Tests for \(R^2\) and Residuals**
    *   Compare calculated \(R^2\) and residuals with known values for simple cases.
*   **Ticket SP2-212: Unit Tests for Standard Error Calculation**
    *   For a simple case with known noise, compare SEs to those from `lm()` on the final linearized model for a single voxel.
    *   Test `compute_se = FALSE` flag.
*   **Ticket SP2-213: Unit Tests for S3 Methods**
    *   Verify `summary`, `coef`, `fitted`, `residuals`, `plot` methods run without error and produce outputs of expected type/dimension.

---

**Deliverables for End of Sprint 2:**
*   `estimate_parametric_hrf` now supports iterative global re-centering and "data_median" seeding.
*   The `.parametric_engine` calculates and returns \(R^2\), residuals, and (optionally) standard errors for \(\hat\theta_v\).
*   The `parametric_hrf_fit` object stores these additional diagnostics.
*   Functional S3 methods: `print`, `summary`, `coef`, `fitted`, `residuals`, `plot`.
*   Unit tests covering all new and enhanced functionality.

This sprint significantly advances the usability and completeness of the parametric estimation, providing key diagnostics and uncertainty measures. The K-means and tiered refinement steps from the original proposal would logically follow in a Sprint 3 if further robustness and adaptation to heterogeneous HRF shapes are desired.