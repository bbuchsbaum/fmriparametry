# Sprint 3: Advanced Robustness Features and Production Readiness

## Sprint Overview

**Duration:** 2 weeks  
**Start Date:** [TBD]  
**End Date:** [TBD]  
**Prerequisites:** Successful completion of Sprint 1 and Sprint 2

### Sprint Goal
Complete the `fmriparametric` package implementation by adding K-means based spatial re-centering, tiered refinement strategies, and comprehensive parallelization. Transform the package into a production-ready tool capable of handling challenging real-world fMRI data with heterogeneous HRF shapes.

### Sprint Objectives
1. Implement K-means based spatial re-centering for heterogeneous data
2. Add tiered refinement queue for challenging voxels
3. Implement local voxel-specific re-centering and Gauss-Newton optimization
4. Complete parallelization across all voxel-wise operations
5. Finalize all S3 methods and visualization tools
6. Comprehensive documentation and performance optimization
7. Production-ready testing and validation

### Building on Previous Sprints
- Leverages iterative global re-centering from Sprint 2
- Extends R², residual, and SE calculations
- Enhances the `parametric_hrf_fit` object with refinement diagnostics
- Builds on existing S3 method infrastructure

### Success Criteria
- [ ] K-means re-centering improves fits for spatially clustered data by >20%
- [ ] Refinement queue correctly identifies and improves difficult voxels
- [ ] Gauss-Newton optimization converges for >90% of hard cases
- [ ] Parallel processing achieves >5x speedup on 8 cores
- [ ] All voxels achieve R² > 0.5 or are flagged as failures
- [ ] Package passes CRAN checks

---

## Sprint 3 Tickets

### Epic 1: K-Means Based Spatial Re-Centering (Est: 2.5 days)

#### Ticket SP3-301: Implement K-Means Re-Centering Engine
**Priority:** Critical  
**Estimate:** 16 hours  
**Dependencies:** SP2-201 (global re-centering)

**Description:**  
Implement K-means clustering on parameter estimates to identify spatial clusters with distinct HRF characteristics, then re-center within each cluster.

**Implementation Steps:**

1. **Modify `.parametric_engine` to support K-means re-centering:**
   ```r
   .kmeans_recentering <- function(
     theta_hat_voxel,
     R2_voxel,
     Y_proj,
     S_target_proj,
     recenter_kmeans_passes = 2,
     kmeans_k = 5,
     r2_threshold_kmeans = 0.7,
     ...
   )
   ```

2. **K-means clustering workflow:**
   - Select good voxels based on R² threshold
   - Transform parameters for clustering (e.g., log(sigma))
   - Perform K-means clustering with multiple starts
   - Use cluster centers as expansion points
   - Re-fit voxels within each cluster

3. **Cluster processing:**
   - Apply single Taylor pass per cluster
   - Update only if R² improves
   - Track cluster assignments for diagnostics

**Acceptance Criteria:**
- [ ] K-means identifies meaningful parameter clusters
- [ ] Each cluster's centroid used as expansion point
- [ ] Updates only improve fit quality (R²)
- [ ] Handles edge cases (few good voxels, small clusters)
- [ ] Performance scales linearly with voxel count
- [ ] Reproducible results with seed setting

#### Ticket SP3-302: Add K-Means Options to Main Function
**Priority:** High  
**Estimate:** 4 hours  
**Dependencies:** SP3-301

**Description:**  
Integrate K-means re-centering into the main `estimate_parametric_hrf` function.

**Implementation:**
```r
estimate_parametric_hrf <- function(
  ...,
  recenter_kmeans_passes = 2,
  kmeans_k = 5,
  r2_threshold_kmeans = 0.7,
  ...
)
```

**Acceptance Criteria:**
- [ ] New parameters properly documented
- [ ] Integration with existing workflow
- [ ] Default values provide good performance
- [ ] User can disable K-means with passes = 0

### Epic 2: Tiered Refinement Queue System (Est: 3 days)

#### Ticket SP3-303: Implement Refinement Queue Classification
**Priority:** Critical  
**Estimate:** 6 hours  
**Dependencies:** SP3-301

**Description:**  
Classify voxels into refinement tiers based on fit quality and uncertainty metrics.

**Implementation:**
```r
.classify_refinement_queue <- function(
  R2_voxel,
  se_theta_hat_voxel = NULL,
  refinement_opts = list(
    apply_refinement = TRUE,
    r2_threshold_hard = 0.3,
    r2_threshold_moderate = 0.7,
    se_threshold_hard = 0.5,
    se_threshold_moderate = 0.3
  )
)
```

**Queue Definitions:**
- **"easy"**: Good fit (R² > 0.7, low SE) - no refinement needed
- **"moderate_local_recenter"**: Moderate fit - benefits from local re-centering
- **"hard_GN"**: Poor fit - requires Gauss-Newton optimization

**Acceptance Criteria:**
- [ ] Correctly classifies voxels based on R² and SE thresholds
- [ ] Handles missing SE information gracefully
- [ ] Queue assignments are mutually exclusive
- [ ] Summary statistics reported
- [ ] Customizable thresholds via options

#### Ticket SP3-304: Implement Local Re-centering for Moderate Voxels
**Priority:** High  
**Estimate:** 8 hours  
**Dependencies:** SP3-303

**Description:**  
For moderately difficult voxels, perform local re-centering using voxel-specific expansion points.

**Algorithm:**
1. Identify moderate queue voxels
2. Use each voxel's current estimate as local expansion point
3. Perform single Taylor pass for each voxel
4. Update only if R² improves
5. Track refinement success

**Acceptance Criteria:**
- [ ] Uses voxel's current estimate as expansion point
- [ ] Processes only moderate queue voxels
- [ ] Updates only when fit improves
- [ ] Updates queue labels after refinement
- [ ] Reports improvement statistics

#### Ticket SP3-305: Implement Gauss-Newton Refinement for Hard Voxels
**Priority:** High  
**Estimate:** 12 hours  
**Dependencies:** SP3-303

**Description:**  
Implement full nonlinear Gauss-Newton optimization for the most challenging voxels.

**Implementation:**
```r
.gauss_newton_refinement <- function(
  Y_proj,
  S_target_proj,
  scan_times,
  hrf_eval_times,
  theta_hat_voxel,
  queue_labels,
  hrf_interface,
  max_iter_gn = 5,
  tol_gn = 1e-4,
  lambda_ridge = 0.01,
  ...
)
```

**Algorithm Steps:**
1. Identify hard queue voxels
2. For each hard voxel:
   - Initialize with current parameter estimates
   - Iteratively update using Gauss-Newton method
   - Apply parameter bounds at each iteration
   - Re-estimate amplitude with new parameters
   - Check convergence criteria
3. Update results only if converged or improved

**Acceptance Criteria:**
- [ ] Implements proper Gauss-Newton algorithm
- [ ] Handles convergence criteria correctly
- [ ] Re-estimates amplitude at each iteration
- [ ] Applies parameter bounds
- [ ] Tracks convergence status
- [ ] Improves R² for majority of hard cases

#### Ticket SP3-306: Integrate Tiered Refinement into Main Function
**Priority:** High  
**Estimate:** 6 hours  
**Dependencies:** SP3-304, SP3-305

**Description:**  
Integrate the complete refinement pipeline into the main estimation workflow.

**Implementation:**
```r
estimate_parametric_hrf <- function(
  ...,
  refinement_opts = list(
    apply_refinement = TRUE,
    r2_threshold_hard = 0.3,
    r2_threshold_moderate = 0.7,
    max_iter_gn = 5
  ),
  ...
)
```

**Acceptance Criteria:**
- [ ] Complete refinement pipeline integrated
- [ ] Proper error handling for edge cases
- [ ] User can disable refinement
- [ ] Refinement diagnostics stored in output

### Epic 3: Parallelization and Performance (Est: 2 days)

#### Ticket SP3-307: Implement Comprehensive Parallel Processing
**Priority:** High  
**Estimate:** 12 hours  
**Dependencies:** SP3-304, SP3-305

**Description:**  
Implement parallel processing across all computationally intensive voxel-wise operations.

**Implementation Strategy:**
```r
.setup_parallel_backend <- function(n_cores = NULL, verbose = TRUE)

.parallel_voxel_processing <- function(
  voxel_indices,
  process_function,
  ...,
  chunk_size = "auto",
  progress = TRUE
)
```

**Parallelization Points:**
1. K-means cluster processing
2. Moderate voxel local re-centering  
3. Gauss-Newton optimization for hard voxels
4. Standard error calculation

**Acceptance Criteria:**
- [ ] Parallel backend setup/teardown works correctly
- [ ] All voxel-wise operations parallelized
- [ ] Results identical to sequential processing
- [ ] Linear speedup to at least 4 cores
- [ ] Memory usage stays within reasonable bounds
- [ ] Progress reporting works in parallel mode

#### Ticket SP3-308: Add Parallel Options to Main Function
**Priority:** Medium  
**Estimate:** 4 hours  
**Dependencies:** SP3-307

**Description:**  
Add parallel processing options to the main function with sensible defaults.

**Implementation:**
```r
estimate_parametric_hrf <- function(
  ...,
  n_cores = NULL,  # NULL = auto-detect
  progress = TRUE,
  ...
)
```

**Acceptance Criteria:**
- [ ] Auto-detection of available cores
- [ ] Graceful fallback if parallel packages unavailable
- [ ] User can control core count
- [ ] Progress reporting optional

### Epic 4: Enhanced S3 Methods and Visualization (Est: 1.5 days)

#### Ticket SP3-309: Enhance S3 Methods for Refinement Information
**Priority:** Medium  
**Estimate:** 8 hours  
**Dependencies:** SP3-306

**Description:**  
Enhance S3 methods to display refinement information and provide production-ready visualizations.

**Enhanced Methods:**
1. **Summary method** - add refinement statistics
2. **Plot method** - add refinement diagnostics
3. **Print method** - show refinement summary

**Implementation:**
```r
summary.parametric_hrf_fit <- function(object, ...) {
  # Base summary + refinement information
}

plot.parametric_hrf_fit <- function(
  x,
  type = c("hrf", "parameters", "diagnostic", "refinement"),
  ...
)
```

**Acceptance Criteria:**
- [ ] Summary includes refinement statistics
- [ ] Plot method shows refinement diagnostics
- [ ] All methods handle missing refinement info
- [ ] Visualizations are publication-ready
- [ ] Methods maintain backward compatibility

#### Ticket SP3-310: Add Standard Error Recalculation
**Priority:** Medium  
**Estimate:** 4 hours  
**Dependencies:** SP3-305, SP2-204

**Description:**  
Recalculate standard errors for refined voxels using final parameter estimates.

**Implementation:**
```r
.recalculate_se_after_refinement <- function(
  theta_hat_voxel,
  Y_proj,
  queue_labels,
  hrf_interface,
  ...
)
```

**Acceptance Criteria:**
- [ ] Recalculates SEs only for refined voxels
- [ ] Uses final parameter estimates as expansion points
- [ ] Preserves existing SEs for non-refined voxels
- [ ] Computation time proportional to refined voxel count

### Epic 5: Documentation and Testing (Est: 2 days)

#### Ticket SP3-311: Comprehensive Unit Tests
**Priority:** High  
**Estimate:** 12 hours  
**Dependencies:** All Sprint 3 features

**Description:**  
Create comprehensive unit tests for all new functionality.

**Test Categories:**
1. **K-means re-centering tests** - clustered synthetic data
2. **Tiered refinement tests** - challenging voxel scenarios
3. **Parallel processing tests** - result consistency and speedup
4. **Integration tests** - end-to-end workflows

**Acceptance Criteria:**
- [ ] K-means improves accuracy for clustered data
- [ ] Queue classification matches expected difficulty
- [ ] Parallel results identical to sequential
- [ ] Full pipeline completes without errors
- [ ] Test coverage >85%

#### Ticket SP3-312: Complete Documentation
**Priority:** High  
**Estimate:** 4 hours  
**Dependencies:** All Sprint 3 features

**Description:**  
Complete documentation for all advanced features and prepare for CRAN submission.

**Documentation Tasks:**
1. Update main function documentation
2. Create advanced features vignette
3. Add performance guide
4. Update package documentation

**Acceptance Criteria:**
- [ ] All functions fully documented with examples
- [ ] Advanced refinement vignette complete
- [ ] Performance guide with benchmarks
- [ ] Package passes R CMD check with no warnings
- [ ] Examples run in reasonable time

---

## Sprint 3 Deliverables

### Code Deliverables
1. **Advanced Recentering**
   - [x] K-means spatial clustering and recentering
   - [x] Cluster-specific expansion points
   - [x] Multiple recentering passes

2. **Tiered Refinement System**
   - [x] Automatic voxel classification by difficulty
   - [x] Local voxel-specific recentering
   - [x] Full Gauss-Newton optimization
   - [x] Refinement tracking and diagnostics

3. **Performance Optimization**
   - [x] Comprehensive parallelization
   - [x] Memory-efficient processing
   - [x] Progress reporting

4. **Enhanced S3 Methods**
   - [x] Refinement-aware summary and plot methods
   - [x] Advanced diagnostic visualizations
   - [x] Queue analysis tools

5. **Production Readiness**
   - [x] Comprehensive test suite (>85% coverage)
   - [x] Complete documentation
   - [x] Performance benchmarks
   - [x] CRAN-ready package

### Technical Achievements
- K-means improves fits for spatially heterogeneous data by 20-40%
- Tiered refinement recovers >90% of difficult voxels
- Gauss-Newton converges for 85-95% of hard cases
- 5-8x speedup with parallel processing on 8 cores
- Memory usage < 2x data size even with all features

### Performance Targets
- 1000 voxels full pipeline: <1 minute (8 cores)
- 10k voxels full pipeline: <5 minutes (8 cores)
- 100k voxels full pipeline: <30 minutes (8 cores)
- Memory overhead: <50% of data size

### Quality Metrics
- Test coverage: >85%
- R CMD check: 0 errors, 0 warnings, 0 notes
- Documentation: 100% complete
- Examples: All run in <5 seconds

---

## Post-Sprint 3 Summary

### Package Status
The `fmriparametric` package is now a complete, production-ready tool for parametric HRF estimation featuring:

1. **Robust Estimation**
   - Multiple initialization strategies
   - Adaptive refinement for difficult cases
   - Comprehensive convergence handling

2. **Scientific Features**
   - Full uncertainty quantification
   - Spatial heterogeneity handling
   - Diagnostic tools for quality assessment

3. **Performance**
   - Scales to whole-brain analyses
   - Efficient parallel processing
   - Memory-conscious implementation

4. **Usability**
   - Clear documentation and vignettes
   - Intuitive S3 methods
   - Informative visualizations

### Future Enhancements
1. **Additional HRF Models**
   - Double-gamma implementation
   - Flexible basis functions
   - Custom user-defined models

2. **Advanced Features**
   - Bayesian estimation options
   - Multi-session support
   - Real-time processing mode

3. **Integration**
   - BIDS app wrapper
   - fMRIPrep compatibility
   - Report generation