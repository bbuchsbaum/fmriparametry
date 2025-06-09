# `fmriparametric` – Extensible Parametric HRF Estimation

## Executive Summary

The `fmriparametric` R package provides robust and efficient tools for estimating parameters of parametric Hemodynamic Response Function (HRF) models from fMRI data. Using an iterative linear Taylor approximation method, it enables voxel-wise estimation of interpretable HRF parameters (e.g., lag, width, undershoot amplitude) with uncertainty quantification. The package is designed for extensibility, initially focusing on the Lag-Width-Undershoot (LWU) model while providing a framework for incorporating additional parametric models.

**Key Features:**
- Fast voxel-wise parametric HRF estimation using Taylor approximation
- Robust initialization with adaptive seeding strategies
- Standard error estimation for all parameters
- Extensible architecture for multiple HRF models
- Integration with the `fmrireg` ecosystem

---

## 1. Core Philosophy & Goals
`fmriparametric` aims to provide robust and efficient tools for estimating parameters of *pre-defined parametric HRF models* (like the Lag-Width-Undershoot model) on a voxel-wise basis from fMRI data. It will leverage the `fmrireg` package for core fMRI data structures, event modeling, and HRF function definitions.

### Core Goals:
- **Parametric HRF Estimation:** Enable direct estimation of interpretable HRF parameters (e.g., lag, width, undershoot amplitude).
- **Specific Model Focus (Initial):** Initially, focus on the 3-parameter LWU model (`fmrireg::hrf_lwu`).
- **Extensibility:** The core fitting engine should be designed to accommodate other parametric HRF models from `fmrireg` (or user-defined ones that conform to a specific interface) in the future.
- **Efficiency:** Employ an iterative linear Taylor approximation method for computational speed across many voxels.
- **Robustness:** Incorporate adaptive seeding for the Taylor expansion point \(\theta_0\) and optional tiered refinement strategies for challenging voxels.
- **Uncertainty Quantification:** Provide standard errors for the estimated HRF parameters.
- **Leverage `fmrireg`:** Maximize use of `fmrireg` for:
    - HRF function definitions (e.g., `fmrireg::hrf_lwu`).
    - Taylor basis construction (e.g., `fmrireg::hrf_basis_lwu`).
    - Event modeling (`fmrireg::event_model`).
    - Data representation (`fmrireg::fmri_dataset`, `fmrireg::matrix_dataset`).
    - Confound regression utilities if applicable.
- **User-Friendliness:** Offer a clear interface for specifying models and interpreting results.

## 2. The Core Method: Iterative Linear Taylor Approximation

The parameters \(\theta = (\theta_1, \ldots, \theta_P)\) of a given \(P\)-parameter HRF model \(h(t; \theta)\) are estimated voxel-wise by linearizing \(h(t; \theta)\) around an expansion point \(\theta_0\).

### 2.1 Example Parametric HRF Model: Lag-Width-Undershoot (LWU)
The primary initial target for this package is the LWU model, which is assumed to be available from the `fmrireg` package (as per `repomix-output.txt` under `R/hrf-functions.R`).
#### Functional Form:
 \(h(t;\,\tau,\sigma,\rho)= e^{-\frac{(t-\tau)^2}{2\sigma^{2}}} - \rho\,e^{-\frac{\bigl(t-\tau-2\sigma\bigr)^2}{2(1.6\sigma)^{2}}}\)

#### Parameters: \(\theta_{LWU} = (\tau, \sigma, \rho)\). So, for LWU, \(P=3\).

#### Normalization: The `fmrireg::hrf_lwu` function includes a `normalize` argument (e.g., "none", "height", "area"). For the internal steps of the Taylor expansion fitting process, `normalize = "none"` will be used for the HRF components and their derivatives. The overall amplitude of the fitted HRF will be captured by a separate \(\beta_0\) parameter.

#### Parameter Bounds: `fmrireg::hrf_lwu` itself has inherent or recommended bounds (e.g., \(\sigma > 0.05\), \(0 \le \rho \le 1.5\)). These, along with user-specified `theta_bounds`, will be respected during the estimation process.

### 2.2 Taylor Basis Construction for LWU Model
The Taylor expansion requires the HRF value and its partial derivatives with respect to its parameters, all evaluated at the current expansion point \(\theta_0\). The `fmrireg::hrf_basis_lwu` function (assumed available from `fmrireg`, as per `repomix-output.txt` under `R/hrf-functions.R`) provides this for the LWU model.
For a given expansion point \(\theta_0 = (\tau_0, \sigma_0, \rho_0)\) and a set of time points \(t_{hrf}\) for HRF evaluation, it constructs:
    \(X_{taylor\_basis}(\theta_0; t_{hrf}) = [h(t_{hrf};\theta_0), \frac{\partial h}{\partial \tau}|_{\theta_0}, \frac{\partial h}{\partial \sigma}|_{\theta_0}, \frac{\partial h}{\partial \rho}|_{\theta_0}] \in \mathbb{R}^{T_{hrf} \times (3+1)}\).

The `normalize_primary` argument of `fmrireg::hrf_basis_lwu` will be set to `"none"` to ensure that the derivatives are scaled appropriately relative to the primary HRF component \(h(t_{hrf};\theta_0)\) for the linear approximation.

### 2.3 Generic Parametric HRF Interface

**Design Note:** This interface enables extensibility to support multiple parametric HRF models beyond LWU.
To support extensibility beyond the LWU model, any parametric HRF intended for use with this engine must adhere to a defined interface, or be wrapped by functions that provide this interface. This interface consists of:

#### a. `hrf_function(t, params_vector, ...)`
- **Purpose:** Evaluates the parametric HRF.
- **Inputs:**
    - `t`: Numeric vector of time points.
    - `params_vector`: Numeric vector containing the current values of the HRF's \(P\) parameters, in a consistent order.
    - `...`: Any other fixed arguments required by the specific HRF (e.g., model constants not being estimated).
- **Output:** A numeric vector (if \(P_{basis}=1\)) or matrix (if \(P_{basis}>1\), e.g., for multi-component parametric HRFs) of HRF values at time points `t`. For the Taylor expansion, we typically assume \(P_{basis}=1\) for the primary \(h(t;\theta)\) component.
- **Example:** `hrf_lwu(t, tau=params_vector[1], sigma=params_vector[2], rho=params_vector[3], normalize="none")`.

#### b. `hrf_taylor_basis_function(params_vector0, t_hrf_eval, ...)`
- **Purpose:** Constructs the Taylor expansion basis (HRF value and partial derivatives).
- **Inputs:**
    - `params_vector0`: Numeric vector of the \(P\) parameter values at the expansion point \(\theta_0\).
    - `t_hrf_eval`: Numeric vector of time points for evaluating the basis.
    - `...`: Other fixed arguments for the HRF.
- **Output:** A numeric matrix of dimension `length(t_hrf_eval) x (P+1)`.
    - Column 1: \(h(t_{hrf\_eval}; \theta_0)\) (the HRF evaluated at \(\theta_0\)).
    - Column 2 to P+1: Partial derivatives \(\partial h / \partial \theta_j\) evaluated at \(\theta_0\), for each of the \(P\) parameters.
- **Example:** `fmrireg::hrf_basis_lwu(theta0=params_vector0, t=t_hrf_eval, normalize_primary="none")`.

#### c. `hrf_parameter_names(hrf_object_or_name)`
- **Purpose:** Returns the ordered names of the \(P\) parameters.
- **Input:** The HRF object itself or its registered name.
- **Output:** A character vector of length \(P\) (e.g., for LWU: `c("tau", "sigma", "rho")`).

#### d. `hrf_default_seed(hrf_object_or_name)`
- **Purpose:** Returns a default starting seed \(\theta_0\) for the parameters.
- **Input:** The HRF object or its registered name.
- **Output:** A numeric vector of length \(P\).

#### e. `hrf_default_bounds(hrf_object_or_name)`
- **Purpose:** Returns default lower and upper bounds for the parameters.
- **Input:** The HRF object or its registered name.
- **Output:** A list `list(lower = numeric_vector_P, upper = numeric_vector_P)`.

**Implementation Note:** This interface could be implemented via S3 methods dispatched on the `fmrireg::HRF` object if its class structure is extended to signify "parametric" nature and store necessary metadata. Alternatively, `fmriparametric` could maintain a registry mapping HRF names/objects to these interface functions.

## 3. Key Functions in `fmriparametric`

### 3.1 Main User-Facing Function: `estimate_parametric_hrf`

This is the primary function users will interact with to perform parametric HRF estimation.

#### Function Signature:
```r
estimate_parametric_hrf(
  fmri_data,
  event_model,
  parametric_hrf = "lwu",
  theta_seed = NULL,
  theta_bounds = NULL,
  max_iter = 3,
  relative_tolerance = 0.01,
  absolute_tolerance = 0.001,
  adaptive_seeding = TRUE,
  kmeans_recentering = TRUE,
  tiered_refinement = FALSE,
  mask = NULL,
  verbose = FALSE
)
```

#### Workflow:

1. **Resolve Parametric HRF & Parameters:**
   - Determine the specific parametric HRF model (e.g., LWU) and its properties (`P`, parameter names, default seed/bounds) using the interface functions.

2. **Prepare Inputs** (via `.prepare_parametric_inputs` internal helper):
   - Extract and project BOLD data
   - Construct event design matrices
   - Handle confound regression if needed
   - Generate timing vectors for HRF evaluation

3. **Core Parametric HRF Estimation** (via `.parametric_engine` internal helper):
   - Initialize with adaptive seeding strategies
   - Perform iterative Taylor approximation updates
   - Apply optional K-means recentering for robust initialization
   - Implement tiered refinement for challenging voxels
   - Calculate standard errors via information matrix

4. **Construct Output Object** (`parametric_hrf_fit` S3 class):
   - Store estimated parameters and standard errors
   - Include convergence diagnostics
   - Preserve model metadata and fitting options
   - Provide methods for visualization and summary

### 3.2 Internal Helper Functions

#### `.prepare_parametric_inputs`
Prepares data structures for the parametric fitting engine:
- Projects BOLD data and handles masking
- Constructs stimulus design matrices
- Manages confound regression
- Returns standardized inputs for the engine

#### `.parametric_engine`
The core estimation workhorse that:
- Implements the iterative Taylor approximation algorithm
- Handles multiple initialization strategies
- Performs voxel-wise optimization with convergence checks
- Calculates parameter uncertainties

### 3.3 S3 Methods for `parametric_hrf_fit`

```r
# Summary method
summary.parametric_hrf_fit <- function(object, ...) 

# Plotting method  
plot.parametric_hrf_fit <- function(x, voxel_indices = NULL, ...)

# Extraction methods
coef.parametric_hrf_fit <- function(object, ...)
vcov.parametric_hrf_fit <- function(object, voxel_index, ...)

# Prediction method
predict.parametric_hrf_fit <- function(object, newdata = NULL, ...)
```

## 4. Technical Implementation Details

### 4.1 Iterative Taylor Approximation Algorithm

For each voxel, the algorithm:

1. **Initialize** θ₀ using one of several strategies:
   - User-provided seed
   - Default model values
   - Data-driven initialization (peak detection)
   - K-means clustering of initial estimates

2. **Iterate** until convergence:
   ```
   a. Construct Taylor basis X_taylor at current θ₀
   b. Solve linear system: β = (X'X)⁻¹X'y
   c. Update parameters: θ_new = θ₀ + [β₂, β₃, ..., β_P+1]'
   d. Apply bounds: θ_new = clip(θ_new, bounds)
   e. Check convergence: ||θ_new - θ₀|| < tol
   f. Update θ₀ = θ_new
   ```

3. **Compute uncertainties** via the information matrix

### 4.2 Robustness Features

#### Adaptive Seeding
- Peak-based initialization for lag parameter
- Width estimation from FWHM
- Undershoot detection from negative lobes

#### K-means Recentering
- Cluster voxels by initial parameter estimates
- Use cluster centers as refined seeds
- Reduces sensitivity to outliers

#### Tiered Refinement
- Identify poorly converged voxels
- Re-fit with alternative strategies
- Multiple initialization attempts

## 5. Dependencies

### Required Packages:
- **`fmrireg`** (>= 0.2.0): Core fMRI functionality
- **`Matrix`**: Sparse matrix operations
- **`RcppEigen`**: Fast linear algebra

### Suggested Packages:
- **`future`**: Parallel processing
- **`progressr`**: Progress reporting
- **`ggplot2`**: Visualization

## 6. Development Roadmap

### Phase 1: Core Infrastructure (Weeks 1-2)
- Set up package structure and testing framework
- Implement generic HRF interface
- Create data preparation utilities

### Phase 2: LWU Implementation (Weeks 3-4)
- Implement Taylor approximation for LWU model
- Add adaptive seeding strategies
- Create basic S3 methods

### Phase 3: Robustness & Performance (Weeks 5-6)
- Add K-means recentering
- Implement tiered refinement
- Optimize computational bottlenecks
- Add parallel processing support

### Phase 4: Polish & Documentation (Week 7)
- Comprehensive testing suite
- User documentation and vignettes
- Performance benchmarking
- CRAN submission preparation

## 7. Example Usage

```r
library(fmriparametric)
library(fmrireg)

# Load fMRI data
fmri_data <- fmrireg::read_fmri("subject01_task.nii.gz")

# Define event model
events <- data.frame(
  onset = c(10, 40, 70, 100),
  duration = c(2, 2, 2, 2),
  trial_type = "stimulus"
)
event_model <- fmrireg::event_model(
  onset ~ trial_type,
  data = events,
  sampling_rate = 1/2  # TR = 2s
)

# Estimate parametric HRF
fit <- estimate_parametric_hrf(
  fmri_data = fmri_data,
  event_model = event_model,
  parametric_hrf = "lwu",
  adaptive_seeding = TRUE,
  kmeans_recentering = TRUE
)

# Examine results
summary(fit)
plot(fit, voxel_indices = c(1000, 2000, 3000))

# Extract parameters for specific voxel
voxel_params <- coef(fit)[1000, ]
voxel_se <- sqrt(diag(vcov(fit, voxel_index = 1000)))
```

## 8. Potential Challenges & Mitigation

### Technical Challenges:
1. **Convergence issues**: Mitigated by robust initialization and bounds
2. **Computational speed**: Addressed via vectorization and parallelization
3. **Memory usage**: Managed through chunked processing for large datasets

### User Experience:
1. **Parameter interpretation**: Clear documentation with physiological context
2. **Model selection**: Guidance on choosing appropriate HRF models
3. **Quality control**: Built-in diagnostics and visualization tools

## 9. Future Extensions

1. **Additional HRF Models**:
   - Double-gamma model
   - FIR-constrained models
   - Physiologically-informed models

2. **Advanced Features**:
   - Multi-run/session estimation
   - Group-level analysis tools
   - Bayesian estimation options
   - Integration with preprocessing pipelines

3. **Performance Enhancements**:
   - GPU acceleration for large datasets
   - Approximate methods for real-time analysis
   - Sparse estimation for efficiency

## 10. Conclusion

The `fmriparametric` package fills a critical gap in the fMRI analysis ecosystem by providing efficient, robust tools for parametric HRF estimation. By building on the `fmrireg` foundation and implementing state-of-the-art optimization techniques, it enables researchers to extract physiologically meaningful parameters from their fMRI data while maintaining computational efficiency for whole-brain analyses.

The extensible architecture ensures that the package can grow to accommodate new HRF models and analysis techniques as the field evolves, making it a valuable long-term investment for the neuroimaging community.

## find information about fmrireg and neuroim2 in cheatsheets

data-raw/fmrireg_cheatsheet.md
data-raw/neuroim2_cheatsheet.md