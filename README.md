# fmriparametric

<!-- badges: start -->
[![R-CMD-check](https://github.com/bbuchsbaum/fmriparametry/workflows/R-CMD-check/badge.svg)](https://github.com/bbuchsbaum/fmriparametry/actions)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

`fmriparametric` provides robust and efficient tools for estimating parameters of parametric Hemodynamic Response Function (HRF) models from fMRI data. Using an iterative linear Taylor approximation method, it enables voxel-wise estimation of interpretable HRF parameters with uncertainty quantification.

## Features

- **Fast parameter estimation**: Efficient implementation using sparse matrices and vectorized operations
- **Flexible HRF models**: Currently supports the Lag-Width-Undershoot (LWU) model with extensible architecture
- **Robust estimation**: Ridge regularization and parameter bounds enforcement
- **Uncertainty quantification**: Standard error estimation via delta method
- **Integration**: Works seamlessly with the `fmrireg` ecosystem

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("bbuchsbaum/fmriparametry")
```

## Quick Start

```r
library(fmriparametric)
library(fmrireg)

# Load your fMRI data (timepoints x voxels)
fmri_data <- matrix(rnorm(200 * 100), nrow = 200, ncol = 100)

# Create event model (e.g., stimulus onsets)
event_model <- matrix(0, nrow = 200, ncol = 1)
event_model[c(20, 60, 100, 140, 180), 1] <- 1

# Estimate HRF parameters
fit <- estimate_parametric_hrf(
  fmri_data = fmri_data,
  event_model = event_model,
  parametric_hrf = "lwu"
)

# View results
print(fit)
summary(fit)

# Extract parameters
params <- coef(fit)  # Returns tau, sigma, rho for each voxel
```

### Single-Voxel Sanity Check

For quick verification that the estimation machinery is functioning, the
package exposes `single_voxel_sanity_check()`. It fits a single time series
against a one-condition design.

```r
y <- rnorm(40)
onsets <- rep(0, 40); onsets[c(10, 30)] <- 1
res <- single_voxel_sanity_check(y, onsets)
```

## Documentation

Full function reference and articles are available at
<https://bbuchsbaum.github.io/fmriparametry>.

## HRF Models

### Lag-Width-Undershoot (LWU) Model

The LWU model parameterizes the HRF with three interpretable parameters:

- **τ (tau)**: Lag of the peak response (seconds)
- **σ (sigma)**: Width/duration of the response (seconds)  
- **ρ (rho)**: Undershoot amplitude ratio

Default bounds:
- τ: [0, 20] seconds
- σ: [0.05, 10] seconds
- ρ: [0, 1.5]

## Advanced Usage

### Custom Seeds and Bounds

```r
# Custom parameter seeds
custom_seed <- c(tau = 5, sigma = 2, rho = 0.3)

# Custom bounds
custom_bounds <- list(
  lower = c(2, 0.5, 0),
  upper = c(15, 8, 1.2)
)

fit <- estimate_parametric_hrf(
  fmri_data = fmri_data,
  event_model = event_model,
  parametric_hrf = "lwu",
  theta_seed = custom_seed,
  theta_bounds = custom_bounds
)
```

### Controlling Global Refinement

By default, Stage 3 global iterative refinement runs automatically. To disable
it across all calls, set the package option `fmriparametric.refine_global` to
`FALSE`:

```r
options(fmriparametric.refine_global = FALSE)
```

Re-enable by setting the option back to `TRUE`.

### Working with fmrireg Objects

```r
# Using fmrireg data structures
fmri_dataset <- fmrireg::fmri_dataset(
  scans = fmri_data,
  mask = rep(TRUE, 100),
  TR = 2.0
)

event_model <- fmrireg::event_model(
  onset = c(20, 60, 100, 140, 180),
  blockids = rep(1, 5),
  durations = rep(0, 5)
)

fit <- estimate_parametric_hrf(
  fmri_data = fmri_dataset,
  event_model = event_model,
  parametric_hrf = "lwu"
)
```

### Using the Ultimate Wrapper

For convenience, the package also provides
`estimate_parametric_hrf_ultimate()`. This helper mirrors the main
function but sets sensible defaults for iterative refinement and other
advanced options.

```r
fit <- estimate_parametric_hrf_ultimate(
  fmri_data = fmri_dataset,
  event_model = event_model,
  iterative_recentering = FALSE,
  verbose = FALSE
)
```

## Performance

The package is optimized for whole-brain analysis with:

- Vectorized operations across voxels
- Sparse matrix computations
- Efficient convolution using FFT
- Optional parallel processing (future releases)

Typical performance:
- 200 timepoints × 50 voxels: ~0.15 seconds
- 500 timepoints × 1000 voxels: ~3-5 seconds

## Methods

The package implements parametric HRF estimation using Taylor approximation:

```
h(t; θ) ≈ h(t; θ₀) + Σᵢ ∂h/∂θᵢ|θ₀ (θᵢ - θ₀ᵢ)
```

This linearization allows efficient parameter estimation through standard linear regression with appropriate regularization.

## Vignettes

Detailed walk-throughs of package features are provided in the package vignettes:

- `vignettes/introduction.Rmd` – getting started and basic concepts
- `vignettes/advanced-features.Rmd` – advanced options and large-scale analyses
- `vignettes/performance-guide.Rmd` – benchmarking and tuning tips

Use `browseVignettes("fmriparametric")` within R to view them after installation.

## Contributing

Contributions are welcome! Please file bug reports and feature requests as GitHub issues.

## License

MIT © Bradley Buchsbaum

## References

- Glover, G. H. (1999). Deconvolution of impulse response in event-related BOLD fMRI. NeuroImage, 9(4), 416-429.
- Lindquist, M. A., & Wager, T. D. (2007). Validity and power in hemodynamic response modeling: A comparison study and a new approach. Human Brain Mapping, 28(8), 764-784.