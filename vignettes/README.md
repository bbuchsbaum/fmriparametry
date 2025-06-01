# fmriparametric Vignettes

The `fmriparametric` package includes three vignettes to help you understand and use the package:

## 1. Quick Start Guide
**File:** `quick-start.Rmd`

A brief introduction to get you up and running quickly. Covers:
- Installation
- Basic usage with minimal code
- Key options and parameters
- Tips for success

Perfect for users who want to start analyzing data immediately.

## 2. Understanding Parametric HRF Estimation
**File:** `parametric-hrf-estimation.Rmd`

A comprehensive, pedagogical guide that explains:
- Why parametric HRF models are useful
- How the Taylor approximation method works
- Visual explanations and intuitions
- Complete worked examples
- Interpretation of results
- Common issues and solutions

Ideal for users who want to understand the method deeply.

## 3. Technical Details
**File:** `technical-details.Rmd`

Mathematical and implementation details for advanced users:
- Complete mathematical derivations
- Derivatives for the LWU model
- Numerical considerations
- Performance optimizations
- Convergence properties

For developers and researchers who need the mathematical foundations.

## Building the Vignettes

To build all vignettes:

```r
devtools::build_vignettes()
```

To view a specific vignette after building:

```r
vignette("quick-start", package = "fmriparametric")
vignette("parametric-hrf-estimation", package = "fmriparametric")
vignette("technical-details", package = "fmriparametric")
```
