# fmriparametric 0.2.0

## Major improvements

* Complete package stabilization and deployment readiness
* Fixed critical stack overflow issues in parametric engine
* Implemented robust, dependency-free simple engine for core functionality
* Significantly improved performance (50x5 voxels in 0.046s, 200x50 voxels in 0.156s)

## Bug fixes

* Fixed dimensional errors in global refinement code
* Resolved function signature mismatches across deprecated versions
* Fixed parameter bounds validation for LWU model
* Corrected numerical robustness issues in engineering standards

## Documentation

* Added comprehensive README with examples and performance benchmarks
* Fixed all roxygen2 documentation warnings
* Updated package documentation to modern standards
* Added proper NAMESPACE management via roxygen2

## Testing

* Added end-to-end workflow validation tests
* Created comprehensive integration test suite
* Fixed critical test failures in engineering standards
* Improved test stability across platforms

## Internal changes

* Consolidated multiple function versions into single implementations
* Added global variable declarations to suppress R CMD check notes
* Updated .Rbuildignore to exclude development files
* Improved error handling and input validation

# fmriparametric 0.1.0

* Initial release
* Basic implementation of parametric HRF estimation using LWU model
* Core Taylor approximation algorithm
* S3 methods for parametric_hrf_fit class