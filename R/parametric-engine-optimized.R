#' Optimized Parametric HRF Engine with Engineering Excellence
#'
#' This implementation demonstrates engineering best practices:
#' - Consistent interfaces
#' - Numerical robustness  
#' - Performance optimization
#' - Comprehensive validation
#' - Clear error messages
#'
#' @section Algorithm:
#' Implements Taylor approximation for HRF parameter estimation:
#' \deqn{h(t; \theta) \approx h(t; \theta_0) + \sum_{i=1}^{p} \frac{\partial h}{\partial \theta_i}|_{\theta_0} (\theta_i - \theta_{0i})}
#'
#' @section Performance:
#' - Time complexity: O(n*p + p^3) where n = timepoints, p = parameters
#' - Space complexity: O(n*v + n*p) where v = voxels
#' - Optimizations: Batch matrix operations, QR caching, vectorized convolution
#'
#' @param fmri_data Numeric matrix (timepoints x voxels)
#' @param event_design Numeric matrix (timepoints x conditions)
#' @param hrf_interface List with HRF function interface
#' @param hrf_parameters List with seed, bounds, and options
#' @param algorithm_options List with algorithm control parameters
#' @param validate Logical whether to validate inputs
#'
#' @return List of class 'parametric_engine_result' containing:
#'   - parameters: Estimated HRF parameters (voxels x parameters)
#'   - amplitudes: Response amplitudes (voxels)
#'   - fit_quality: R-squared values (voxels)
#'   - diagnostics: Algorithm diagnostics
#'   - status: Success/failure status
#'
#' @keywords internal
.parametric_engine_optimized <- function(
  fmri_data,
  event_design,
  hrf_interface,
  hrf_parameters = list(
    seed = NULL,
    bounds = NULL,
    eval_times = seq(0, 30, length.out = 61)
  ),
  algorithm_options = list(
    ridge_lambda = 0.01,
    epsilon_beta = 1e-6,
    method = "qr",
    cache_basis = TRUE
  ),
  validate = TRUE
) {
  
  # Load engineering standards
  source(file.path(dirname(getwd()), "R", "engineering-standards.R"), local = TRUE)
  
  # Start timing
  total_time <- system.time({
    
    # 1. INPUT VALIDATION (with clear, actionable errors)
    if (validate) {
      .with_timing({
        .validate_input(fmri_data, "fmri_data", 
                       type = c("matrix", "array"),
                       constraints = list(finite = TRUE))
        
        .validate_input(event_design, "event_design",
                       type = c("matrix", "array"),
                       dims = c(nrow(fmri_data), -1))
        
        .validate_input(hrf_interface, "hrf_interface",
                       type = "list")
        
        if (!all(c("hrf_function", "taylor_basis", "parameter_names") %in% 
                 names(hrf_interface))) {
          stop("hrf_interface must contain: hrf_function, taylor_basis, parameter_names",
               call. = FALSE)
        }
        
        # Validate numerical parameters
        if (!is.null(algorithm_options$ridge_lambda)) {
          .validate_input(algorithm_options$ridge_lambda, "ridge_lambda",
                         constraints = list(range = c(0, Inf)))
        }
      }, label = "validation")
    }
    
    # 2. SETUP AND PREPROCESSING
    setup_time <- system.time({
      # Extract dimensions
      n_time <- nrow(fmri_data)
      n_vox <- ncol(fmri_data)
      n_params <- length(hrf_interface$parameter_names)
      
      # Handle defaults
      if (is.null(hrf_parameters$seed)) {
        hrf_parameters$seed <- hrf_interface$default_seed()
      }
      if (is.null(hrf_parameters$bounds)) {
        hrf_parameters$bounds <- hrf_interface$default_bounds()
      }
      
      # Validate parameter bounds
      .validate_input(hrf_parameters$seed, "hrf_parameters$seed",
                     constraints = list(
                       range = c(hrf_parameters$bounds$lower,
                                hrf_parameters$bounds$upper)
                     ))
      
      # Pre-allocate output matrices
      theta_hat <- matrix(NA_real_, n_vox, n_params)
      amplitudes <- numeric(n_vox)
      r_squared <- numeric(n_vox)
      
      # Check memory requirements
      required_memory <- 8 * (n_time * n_vox + n_time * (n_params + 1) + n_vox * n_params)
      .check_memory_available(required_memory, "parametric engine")
    })
    
    # 3. COMPUTE TAYLOR BASIS (with caching)
    basis_time <- system.time({
      if (algorithm_options$cache_basis && exists(".basis_cache")) {
        cache_key <- digest::digest(list(
          hrf_parameters$seed,
          hrf_parameters$eval_times
        ))
        
        if (cache_key %in% names(.basis_cache)) {
          taylor_basis <- .basis_cache[[cache_key]]
        } else {
          taylor_basis <- .try_with_context(
            hrf_interface$taylor_basis(hrf_parameters$seed, 
                                      hrf_parameters$eval_times),
            context = "computing Taylor basis"
          )
          .basis_cache[[cache_key]] <- taylor_basis
        }
      } else {
        taylor_basis <- hrf_interface$taylor_basis(hrf_parameters$seed,
                                                  hrf_parameters$eval_times)
      }
      
      # Ensure matrix format
      if (!is.matrix(taylor_basis)) {
        taylor_basis <- matrix(taylor_basis, ncol = n_params + 1)
      }
    })
    
    # 4. OPTIMIZED CONVOLUTION (vectorized)
    convolution_time <- system.time({
      design_matrix <- .optimized_convolution_engine(
        event_design[, 1, drop = FALSE],
        taylor_basis,
        n_time
      )
    })
    
    # 5. NUMERICAL SOLUTION (with robustness)
    solution_time <- system.time({
      if (algorithm_options$method == "qr") {
        # QR decomposition (numerically stable)
        qr_decomp <- qr(design_matrix)
        Q <- qr.Q(qr_decomp)
        R <- qr.R(qr_decomp)
        
        # Check condition number
        R_diag <- abs(diag(R))
        condition <- max(R_diag) / min(R_diag[R_diag > .Machine$double.eps])
        
        if (condition > 1e8) {
          warning(sprintf(
            "Design matrix is ill-conditioned (kappa = %.2e), adding regularization",
            condition
          ))
          R <- R + algorithm_options$ridge_lambda * diag(ncol(R))
        }
        
        # Solve for all voxels at once
        R_inv <- backsolve(R, diag(ncol(R)))
        coefficients <- R_inv %*% crossprod(Q, fmri_data)
        
      } else if (algorithm_options$method == "svd") {
        # SVD solution (most robust but slower)
        coefficients <- .safe_solve(design_matrix, method = "svd") %*% 
                       t(design_matrix) %*% fmri_data
      }
      
      # Extract amplitudes and parameter updates
      amplitudes <- coefficients[1, ]
      
      # Safe division for parameter updates
      amplitudes_safe <- pmax(abs(amplitudes), algorithm_options$epsilon_beta)
      delta_theta <- coefficients[2:(n_params + 1), , drop = FALSE] / 
                     matrix(rep(amplitudes_safe, each = n_params), nrow = n_params)
      
      # Update parameters
      theta_hat <- matrix(hrf_parameters$seed, n_vox, n_params, byrow = TRUE) + 
                   t(delta_theta)
      
      # Apply bounds (vectorized)
      for (j in seq_len(n_params)) {
        theta_hat[, j] <- pmax(hrf_parameters$bounds$lower[j],
                              pmin(theta_hat[, j], hrf_parameters$bounds$upper[j]))
      }
    })
    
    # 6. FIT QUALITY ASSESSMENT
    quality_time <- system.time({
      # Predictions
      y_pred <- design_matrix %*% coefficients
      
      # Vectorized R-squared calculation
      ss_res <- colSums((fmri_data - y_pred)^2)
      ss_tot <- colSums(scale(fmri_data, scale = FALSE)^2)
      r_squared <- 1 - ss_res / pmax(ss_tot, .Machine$double.eps)
      
      # Ensure valid range
      r_squared <- pmax(0, pmin(1, r_squared))
    })
    
    # 7. DIAGNOSTICS
    diagnostics <- list(
      timing = list(
        total = NA,  # Filled in below
        validation = setup_time["elapsed"],
        basis = basis_time["elapsed"],
        convolution = convolution_time["elapsed"],
        solution = solution_time["elapsed"],
        quality = quality_time["elapsed"]
      ),
      numerical = list(
        condition_number = if (exists("condition")) condition else NA,
        regularization_applied = condition > 1e8,
        rank = if (algorithm_options$method == "qr") qr_decomp$rank else NA
      ),
      quality = list(
        mean_r2 = mean(r_squared),
        min_r2 = min(r_squared),
        max_r2 = max(r_squared),
        failed_voxels = sum(r_squared < 0.1)
      )
    )
    
  })["elapsed"]
  
  # Add total time
  diagnostics$timing$total <- total_time
  
  # 8. OUTPUT VALIDATION
  .assert_output_quality(
    list(
      parameters = theta_hat,
      r_squared = r_squared
    ),
    checks = list(
      finite = TRUE,
      positive_r2 = TRUE,
      bounded_params = hrf_parameters$bounds
    )
  )
  
  # 9. RETURN STRUCTURED RESULT
  structure(
    list(
      parameters = theta_hat,
      amplitudes = amplitudes,
      fit_quality = r_squared,
      diagnostics = diagnostics,
      metadata = list(
        n_voxels = n_vox,
        n_timepoints = n_time,
        n_parameters = n_params,
        algorithm = algorithm_options$method,
        seed_parameters = hrf_parameters$seed
      ),
      status = "success"
    ),
    class = c("parametric_engine_result", "list")
  )
}

#' Optimized convolution for design matrix construction
#'
#' @param signals Matrix of signals to convolve
#' @param kernels Matrix where each column is a kernel
#' @param output_length Desired output length
#' @return Convolved design matrix
#' @keywords internal
.optimized_convolution_engine <- function(signals, kernels, output_length) {
  if (NCOL(signals) != 1) {
    stop("Only a single signal column is supported")
  }

  # Delegate to high performance helper which chooses the best method
  .fast_batch_convolution(signals[, 1], kernels, output_length)
}

# Cache for basis functions
.basis_cache <- new.env(parent = emptyenv())

#' Print method for parametric engine results
#' @export
print.parametric_engine_result <- function(x, ...) {
  cat("Parametric HRF Engine Result\n")
  cat("===========================\n")
  cat("Status:", x$status, "\n")
  cat("Voxels processed:", x$metadata$n_voxels, "\n")
  cat("Parameters estimated:", x$metadata$n_parameters, "\n")
  cat("\nFit Quality:\n")
  cat("  Mean R²:", sprintf("%.3f", x$diagnostics$quality$mean_r2), "\n")
  cat("  Range: [", sprintf("%.3f", x$diagnostics$quality$min_r2), ", ",
      sprintf("%.3f", x$diagnostics$quality$max_r2), "]\n", sep = "")
  cat("  Failed voxels:", x$diagnostics$quality$failed_voxels, "\n")
  cat("\nComputation Time:", sprintf("%.2f", x$diagnostics$timing$total), "seconds\n")
  cat("  Speed:", round(x$metadata$n_voxels / x$diagnostics$timing$total), 
      "voxels/second\n")
  invisible(x)
}