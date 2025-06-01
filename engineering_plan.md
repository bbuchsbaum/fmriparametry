# Engineering Excellence Plan for fmriparametric

## Mission Statement

Transform `fmriparametric` from "middling" to **exemplary** through rigorous engineering practices, algorithmic excellence, and uncompromising code quality. We will demonstrate that computational neuroscience deserves the same engineering standards as mission-critical systems.

## Current State Assessment

### Identified Weaknesses
1. **Inconsistent interfaces** - Mixed parameter naming, varying return structures
2. **Code duplication** - Similar patterns repeated across modules
3. **Unclear contracts** - Functions lack clear pre/post conditions
4. **Performance gaps** - Unoptimized matrix operations, redundant computations
5. **Error handling** - Mix of strategies without clear hierarchy
6. **Documentation debt** - Implementation details poorly explained

### Our Response
We reject "middling." We will systematically address each weakness with engineering rigor.

---

## I. INTERFACE EXCELLENCE

### 1.1 Consistent Parameter Naming Convention

**Current State**: Inconsistent (`Y_proj`, `fmri_data`, `Y`)
**Target State**: Unified, predictable naming

```r
# STANDARD: noun_modifier pattern
# data
fmri_data      # Never: Y, Y_proj, data
event_design   # Never: S, S_target_proj, event_model
mask_data      # Never: mask, brain_mask

# parameters  
hrf_parameters # Never: theta, params, theta_hat
hrf_bounds     # Never: theta_bounds, bounds
hrf_seed       # Never: theta_seed, seed

# options
ridge_lambda   # Never: lambda_ridge, lambda
n_iterations   # Never: recenter_global_passes, iterations
convergence_tol # Never: recenter_epsilon, epsilon
```

### 1.2 Function Signature Standards

**Every function follows this pattern:**
```r
.private_function_name <- function(
  # Required data arguments (no defaults)
  data_arg1,
  data_arg2,
  
  # Required configuration (no defaults)
  config_arg1,
  
  # Optional parameters (with defaults)
  param1 = default1,
  param2 = default2,
  
  # Control flags (always with defaults)
  verbose = FALSE,
  validate = TRUE
) {
  # Contract validation
  if (validate) {
    .validate_inputs(...)
  }
  
  # Core logic
  result <- ...
  
  # Post-condition check
  .assert_valid_output(result)
  
  return(result)
}
```

### 1.3 Return Value Contracts

**Every function returns a documented structure:**

```r
# Level 1: Simple returns (internal functions)
# Always return named lists with fixed structure
list(
  result = ...,      # Primary result
  diagnostic = ...,  # Optional diagnostic info
  status = "success" # Always include status
)

# Level 2: S3 objects (user-facing)
# Always use constructor functions
new_parametric_fit <- function(parameters, metadata) {
  structure(
    list(
      # Data slots (never NULL)
      parameters = parameters,
      amplitudes = amplitudes,
      
      # Metadata (never NULL)
      metadata = list(
        version = "3.0.0",
        timestamp = Sys.time(),
        ...
      ),
      
      # Optional slots (can be NULL but must exist)
      diagnostics = NULL,
      convergence = NULL
    ),
    class = c("parametric_hrf_fit", "list")
  )
}
```

---

## II. ALGORITHMIC IMPECCABILITY

### 2.1 Numerical Precision Standards

**Every numerical operation must be robust:**

```r
# NEVER this:
delta <- (f(x + h) - f(x - h)) / (2 * h)

# ALWAYS this:
.numerical_derivative <- function(f, x, h = NULL, method = "central") {
  # Automatic step size selection
  if (is.null(h)) {
    h <- .optimal_step_size(x, method)
  }
  
  # Guard against catastrophic cancellation
  x_plus <- x + h
  x_minus <- x - h
  h_actual <- x_plus - x_minus  # Account for floating point
  
  # Compute with error estimation
  if (method == "central") {
    f_plus <- f(x_plus)
    f_minus <- f(x_minus)
    
    # Check for numerical issues
    if (!is.finite(f_plus) || !is.finite(f_minus)) {
      # Fall back to one-sided
      return(.numerical_derivative(f, x, h, method = "forward"))
    }
    
    derivative <- (f_plus - f_minus) / h_actual
    
    # Richardson extrapolation for higher accuracy
    if (isTRUE(getOption("fmriparametric.high_precision"))) {
      derivative <- .richardson_extrapolation(f, x, h, derivative)
    }
  }
  
  list(
    value = derivative,
    error_estimate = .truncation_error(h, method),
    method_used = method
  )
}
```

### 2.2 Matrix Operation Standards

**Eliminate redundant computations:**

```r
# NEVER this:
for (v in 1:n_voxels) {
  X_v <- construct_design_matrix(...)
  beta_v <- solve(t(X_v) %*% X_v) %*% t(X_v) %*% Y[, v]
}

# ALWAYS this:
.batch_regression <- function(Y, X, method = "qr", stabilize = TRUE) {
  # Pre-compute once
  if (method == "qr") {
    qr_X <- qr(X)
    Q <- qr.Q(qr_X)
    R <- qr.R(qr_X)
    
    if (stabilize) {
      # Condition number check
      condition <- kappa(R)
      if (condition > 1e8) {
        # Add adaptive regularization
        lambda <- .adaptive_ridge_penalty(R)
        R <- R + lambda * diag(ncol(R))
      }
    }
    
    # Solve for all voxels at once
    R_inv <- backsolve(R, diag(ncol(R)))
    coefficients <- R_inv %*% crossprod(Q, Y)
    
    # Compute diagnostics efficiently
    fitted <- Q %*% (crossprod(Q, Y))
    residuals <- Y - fitted
    
    # Vectorized R-squared
    ss_res <- colSums(residuals^2)
    ss_tot <- colSums(scale(Y, scale = FALSE)^2)
    r_squared <- 1 - ss_res / ss_tot
  }
  
  list(
    coefficients = coefficients,
    fitted = fitted,
    residuals = residuals,
    r_squared = r_squared,
    condition_number = condition,
    method = method
  )
}
```

### 2.3 Convolution Optimization

**Current**: Loop-based, redundant FFT operations
**Target**: Vectorized, cached FFT

```r
.optimized_convolution <- function(signals, kernel, method = "auto") {
  n_signals <- ncol(signals)
  n_time <- nrow(signals)
  kernel_length <- length(kernel)
  
  # Method selection based on problem size
  if (method == "auto") {
    # Empirically determined crossover points
    if (n_time < 100 || kernel_length < 10) {
      method <- "direct"
    } else if (n_signals > 100) {
      method <- "fft_batch"
    } else {
      method <- "fft"
    }
  }
  
  if (method == "fft_batch") {
    # Pad for FFT
    n_fft <- nextn(n_time + kernel_length - 1, factors = 2)
    
    # Pre-compute kernel FFT once
    kernel_fft <- fft(c(kernel, rep(0, n_fft - kernel_length)))
    
    # Batch FFT for all signals
    signals_padded <- rbind(signals, matrix(0, n_fft - n_time, n_signals))
    signals_fft <- mvfft(signals_padded)
    
    # Multiply in frequency domain
    conv_fft <- signals_fft * kernel_fft
    
    # Inverse FFT and trim
    conv_full <- Re(mvfft(conv_fft, inverse = TRUE) / n_fft)
    conv_trimmed <- conv_full[1:n_time, , drop = FALSE]
    
    return(conv_trimmed)
  }
  
  # Fallback methods...
}
```

---

## III. ROBUSTNESS ENGINEERING

### 3.1 Input Validation Hierarchy

```r
# Level 1: Type validation (fast)
.validate_type <- function(x, expected_type, arg_name) {
  if (!inherits(x, expected_type)) {
    stop(sprintf(
      "Argument '%s' must be of type '%s', got '%s'",
      arg_name, expected_type, class(x)[1]
    ), call. = FALSE)
  }
}

# Level 2: Structure validation (medium)
.validate_structure <- function(x, requirements, arg_name) {
  # Dimensions
  if (!is.null(requirements$dims)) {
    if (!all(dim(x) == requirements$dims)) {
      stop(sprintf(
        "Argument '%s' must have dimensions %s, got %s",
        arg_name, 
        paste(requirements$dims, collapse = " x "),
        paste(dim(x), collapse = " x ")
      ), call. = FALSE)
    }
  }
  
  # Additional checks...
}

# Level 3: Content validation (slow, optional)
.validate_content <- function(x, constraints, arg_name) {
  # Numerical constraints
  if (!is.null(constraints$range)) {
    out_of_range <- x < constraints$range[1] | x > constraints$range[2]
    if (any(out_of_range, na.rm = TRUE)) {
      stop(sprintf(
        "Argument '%s' contains values outside range [%g, %g]",
        arg_name, constraints$range[1], constraints$range[2]
      ), call. = FALSE)
    }
  }
  
  # Statistical constraints
  if (!is.null(constraints$properties)) {
    for (prop in constraints$properties) {
      if (!.check_property(x, prop)) {
        stop(sprintf(
          "Argument '%s' violates constraint: %s",
          arg_name, prop
        ), call. = FALSE)
      }
    }
  }
}
```

### 3.2 Defensive Programming Standards

```r
# Every function that could fail gets a safety wrapper
.safe_operation <- function(operation, fallback = NULL, context = "") {
  tryCatch(
    {
      result <- operation()
      # Verify result is valid
      if (is.null(result) || any(!is.finite(result))) {
        stop("Invalid result from operation")
      }
      result
    },
    error = function(e) {
      # Log detailed context
      .log_error(
        operation = context,
        error = e$message,
        traceback = sys.calls()
      )
      
      # Attempt fallback
      if (!is.null(fallback)) {
        tryCatch(
          fallback(),
          error = function(e2) {
            stop(sprintf(
              "Primary operation failed: %s\nFallback also failed: %s",
              e$message, e2$message
            ), call. = FALSE)
          }
        )
      } else {
        stop(sprintf(
          "Operation '%s' failed: %s",
          context, e$message
        ), call. = FALSE)
      }
    }
  )
}
```

---

## IV. PERFORMANCE ENGINEERING

### 4.1 Profiling-Driven Optimization

```r
# Built-in profiling for critical functions
.with_profiling <- function(expr, label) {
  if (getOption("fmriparametric.profile", FALSE)) {
    start_time <- Sys.time()
    start_mem <- gc(reset = TRUE)
    
    result <- force(expr)
    
    end_mem <- gc()
    end_time <- Sys.time()
    
    .record_performance(
      label = label,
      time = as.numeric(end_time - start_time, units = "secs"),
      memory = sum(end_mem[, 2] - start_mem[, 2])
    )
    
    result
  } else {
    expr
  }
}
```

### 4.2 Memory-Efficient Algorithms

```r
# Chunked processing for large datasets
.chunked_apply <- function(data, fun, chunk_size = NULL, 
                          combine = rbind, progress = TRUE) {
  n_total <- nrow(data)
  
  # Intelligent chunk size selection
  if (is.null(chunk_size)) {
    available_memory <- .get_available_memory()
    data_size <- object.size(data)
    overhead_factor <- 3  # Conservative estimate
    chunk_size <- floor(n_total * available_memory / (data_size * overhead_factor))
    chunk_size <- max(100, min(chunk_size, 10000))
  }
  
  n_chunks <- ceiling(n_total / chunk_size)
  results <- vector("list", n_chunks)
  
  if (progress) {
    pb <- .create_progress_bar(n_chunks)
  }
  
  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_total)
    
    chunk <- data[start_idx:end_idx, , drop = FALSE]
    
    # Process chunk with automatic cleanup
    results[[i]] <- .with_cleanup(
      fun(chunk),
      cleanup = function() {
        if (i %% 10 == 0) gc()  # Periodic garbage collection
      }
    )
    
    if (progress) .update_progress_bar(pb, i)
  }
  
  # Efficient combination
  do.call(combine, results)
}
```

### 4.3 Cache-Aware Algorithms

```r
# Cache frequently computed values
.cached_hrf_basis <- local({
  cache <- new.env(parent = emptyenv())
  
  function(theta, t, cache_key = NULL) {
    # Generate cache key if not provided
    if (is.null(cache_key)) {
      cache_key <- digest::digest(list(theta, t))
    }
    
    # Check cache
    if (exists(cache_key, envir = cache)) {
      .increment_cache_hits()
      return(get(cache_key, envir = cache))
    }
    
    # Compute
    result <- .compute_hrf_basis(theta, t)
    
    # Cache with size limit
    cache_size <- length(ls(cache))
    if (cache_size > getOption("fmriparametric.cache_size", 100)) {
      # LRU eviction
      .evict_lru_cache_entry(cache)
    }
    
    assign(cache_key, result, envir = cache)
    result
  }
})
```

---

## V. CODE QUALITY STANDARDS

### 5.1 Function Complexity Limits

**Rule**: No function exceeds cyclomatic complexity of 10

```r
# BAD: Complex monolithic function
estimate_hrf <- function(...) {
  # 200 lines of nested if-else and loops
}

# GOOD: Decomposed into logical units
estimate_hrf <- function(fmri_data, event_design, options = list()) {
  # Validate
  validated_inputs <- .validate_hrf_inputs(fmri_data, event_design, options)
  
  # Prepare
  prepared_data <- .prepare_hrf_data(validated_inputs)
  
  # Estimate
  initial_estimate <- .compute_initial_estimate(prepared_data)
  
  # Refine
  refined_estimate <- .refine_estimate(initial_estimate, prepared_data, options)
  
  # Package results
  .create_hrf_result(refined_estimate, metadata = prepared_data$metadata)
}
```

### 5.2 Documentation Standards

```r
#' Brief description (one line)
#'
#' Detailed description explaining the algorithm, assumptions,
#' and theoretical background. Include references.
#'
#' @section Algorithm:
#' Mathematical description of the algorithm:
#' \deqn{h(t; \theta) = \sum_{k=0}^{K} \frac{1}{k!} \frac{\partial^k h}{\partial \theta^k}|_{\theta_0} (\theta - \theta_0)^k}
#'
#' @section Performance:
#' Time complexity: O(n * p^2) where n = timepoints, p = parameters
#' Space complexity: O(n * v) where v = voxels
#'
#' @param x Description including type, constraints, and units
#' @param validate (logical) Whether to validate inputs. Set FALSE only
#'   when calling from validated parent function.
#'
#' @return Description of return value structure
#'
#' @examples
#' # Example with expected output
#' result <- function_name(x = valid_input)
#' stopifnot(all.equal(result$value, expected_value))
#'
#' @references
#' Author, A. (2024). Title. Journal, 1(1), 1-10.
```

### 5.3 Testing Standards

```r
# Every function gets:
# 1. Unit tests (test-function_name.R)
# 2. Property tests
# 3. Edge case tests
# 4. Performance regression tests

test_that("function satisfies mathematical properties", {
  # Property: Taylor approximation converges to true function
  for (n_terms in 1:4) {
    theta0 <- c(5, 2, 0.3)
    theta_test <- theta0 + rnorm(3, sd = 0.1)
    
    approx <- taylor_approximation(theta_test, theta0, n_terms)
    true_val <- true_function(theta_test)
    
    error <- abs(approx - true_val)
    
    # Error should decrease with more terms
    if (n_terms > 1) {
      expect_true(error < previous_error)
    }
    previous_error <- error
  }
})

test_that("function handles edge cases correctly", {
  # Systematic edge case testing
  edge_cases <- list(
    "empty_data" = matrix(nrow = 0, ncol = 0),
    "single_point" = matrix(1),
    "constant_data" = matrix(1, 100, 10),
    "extreme_values" = matrix(c(-Inf, Inf, NA, NaN), 2, 2)
  )
  
  for (case_name in names(edge_cases)) {
    result <- safely_call(function_name, edge_cases[[case_name]])
    expect_true(
      !is.null(result),
      info = sprintf("Failed on edge case: %s", case_name)
    )
  }
})
```

---

## VI. IMPLEMENTATION TIMELINE

### Phase 1: Interface Standardization (Week 1)
- [ ] Audit all function signatures
- [ ] Implement consistent naming convention
- [ ] Standardize return structures
- [ ] Update all internal callers

### Phase 2: Algorithm Refinement (Week 2)
- [ ] Implement numerical safeguards
- [ ] Optimize matrix operations
- [ ] Add caching layer
- [ ] Benchmark improvements

### Phase 3: Robustness Hardening (Week 3)
- [ ] Implement validation hierarchy
- [ ] Add defensive wrappers
- [ ] Comprehensive error handling
- [ ] Stress testing

### Phase 4: Performance Optimization (Week 4)
- [ ] Profile critical paths
- [ ] Implement chunked processing
- [ ] Optimize memory usage
- [ ] Parallel processing refinement

### Phase 5: Quality Assurance (Week 5)
- [ ] Complete documentation
- [ ] Achieve 95% test coverage
- [ ] Performance regression suite
- [ ] Code review and refactoring

---

## VII. SUCCESS METRICS

### Code Quality Metrics
- **Cyclomatic complexity**: All functions ≤ 10
- **Code coverage**: ≥ 95%
- **Documentation coverage**: 100%
- **Static analysis**: Zero warnings (lintr, goodpractice)

### Performance Metrics
- **Speed**: 2x faster than current implementation
- **Memory**: 50% reduction in peak usage
- **Scalability**: Linear scaling to 1M voxels

### Algorithmic Metrics
- **Numerical accuracy**: Machine precision where applicable
- **Convergence**: Guaranteed for well-posed problems
- **Stability**: Condition number monitoring and adaptation

### User Experience Metrics
- **API consistency**: Zero breaking changes after v3.1
- **Error messages**: 100% actionable
- **Examples**: Every function has runnable example

---

## VIII. ENGINEERING PRINCIPLES

1. **No Compromise on Quality**: Every line of code meets our standards
2. **Algorithmic Correctness First**: Performance never at the cost of correctness
3. **Fail Fast and Loud**: Bad inputs detected immediately with clear messages
4. **Reproducibility**: Same inputs ALWAYS produce same outputs
5. **Defensive by Default**: Trust nothing, verify everything

---

## CONCLUSION

This plan transforms `fmriparametric` from "middling" to **exceptional**. We will demonstrate that scientific computing deserves the same engineering excellence as any critical system.

The days of accepting "good enough" in scientific software are over. We are engineers. We will engineer accordingly.

**Status**: READY TO EXECUTE
**Commitment**: TOTAL
**Outcome**: EXCELLENCE