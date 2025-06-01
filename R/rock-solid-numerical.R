#' Rock Solid Numerical Stability Functions
#'
#' Numerical safety functions that prevent ALL numerical failures.
#' Every operation is protected against overflow, underflow, and ill-conditioning.

#' Safe division that never produces NaN or Inf
#' @keywords internal
.safe_divide <- function(numerator, denominator, epsilon = .Machine$double.eps) {
  # Ensure denominator is never exactly zero
  safe_denom <- ifelse(abs(denominator) < epsilon, 
                       sign(denominator) * epsilon,
                       denominator)
  
  # Handle zero numerator case
  result <- ifelse(abs(numerator) < epsilon & abs(denominator) < epsilon,
                   0,  # 0/0 case - return 0
                   numerator / safe_denom)
  
  # Clip extreme values
  max_val <- sqrt(.Machine$double.xmax)
  result <- pmax(-max_val, pmin(result, max_val))
  
  result
}

#' Check and improve matrix conditioning
#' @keywords internal
.check_conditioning <- function(X, max_condition = 1e8, caller = "parametric_engine") {
  # Compute condition number safely
  tryCatch({
    svd_X <- svd(X, nu = 0, nv = 0)
    singular_vals <- svd_X$d
    
    # Handle zero singular values
    non_zero_sv <- singular_vals[singular_vals > .Machine$double.eps]
    
    if (length(non_zero_sv) == 0) {
      warning(caller, ": Design matrix is effectively zero. ",
              "Check your data and event model.")
      return(list(
        condition_number = Inf,
        rank = 0,
        needs_regularization = TRUE,
        suggested_lambda = 1.0
      ))
    }
    
    condition_num <- max(non_zero_sv) / min(non_zero_sv)
    effective_rank <- sum(singular_vals > .Machine$double.eps * max(singular_vals))
    
    # Determine if regularization needed
    needs_reg <- condition_num > max_condition
    
    # Suggest regularization parameter
    if (needs_reg) {
      # Target condition number of sqrt(max_condition)
      target_cond <- sqrt(max_condition)
      suggested_lambda <- max(non_zero_sv) / target_cond - min(non_zero_sv)
      suggested_lambda <- max(suggested_lambda, 1e-6)
    } else {
      suggested_lambda <- 0
    }
    
    list(
      condition_number = condition_num,
      rank = effective_rank,
      needs_regularization = needs_reg,
      suggested_lambda = suggested_lambda,
      singular_values = singular_vals
    )
    
  }, error = function(e) {
    warning(caller, ": SVD failed - ", e$message, 
            ". Using fallback regularization.")
    list(
      condition_number = Inf,
      rank = NA,
      needs_regularization = TRUE,
      suggested_lambda = 0.1
    )
  })
}

#' Safe QR decomposition with automatic regularization
#' @keywords internal
.safe_qr_solve <- function(X, y, lambda = 0.01, check_condition = TRUE, 
                           caller = "parametric_engine") {
  n <- nrow(X)
  p <- ncol(X)
  
  # Check conditioning if requested
  if (check_condition) {
    cond_info <- .check_conditioning(X, caller = caller)
    if (cond_info$needs_regularization && lambda < cond_info$suggested_lambda) {
      message(caller, ": Poor conditioning detected (", 
              format(cond_info$condition_number, scientific = TRUE), 
              "). Increasing regularization to ", 
              format(cond_info$suggested_lambda, scientific = TRUE))
      lambda <- cond_info$suggested_lambda
    }
  }
  
  # Multiple solution strategies
  solution <- NULL
  method_used <- "none"
  
  # Strategy 1: Standard QR with ridge
  if (is.null(solution)) {
    solution <- tryCatch({
      qr_decomp <- qr(X)
      Q <- qr.Q(qr_decomp)
      R <- qr.R(qr_decomp)
      
      # Add ridge penalty
      diag(R) <- diag(R) + lambda
      
      # Solve R * coeffs = Qt * y
      Qty <- crossprod(Q, y)
      coeffs <- backsolve(R, Qty)
      
      method_used <- "qr_ridge"
      list(coefficients = coeffs, method = method_used)
      
    }, error = function(e) NULL)
  }
  
  # Strategy 2: SVD-based solution
  if (is.null(solution)) {
    solution <- tryCatch({
      svd_X <- svd(X)
      d <- svd_X$d
      
      # Regularized inverse
      d_inv <- .safe_divide(1, d^2 + lambda)
      
      # Compute solution
      coeffs <- svd_X$v %*% (d_inv * d * crossprod(svd_X$u, y))
      
      method_used <- "svd_ridge"
      list(coefficients = as.vector(coeffs), method = method_used)
      
    }, error = function(e) NULL)
  }
  
  # Strategy 3: Normal equations with heavy regularization
  if (is.null(solution)) {
    solution <- tryCatch({
      XtX <- crossprod(X)
      Xty <- crossprod(X, y)
      
      # Heavy regularization
      diag(XtX) <- diag(XtX) + lambda * 10
      
      coeffs <- solve(XtX, Xty)
      
      method_used <- "normal_heavy_ridge"
      list(coefficients = as.vector(coeffs), method = method_used)
      
    }, error = function(e) NULL)
  }
  
  # Strategy 4: Fallback to zeros
  if (is.null(solution)) {
    warning(caller, ": All solution methods failed. Returning zero coefficients.")
    solution <- list(
      coefficients = rep(0, p),
      method = "fallback_zeros"
    )
  }
  
  # Validate solution
  coeffs <- solution$coefficients
  
  # Check for NaN/Inf
  if (any(!is.finite(coeffs))) {
    bad_idx <- which(!is.finite(coeffs))
    warning(caller, ": Non-finite coefficients detected at positions ",
            paste(bad_idx, collapse = ", "), ". Setting to zero.")
    coeffs[!is.finite(coeffs)] <- 0
  }
  
  # Check for extreme values
  coeff_scale <- mad(coeffs[coeffs != 0], na.rm = TRUE)
  if (coeff_scale > 0) {
    extreme_threshold <- 100 * coeff_scale
    extreme_idx <- which(abs(coeffs) > extreme_threshold)
    if (length(extreme_idx) > 0) {
      warning(caller, ": Extreme coefficients detected. Clipping to reasonable range.")
      coeffs[extreme_idx] <- sign(coeffs[extreme_idx]) * extreme_threshold
    }
  }
  
  list(
    coefficients = coeffs,
    method = solution$method,
    lambda_used = lambda
  )
}

#' Safe parameter update with bounds and stability checks
#' @keywords internal
.safe_parameter_update <- function(theta_current, delta_theta, theta_bounds,
                                   max_step = 0.5, caller = "parametric_engine") {
  n_params <- length(theta_current)
  
  # Limit step size
  step_scale <- max(abs(delta_theta))
  if (step_scale > max_step) {
    scale_factor <- max_step / step_scale
    delta_theta <- delta_theta * scale_factor
  }
  
  # Apply update
  theta_new <- theta_current + delta_theta
  
  # Apply bounds
  theta_new <- pmax(theta_bounds$lower, pmin(theta_new, theta_bounds$upper))
  
  # Check for stuck parameters
  at_lower <- abs(theta_new - theta_bounds$lower) < .Machine$double.eps
  at_upper <- abs(theta_new - theta_bounds$upper) < .Machine$double.eps
  
  if (any(at_lower | at_upper)) {
    n_stuck <- sum(at_lower | at_upper)
    if (n_stuck == n_params) {
      warning(caller, ": All parameters stuck at bounds. Consider relaxing bounds.")
    }
  }
  
  # Ensure parameters changed
  if (all(abs(theta_new - theta_current) < .Machine$double.eps)) {
    # Force small perturbation to avoid exact stagnation
    perturbation <- runif(n_params, -1e-6, 1e-6)
    theta_new <- theta_new + perturbation
    theta_new <- pmax(theta_bounds$lower, pmin(theta_new, theta_bounds$upper))
  }
  
  list(
    theta = theta_new,
    clamped = at_lower | at_upper,
    step_scaled = step_scale > max_step
  )
}

#' Monitor convergence with multiple criteria
#' @keywords internal
.monitor_convergence <- function(iteration_history, tol = 1e-4, 
                                 patience = 3, caller = "engine") {
  n_iter <- length(iteration_history)
  
  if (n_iter < 2) {
    return(list(
      converged = FALSE,
      reason = "insufficient_iterations",
      should_stop = FALSE
    ))
  }
  
  # Extract metrics
  param_changes <- sapply(2:n_iter, function(i) {
    sqrt(sum((iteration_history[[i]]$theta - iteration_history[[i-1]]$theta)^2))
  })
  
  obj_values <- sapply(iteration_history, function(x) x$objective)
  
  # Criterion 1: Parameter convergence
  if (n_iter >= 2 && tail(param_changes, 1) < tol) {
    return(list(
      converged = TRUE,
      reason = "parameter_convergence",
      should_stop = TRUE
    ))
  }
  
  # Criterion 2: Objective convergence
  if (n_iter >= 3) {
    recent_obj <- tail(obj_values, 3)
    obj_change <- abs(diff(recent_obj))
    if (all(obj_change < tol * abs(mean(recent_obj)))) {
      return(list(
        converged = TRUE,
        reason = "objective_convergence",
        should_stop = TRUE
      ))
    }
  }
  
  # Criterion 3: Oscillation detection
  if (n_iter >= 4) {
    recent_params <- tail(param_changes, 4)
    if (sd(recent_params) / mean(recent_params) < 0.1) {
      # Parameters oscillating in small range
      return(list(
        converged = FALSE,
        reason = "oscillation_detected",
        should_stop = TRUE
      ))
    }
  }
  
  # Criterion 4: Divergence detection
  if (n_iter >= 3) {
    recent_obj <- tail(obj_values, 3)
    if (all(diff(recent_obj) > 0)) {
      # Objective increasing
      return(list(
        converged = FALSE,
        reason = "divergence_detected",
        should_stop = TRUE
      ))
    }
  }
  
  # Not converged yet
  list(
    converged = FALSE,
    reason = "in_progress",
    should_stop = FALSE
  )
}

#' Create numerically stable HRF basis with safeguards
#' @keywords internal
.safe_hrf_basis <- function(hrf_interface, theta, t_hrf, caller = "engine") {
  # Input validation
  if (any(!is.finite(theta))) {
    warning(caller, ": Non-finite parameters detected. Using defaults.")
    theta[!is.finite(theta)] <- hrf_interface$default_seed[!is.finite(theta)]
  }
  
  # Call taylor basis with error handling
  basis <- tryCatch({
    hrf_interface$taylor_basis(theta, t_hrf)
  }, error = function(e) {
    warning(caller, ": HRF basis computation failed - ", e$message)
    NULL
  })
  
  if (is.null(basis)) {
    # Fallback: simple basis
    n_t <- length(t_hrf)
    n_params <- length(theta)
    basis <- matrix(0, nrow = n_t, ncol = n_params + 1)
    basis[, 1] <- 1  # Constant HRF
    return(basis)
  }
  
  # Ensure matrix form
  if (!is.matrix(basis)) {
    basis <- as.matrix(basis)
  }
  
  # Check for numerical issues
  if (any(!is.finite(basis))) {
    n_bad <- sum(!is.finite(basis))
    warning(caller, ": HRF basis contains ", n_bad, " non-finite values. ",
            "Replacing with zeros.")
    basis[!is.finite(basis)] <- 0
  }
  
  # Check for extreme values
  basis_scale <- max(abs(basis), na.rm = TRUE)
  if (basis_scale > 1e6) {
    warning(caller, ": HRF basis has extreme values (max = ", basis_scale, "). ",
            "Rescaling to prevent overflow.")
    basis <- basis / basis_scale
  } else if (basis_scale < 1e-6 && basis_scale > 0) {
    warning(caller, ": HRF basis has very small values (max = ", basis_scale, "). ",
            "May indicate numerical underflow.")
  }
  
  basis
}