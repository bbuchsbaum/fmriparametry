#' Gauss-Newton refinement for hard voxels
#'
#' Implements full nonlinear Gauss-Newton optimization for the most challenging
#' voxels that don't respond well to Taylor approximation methods. This is the
#' most computationally intensive refinement but can recover good fits for
#' difficult cases.
#'
#' @param theta_hat_voxel Matrix of current parameter estimates (voxels x parameters)
#' @param r2_voxel Numeric vector of current R-squared values
#' @param Y_proj Numeric matrix of projected BOLD data (timepoints x voxels)
#' @param S_target_proj Numeric matrix of projected stimulus design
#' @param scan_times Numeric vector of scan acquisition times
#' @param hrf_eval_times Numeric vector of HRF evaluation time points
#' @param hrf_interface List with HRF interface functions
#' @param theta_bounds List with elements `lower` and `upper`
#' @param queue_labels Character vector of refinement queue assignments
#' @param max_iter_gn Maximum iterations for Gauss-Newton
#' @param tol_gn Convergence tolerance
#' @param lambda_ridge Ridge penalty
#' @param step_size Initial step size for line search
#' @param verbose Logical whether to print progress
#'
#' @return List with updated estimates and convergence information
#' @keywords internal
.gauss_newton_refinement <- function(
  theta_hat_voxel,
  r2_voxel,
  Y_proj,
  S_target_proj,
  scan_times,
  hrf_eval_times,
  hrf_interface,
  theta_bounds,
  queue_labels,
  max_iter_gn = 5,
  tol_gn = 1e-4,
  lambda_ridge = 0.01,
  step_size = 1.0,
  verbose = FALSE
) {
  n_vox <- nrow(theta_hat_voxel)
  n_params <- ncol(theta_hat_voxel)
  n_time <- nrow(Y_proj)
  
  # Identify hard voxels
  idx_hard <- which(queue_labels == "hard_GN")
  n_hard <- length(idx_hard)
  
  if (n_hard == 0) {
    if (verbose) cat("No hard voxels to refine\n")
    return(list(
      theta_hat = theta_hat_voxel,
      r2 = r2_voxel,
      n_refined = 0,
      n_converged = 0,
      n_improved = 0,
      convergence_status = character(0)
    ))
  }
  
  if (verbose) {
    cat("\nGauss-Newton refinement for", n_hard, "hard voxels\n")
    cat("  Max iterations:", max_iter_gn, "\n")
    cat("  Convergence tolerance:", tol_gn, "\n")
  }
  
  # Store original values
  theta_hat_orig <- theta_hat_voxel
  r2_orig <- r2_voxel
  
  # Track convergence
  convergence_status <- character(n_hard)
  n_converged <- 0
  n_improved <- 0
  iteration_counts <- numeric(n_hard)
  
  # Process each hard voxel
  for (i in seq_along(idx_hard)) {
    v <- idx_hard[i]
    y_v <- Y_proj[, v]
    
    # Initialize with current estimate
    if (v > nrow(theta_hat_voxel)) {
      warning(sprintf("Gauss-Newton: Voxel index %d exceeds matrix rows %d", v, nrow(theta_hat_voxel)))
      convergence_status[i] <- "invalid_index"
      iteration_counts[i] <- 0
      next
    }
    
    theta_current <- theta_hat_voxel[v, , drop = TRUE]
    
    # Validate extraction
    if (length(theta_current) == 0 || length(theta_current) != n_params) {
      warning(sprintf("Gauss-Newton: Invalid theta extraction for voxel %d (got length %d, expected %d)",
                     v, length(theta_current), n_params))
      cat("Debug - theta_hat_voxel dimensions:", dim(theta_hat_voxel), "\n")
      cat("Debug - voxel index:", v, "\n")
      convergence_status[i] <- "invalid_input"
      iteration_counts[i] <- 0
      next
    }
    
    # Ensure theta_current is a numeric vector
    if (is.matrix(theta_current)) {
      theta_current <- as.numeric(theta_current)
    }
    theta_best <- theta_current
    r2_best <- r2_voxel[v]
    
    # Calculate initial objective
    obj_current <- tryCatch({
      .calculate_objective_gn(
        theta_current, y_v, S_target_proj, hrf_eval_times,
        hrf_interface, n_time
      )
    }, error = function(e) {
      warning(sprintf("GN: Error in initial objective for voxel %d: %s", v, e$message))
      cat("Debug - theta_current:", theta_current, "\n")
      cat("Debug - y_v length:", length(y_v), "\n")
      cat("Debug - S_target_proj dim:", dim(S_target_proj), "\n")
      Inf
    })

    if (is.infinite(obj_current)) {
      convergence_status[i] <- "singular_system"
      iteration_counts[i] <- 0
      next
    }
    
    converged <- FALSE
    iter <- 0
    singular <- FALSE
    
    # Gauss-Newton iterations
    for (iter in seq_len(max_iter_gn)) {
      # Get Jacobian and residuals at current point
      jacob_info <- .get_jacobian_and_residuals(
        theta_current, y_v, S_target_proj, hrf_eval_times,
        hrf_interface, n_time
      )
      
      if (is.null(jacob_info)) {
        convergence_status[i] <- "singular_system"
        break
      }
      
      J <- jacob_info$jacobian
      residuals <- jacob_info$residuals
      
      # Gauss-Newton direction: solve (J'J + lambda*I) * delta = -J' * r
      JtJ <- crossprod(J)
      Jtr <- crossprod(J, residuals)
      
      # Add ridge penalty for stability
      JtJ_ridge <- JtJ + lambda_ridge * diag(n_params)
      
      # Solve for update direction
      delta <- tryCatch({
        -solve(JtJ_ridge, Jtr)
      }, error = function(e) {
        NULL
      })
      
      if (is.null(delta) || length(delta) == 0) {
        convergence_status[i] <- "singular_system"
        break
      }
      
      # Ensure delta is numeric vector of correct length
      if (is.matrix(delta)) {
        delta <- as.numeric(delta)
      }
      
      if (length(delta) != n_params) {
        warning(sprintf("GN: Invalid delta length %d (expected %d)", length(delta), n_params))
        convergence_status[i] <- "numerical_error"
        break
      }
      
      # Line search to find appropriate step size
      alpha <- step_size
      theta_new <- theta_current
      obj_new <- obj_current
      
      for (ls_iter in 1:10) {
        # Proposed update
        # Ensure both are numeric vectors
        if (!is.numeric(theta_current)) theta_current <- as.numeric(theta_current)
        if (!is.numeric(delta)) delta <- as.numeric(delta)
        
        theta_proposal <- theta_current + alpha * delta
        
        # Debug check
        if (length(theta_proposal) == 0) {
          warning(sprintf("GN Debug: theta_proposal is empty! theta_current length=%d, delta length=%d, alpha=%f",
                         length(theta_current), length(delta), alpha))
          cat("theta_current:", theta_current, "\n")
          cat("delta:", delta, "\n")
          cat("class(theta_current):", class(theta_current), "\n")
          cat("class(delta):", class(delta), "\n")
          break
        }
        
        # Validate theta_proposal 
        if (length(theta_proposal) != n_params) {
          warning(sprintf("Gauss-Newton: Invalid theta_proposal length %d (expected %d) in line search",
                         length(theta_proposal), n_params))
          convergence_status[i] <- "numerical_error"
          singular <- TRUE
          break
        }
        
        # Apply bounds
        theta_proposal_bounded <- pmax(theta_bounds$lower, 
                                       pmin(theta_proposal, theta_bounds$upper))
        
        # Check if bounds application caused issues
        if (length(theta_proposal_bounded) != length(theta_proposal)) {
          warning("Bounds application changed vector length!")
          cat("Before bounds:", theta_proposal, "\n")
          cat("After bounds:", theta_proposal_bounded, "\n")
          cat("Bounds lower:", theta_bounds$lower, "\n")
          cat("Bounds upper:", theta_bounds$upper, "\n")
        }
        
        theta_proposal <- theta_proposal_bounded
        
        # Evaluate objective
        obj_proposal <- tryCatch({
          .calculate_objective_gn(
            theta_proposal, y_v, S_target_proj, hrf_eval_times,
            hrf_interface, n_time
          )
        }, error = function(e) {
          warning(sprintf("GN: Error in line search objective: %s", e$message))
          cat("Debug - theta_proposal:", theta_proposal, "\n")
          cat("Debug - theta_current:", theta_current, "\n")
          cat("Debug - delta:", delta, "\n")
          cat("Debug - alpha:", alpha, "\n")
          Inf
        })

        if (is.infinite(obj_proposal)) {
          convergence_status[i] <- "singular_system"
          singular <- TRUE
          break
        }
        
        # Accept if improved
        if (obj_proposal < obj_current) {
          theta_new <- theta_proposal
          obj_new <- obj_proposal
          break
        }
        
        # Reduce step size
        alpha <- alpha * 0.5
        
        if (alpha < 1e-6) {
          break
        }
      }

      if (singular) {
        break
      }
      
      # Check if line search produced an update
      if (all(theta_new == theta_current)) {
        if (alpha < 1e-6) {
          convergence_status[i] <- "line_search_failed"
          break
        } else {
          next
        }
      }

      # Check convergence
      param_change <- sqrt(sum((theta_new - theta_current)^2))
      obj_change <- abs(obj_new - obj_current) / (abs(obj_current) + 1e-10)
      
      if (param_change < tol_gn || obj_change < tol_gn) {
        converged <- TRUE
        convergence_status[i] <- "converged"
        n_converged <- n_converged + 1
        break
      }
      
      # Update for next iteration
      theta_current <- theta_new
      obj_current <- obj_new
      
      # Track best solution
      r2_current <- 1 - obj_current / sum((y_v - mean(y_v))^2)
      if (r2_current > r2_best) {
        theta_best <- theta_current
        r2_best <- r2_current
      }
    }
    
    iteration_counts[i] <- iter
    
    # Set status if not converged
    if (!converged && convergence_status[i] == "") {
      convergence_status[i] <- "max_iterations"
    }
    
    # Update if improved
    if (r2_best > r2_voxel[v]) {
      theta_hat_voxel[v, ] <- theta_best
      r2_voxel[v] <- r2_best
      n_improved <- n_improved + 1
    }
    
    # Progress reporting
    if (verbose && i %% 10 == 0) {
      cat("  Processed", i, "/", n_hard, "voxels,",
          n_converged, "converged,", n_improved, "improved\n")
    }
  }
  
  if (verbose) {
    cat("  Gauss-Newton refinement complete:\n")
    cat("    Refined:", n_hard, "voxels\n")
    cat("    Converged:", n_converged, "(", 
        round(100 * n_converged / n_hard, 1), "%)\n")
    cat("    Improved:", n_improved, "(", 
        round(100 * n_improved / n_hard, 1), "%)\n")
    cat("    Mean iterations:", round(mean(iteration_counts), 1), "\n")
    
    # Convergence summary
    conv_table <- table(convergence_status)
    cat("    Convergence status:\n")
    for (status in names(conv_table)) {
      cat("      ", status, ":", conv_table[status], "\n")
    }
  }
  
  # Update queue labels for successfully refined voxels
  successfully_refined <- idx_hard[which(r2_voxel[idx_hard] > r2_orig[idx_hard])]
  if (length(successfully_refined) > 0) {
    queue_labels[successfully_refined] <- "easy"
  }
  
  # Return results
  list(
    theta_hat = theta_hat_voxel,
    r2 = r2_voxel,
    queue_labels = queue_labels,
    n_refined = n_hard,
    n_converged = n_converged,
    n_improved = n_improved,
    convergence_status = convergence_status,
    iteration_counts = iteration_counts,
    improvement_summary = list(
      mean_r2_improvement = mean(r2_voxel[idx_hard] - r2_orig[idx_hard], na.rm = TRUE),
      max_r2_improvement = max(r2_voxel[idx_hard] - r2_orig[idx_hard], na.rm = TRUE)
    )
  )
}

#' Calculate objective function for Gauss-Newton
#'
#' Returns `Inf` if the HRF predictor has near-zero magnitude, which
#' indicates that the system is singular and the amplitude cannot be
#' estimated reliably.
#'
#' @keywords internal
.calculate_objective_gn <- function(theta, y, S, t_hrf, hrf_interface, n_time) {
  # Validate theta
  if (!is.numeric(theta) || length(theta) != length(hrf_interface$parameter_names)) {
    stop(sprintf("Invalid theta: expected numeric vector of length %d, got %s of length %d",
                 length(hrf_interface$parameter_names), 
                 class(theta)[1], 
                 length(theta)))
  }
  
  # Generate HRF at current parameters
  hrf_vals <- hrf_interface$hrf_function(t_hrf, theta)

  # Convolve with stimulus
  conv_full <- stats::convolve(S[, 1], rev(hrf_vals), type = "open")
  x_pred_raw <- conv_full[seq_len(n_time)]

  denom <- sum(x_pred_raw^2)
  if (denom < 1e-8) {
    return(Inf)
  }

  # Fit amplitude analytically
  beta <- as.numeric(crossprod(x_pred_raw, y)) / denom
  x_pred <- beta * x_pred_raw

  # Return sum of squared residuals
  sum((y - x_pred)^2)
}

#' Get Jacobian matrix and residuals for Gauss-Newton
#'
#' Returns `NULL` if the HRF predictor has near-zero magnitude, which
#' indicates a singular system.
#'
#' @keywords internal
.get_jacobian_and_residuals <- function(theta, y, S, t_hrf, hrf_interface, n_time) {
  n_params <- length(theta)
  
  # Get Taylor basis (HRF and derivatives)
  taylor_basis <- hrf_interface$taylor_basis(theta, t_hrf)
  if (!is.matrix(taylor_basis)) {
    taylor_basis <- matrix(taylor_basis, ncol = n_params + 1)
  }
  
  # Convolve each basis function
  X_conv <- matrix(0, nrow = n_time, ncol = ncol(taylor_basis))
  for (j in seq_len(ncol(taylor_basis))) {
    conv_full <- stats::convolve(S[, 1], rev(taylor_basis[, j]), type = "open")
    X_conv[, j] <- conv_full[seq_len(n_time)]
  }
  
  # Fit amplitude for current HRF
  x_hrf <- X_conv[, 1]
  denom <- sum(x_hrf^2)
  if (denom < 1e-8) {
    return(NULL)
  }
  beta <- as.numeric(crossprod(x_hrf, y)) / denom
  
  # Residuals
  residuals <- y - beta * x_hrf
  
  # Jacobian w.r.t. parameters (chain rule through amplitude)
  # d(residual)/d(theta_k) = -beta * d(x_hrf)/d(theta_k) - x_hrf * d(beta)/d(theta_k)
  jacobian <- matrix(0, nrow = n_time, ncol = n_params)
  
  for (k in seq_len(n_params)) {
    dx_dtheta_k <- X_conv[, k + 1]
    
    # Derivative of beta w.r.t. theta_k
    # beta = (x_hrf' * y) / (x_hrf' * x_hrf)
    # d(beta)/d(theta_k) requires chain rule through denominator
    # d(denom)/d(theta_k) = d(x_hrf' * x_hrf)/d(theta_k) = 2 * x_hrf' * dx_dtheta_k
    dbeta_dtheta_k <- (as.numeric(crossprod(dx_dtheta_k, y)) -
      2 * beta * as.numeric(crossprod(x_hrf, dx_dtheta_k))) / denom
    
    # Full derivative
    jacobian[, k] <- -beta * dx_dtheta_k - dbeta_dtheta_k * x_hrf
  }
  
  list(
    jacobian = jacobian,
    residuals = residuals,
    amplitude = beta
  )
}