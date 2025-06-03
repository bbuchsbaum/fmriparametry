#' Set up parallel backend
#'
#' Configures the parallel processing backend for voxel-wise operations.
#' Supports both future and parallel packages with automatic fallback.
#'
#' @param n_cores Number of cores to use. NULL for auto-detection.
#' @param verbose Logical whether to print setup messages
#' 
#' @return List with parallel configuration info
#' @keywords internal
.setup_parallel_backend <- function(n_cores = NULL, verbose = TRUE) {
  # Check available cores
  available_cores <- parallel::detectCores()
  
  # Determine cores to use
  if (is.null(n_cores)) {
    # Use all but one core
    n_cores <- max(1, available_cores - 1)
  } else {
    n_cores <- min(n_cores, available_cores)
  }
  
  # Single core - no parallel setup needed
  if (n_cores == 1) {
    if (verbose) message("Using sequential processing (1 core)")
    return(list(
      backend = "sequential",
      n_cores = 1,
      cleanup = function() invisible(NULL)
    ))
  }
  
  # Try to use future package first (more flexible)
  if (requireNamespace("future", quietly = TRUE) && 
      requireNamespace("future.apply", quietly = TRUE)) {
    
    if (verbose) message("Setting up future backend with ", n_cores, " cores")
    
    # Set up multicore/multisession plan
    if (.Platform$OS.type == "unix") {
      future::plan(future::multicore, workers = n_cores)
      backend_type <- "future_multicore"
    } else {
      future::plan(future::multisession, workers = n_cores)
      backend_type <- "future_multisession"
    }
    
    return(list(
      backend = backend_type,
      n_cores = n_cores,
      cleanup = function() future::plan(future::sequential)
    ))
  }
  
  # Fall back to parallel package
  if (verbose) message("Setting up parallel backend with ", n_cores, " cores")
  
  if (.Platform$OS.type == "unix") {
    # Use mclapply on Unix-like systems
    return(list(
      backend = "mclapply",
      n_cores = n_cores,
      cleanup = function() invisible(NULL)
    ))
  } else {
    # Use parLapply on Windows
    cl <- parallel::makeCluster(n_cores)
    return(list(
      backend = "parLapply",
      n_cores = n_cores,
      cluster = cl,
      cleanup = function() parallel::stopCluster(cl)
    ))
  }
}

#' Parallel voxel processing
#'
#' Generic function for parallel processing of voxel-wise operations.
#' Handles chunking, progress reporting, and different backend types.
#'
#' @param voxel_indices Indices of voxels to process
#' @param process_function Function to apply to each voxel/chunk
#' @param parallel_config Configuration from .setup_parallel_backend
#' @param chunk_size Size of chunks ("auto" or numeric)
#' @param progress Logical whether to show progress
#' @param ... Additional arguments passed to process_function
#'
#' @return Combined results from all voxels/chunks
#' @keywords internal
.parallel_voxel_processing <- function(
  voxel_indices,
  process_function,
  parallel_config,
  chunk_size = "auto",
  progress = TRUE,
  ...
) {
  n_vox <- length(voxel_indices)
  
  # Determine chunk size
  if (identical(chunk_size, "auto")) {
    # Balance between overhead and memory usage
    chunk_size <- ceiling(n_vox / (parallel_config$n_cores * 10))
    chunk_size <- max(10, min(chunk_size, 1000))
  }
  
  # Create chunks
  n_chunks <- ceiling(n_vox / chunk_size)
  chunks <- split(voxel_indices, ceiling(seq_along(voxel_indices) / chunk_size))
  
  if (progress && n_chunks > 1) {
    message("Processing ", n_vox, " voxels in ", n_chunks, " chunks...")
  }
  
  # Sequential processing
  if (parallel_config$backend == "sequential") {
    results <- lapply(chunks, function(chunk) {
      process_function(chunk, ...)
    })
    return(do.call(c, results))
  }
  
  # Future-based processing
  if (grepl("^future", parallel_config$backend)) {
    results <- future.apply::future_lapply(
      chunks,
      function(chunk) process_function(chunk, ...),
      future.seed = TRUE
    )
    return(do.call(c, results))
  }
  
  # mclapply (Unix)
  if (parallel_config$backend == "mclapply") {
    results <- parallel::mclapply(
      chunks,
      function(chunk) process_function(chunk, ...),
      mc.cores = parallel_config$n_cores
    )
    return(do.call(c, results))
  }
  
  # parLapply (Windows)
  if (parallel_config$backend == "parLapply") {
    # Export necessary variables to cluster
    parallel::clusterExport(
      parallel_config$cluster,
      varlist = "process_function",
      envir = environment()
    )
    
    results <- parallel::parLapply(
      parallel_config$cluster,
      chunks,
      function(chunk) process_function(chunk, ...)
    )
    return(do.call(c, results))
  }
  
  stop("Unknown parallel backend: ", parallel_config$backend)
}

#' Parallel K-means cluster processing
#'
#' Process K-means clusters in parallel, with each cluster getting its own
#' worker for load balancing.
#'
#' @param cluster_data List with cluster assignments and data
#' @param process_cluster_fn Function to process each cluster
#' @param parallel_config Parallel configuration
#' @param verbose Logical for progress messages
#' @param ... Additional arguments for process_cluster_fn
#'
#' @return Combined results from all clusters
#' @keywords internal
.parallel_kmeans_processing <- function(
  cluster_data,
  process_cluster_fn,
  parallel_config,
  verbose = FALSE,
  ...
) {
  n_clusters <- length(unique(cluster_data$cluster_assignments))
  
  if (verbose) {
    message("Processing ", n_clusters, " K-means clusters in parallel...")
  }
  
  # Create cluster tasks
  cluster_tasks <- lapply(seq_len(n_clusters), function(k) {
    list(
      cluster_id = k,
      voxel_indices = which(cluster_data$cluster_assignments == k),
      cluster_center = cluster_data$cluster_centers[k, ]
    )
  })
  
  # Sequential processing
  if (parallel_config$backend == "sequential") {
    results <- lapply(cluster_tasks, function(task) {
      process_cluster_fn(task, ...)
    })
    return(results)
  }
  
  # Future-based processing
  if (grepl("^future", parallel_config$backend)) {
    results <- future.apply::future_lapply(
      cluster_tasks,
      function(task) process_cluster_fn(task, ...),
      future.seed = TRUE
    )
    return(results)
  }
  
  # mclapply (Unix)
  if (parallel_config$backend == "mclapply") {
    results <- parallel::mclapply(
      cluster_tasks,
      function(task) process_cluster_fn(task, ...),
      mc.cores = parallel_config$n_cores
    )
    return(results)
  }
  
  # parLapply (Windows)
  if (parallel_config$backend == "parLapply") {
    parallel::clusterExport(
      parallel_config$cluster,
      varlist = c("process_cluster_fn"),
      envir = environment()
    )
    
    results <- parallel::parLapply(
      parallel_config$cluster,
      cluster_tasks,
      function(task) process_cluster_fn(task, ...)
    )
    return(results)
  }
}

#' Parallel local re-centering
#'
#' Perform local re-centering for moderate voxels in parallel.
#'
#' @param moderate_indices Indices of moderate voxels
#' @param fit_data Current fit data
#' @param prepared_data Prepared input data
#' @param hrf_interface HRF interface functions
#' @param parallel_config Parallel configuration
#' @param ... Additional parameters
#'
#' @return Updated fit results for moderate voxels
#' @keywords internal
.parallel_local_recentering <- function(
  moderate_indices,
  fit_data,
  prepared_data,
  hrf_interface,
  parallel_config,
  theta_bounds,
  lambda_ridge = 0.01,
  verbose = FALSE
) {
  if (length(moderate_indices) == 0) {
    return(list())
  }
  
  if (verbose) {
    message("Parallel local re-centering for ", length(moderate_indices), " voxels...")
  }
  
  # Define processing function for each voxel
  process_voxel <- function(v) {
    # Use voxel's current estimate as expansion point
    theta_local <- fit_data$theta_hat[v, ]
    
    # Single Taylor pass for this voxel
    local_res <- .parametric_engine(
      Y_proj = prepared_data$Y_proj[, v, drop = FALSE],
      S_target_proj = prepared_data$S_target_proj,
      hrf_eval_times = prepared_data$hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_seed = theta_local,
      theta_bounds = theta_bounds,
      lambda_ridge = lambda_ridge,
      verbose = FALSE
    )
    
    # Return update info
    list(
      voxel_idx = v,
      improved = local_res$r_squared[1] > fit_data$r_squared[v],
      theta_new = local_res$theta_hat[1, ],
      beta_new = local_res$beta0[1],
      r2_new = local_res$r_squared[1],
      r2_old = fit_data$r_squared[v]
    )
  }
  
  # Process in parallel
  results <- .parallel_voxel_processing(
    voxel_indices = moderate_indices,
    process_function = process_voxel,
    parallel_config = parallel_config,
    chunk_size = 1,  # One voxel per task for load balancing
    progress = verbose
  )
  
  return(results)
}

#' Parallel Gauss-Newton optimization
#'
#' Perform Gauss-Newton optimization for hard voxels in parallel.
#'
#' @param hard_indices Indices of hard voxels
#' @param fit_data Current fit data
#' @param prepared_data Prepared input data
#' @param hrf_interface HRF interface functions
#' @param parallel_config Parallel configuration
#' @param ... Additional GN parameters
#'
#' @return Updated fit results for hard voxels
#' @keywords internal
.parallel_gauss_newton <- function(
  hard_indices,
  fit_data,
  prepared_data,
  hrf_interface,
  parallel_config,
  theta_bounds,
  max_iter_gn = 5,
  tol_gn = 1e-4,
  lambda_ridge = 0.01,
  verbose = FALSE
) {
  if (length(hard_indices) == 0) {
    return(list())
  }
  
  if (verbose) {
    message("Parallel Gauss-Newton for ", length(hard_indices), " voxels...")
  }
  
  # Define processing function for each voxel
  process_hard_voxel <- function(v) {
    y_v <- prepared_data$Y_proj[, v]
    theta_current <- fit_data$theta_hat[v, ]
    theta_best <- theta_current
    r2_best <- fit_data$r_squared[v]
    
    converged <- FALSE
    iter <- 0
    
    # Gauss-Newton iterations
    for (iter in seq_len(max_iter_gn)) {
      # Get Jacobian and residuals
      jacob_info <- .get_jacobian_and_residuals(
        theta_current, y_v, prepared_data$S_target_proj, 
        prepared_data$hrf_eval_times, hrf_interface, nrow(prepared_data$Y_proj)
      )
      
      if (is.null(jacob_info)) break
      
      # Compute update
      JtJ <- crossprod(jacob_info$jacobian)
      Jtr <- crossprod(jacob_info$jacobian, jacob_info$residuals)
      JtJ_ridge <- JtJ + lambda_ridge * diag(ncol(JtJ))
      
      delta <- tryCatch({
        -solve(JtJ_ridge, Jtr)
      }, error = function(e) NULL)
      
      if (is.null(delta)) break
      
      # Line search
      alpha <- 1.0
      for (ls in 1:10) {
        theta_new <- theta_current + alpha * as.numeric(delta)
        theta_new <- pmax(theta_bounds$lower, pmin(theta_new, theta_bounds$upper))
        
        obj_new <- .calculate_objective_gn(
          theta_new, y_v, prepared_data$S_target_proj,
          prepared_data$hrf_eval_times, hrf_interface, nrow(prepared_data$Y_proj)
        )
        
        obj_current <- sum(jacob_info$residuals^2)
        if (obj_new < obj_current) break
        alpha <- alpha * 0.5
      }
      
      # Check convergence
      if (sqrt(sum((theta_new - theta_current)^2)) < tol_gn) {
        converged <- TRUE
        break
      }
      
      theta_current <- theta_new
    }
    
    # Calculate final R²
    hrf_vals <- hrf_interface$hrf_function(prepared_data$hrf_eval_times, theta_current)
    conv_full <- stats::convolve(prepared_data$S_target_proj[, 1], rev(hrf_vals), type = "open")
    x_pred <- conv_full[seq_len(nrow(prepared_data$Y_proj))]
    beta_new <- sum(x_pred * y_v) / sum(x_pred^2)
    r2_new <- 1 - sum((y_v - beta_new * x_pred)^2) / sum((y_v - mean(y_v))^2)
    
    list(
      voxel_idx = v,
      converged = converged,
      iterations = iter,
      improved = r2_new > r2_best,
      theta_new = theta_current,
      beta_new = beta_new,
      r2_new = r2_new,
      r2_old = r2_best
    )
  }
  
  # Process in parallel
  results <- .parallel_voxel_processing(
    voxel_indices = hard_indices,
    process_function = process_hard_voxel,
    parallel_config = parallel_config,
    chunk_size = 1,  # One voxel per task
    progress = verbose
  )
  
  return(results)
}

#' Parallel standard error calculation
#'
#' Calculate standard errors for multiple voxels in parallel using the
#' Delta method.
#'
#' @param voxel_indices Indices of voxels to process
#' @param fit_data Current fit data
#' @param prepared_data Prepared input data
#' @param hrf_interface HRF interface functions
#' @param parallel_config Parallel configuration
#' @param lambda_ridge Ridge penalty
#' @param verbose Logical for progress
#'
#' @return Matrix of standard errors
#' @keywords internal
.parallel_se_calculation <- function(
  voxel_indices,
  fit_data,
  prepared_data,
  hrf_interface,
  parallel_config,
  lambda_ridge = 0.01,
  verbose = FALSE
) {
  if (length(voxel_indices) == 0) {
    return(matrix(NA, 0, ncol(fit_data$theta_hat)))
  }
  
  if (verbose) {
    message("Parallel SE calculation for ", length(voxel_indices), " voxels...")
  }
  
  n_params <- ncol(fit_data$theta_hat)
  
  # Define processing function
  calculate_se_voxel <- function(v) {
    theta_v <- fit_data$theta_hat[v, ]
    
    # Get Taylor basis
    taylor_basis <- hrf_interface$taylor_basis(theta_v, prepared_data$hrf_eval_times)
    if (!is.matrix(taylor_basis)) {
      taylor_basis <- matrix(taylor_basis, ncol = n_params + 1)
    }
    
    # Convolve derivatives
    X_deriv <- matrix(0, nrow = nrow(prepared_data$Y_proj), ncol = n_params)
    for (j in seq_len(n_params)) {
      conv_full <- stats::convolve(
        prepared_data$S_target_proj[, 1], 
        rev(taylor_basis[, j + 1]), 
        type = "open"
      )
      X_deriv[, j] <- conv_full[seq_len(nrow(prepared_data$Y_proj))]
    }
    
    # Estimate error variance
    var_beta <- if (!is.null(fit_data$residuals)) {
      sum(fit_data$residuals[, v]^2) / (nrow(prepared_data$Y_proj) - n_params - 1)
    } else {
      1
    }
    
    # Calculate covariance matrix
    cov_theta <- tryCatch({
      var_beta * solve(crossprod(X_deriv) + lambda_ridge * diag(n_params))
    }, error = function(e) {
      matrix(NA, n_params, n_params)
    })
    
    sqrt(diag(cov_theta))
  }
  
  # Process in parallel
  se_list <- .parallel_voxel_processing(
    voxel_indices = voxel_indices,
    process_function = calculate_se_voxel,
    parallel_config = parallel_config,
    chunk_size = "auto",
    progress = verbose
  )
  
  # Convert to matrix
  do.call(rbind, se_list)
}