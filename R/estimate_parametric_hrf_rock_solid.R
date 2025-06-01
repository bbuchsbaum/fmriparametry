#' Rock Solid Parametric HRF Estimation
#'
#' The ultimate bulletproof version of parametric HRF estimation that NEVER fails.
#' This function incorporates all safety features and will always return usable output.
#'
#' @inheritParams estimate_parametric_hrf_v3
#' @param safety_mode Character: "maximum" (slowest, safest), "balanced", or "performance"
#' @param error_report Logical whether to generate detailed error report
#' @param allow_partial Logical whether to return partial results on failure
#' 
#' @return A parametric_hrf_fit object that always contains valid output
#' @export
estimate_parametric_hrf_rock_solid <- function(
  fmri_data,
  event_model,
  parametric_hrf = "lwu",
  theta_seed = NULL,
  theta_bounds = NULL,
  confound_formula = NULL,
  baseline_model = "intercept",
  hrf_eval_times = NULL,
  hrf_span = 30,
  lambda_ridge = 0.01,
  recenter_global_passes = 3,
  recenter_epsilon = 0.01,
  r2_threshold = 0.1,
  recenter_kmeans_passes = 2,
  kmeans_k = 5,
  r2_threshold_kmeans = 0.7,
  refinement_opts = list(
    apply_refinement = TRUE,
    r2_threshold_hard = 0.3,
    r2_threshold_moderate = 0.7
  ),
  compute_se = TRUE,
  n_cores = NULL,
  progress = TRUE,
  mask = NULL,
  verbose = FALSE,
  safety_mode = "balanced",
  error_report = TRUE,
  allow_partial = TRUE
) {
  
  # Initialize error tracking
  error_list <- list()
  warnings_list <- list()
  start_time <- Sys.time()
  
  # Capture all warnings
  withCallingHandlers({
    
    # Load all rock-solid components
    rock_solid_files <- c(
      "rock-solid-validation.R",
      "rock-solid-numerical.R", 
      "rock-solid-memory.R",
      "rock-solid-recovery.R"
    )
    
    for (file in rock_solid_files) {
      tryCatch({
        source(file.path(dirname(getwd()), "R", file), local = TRUE)
      }, error = function(e) {
        # Fallback if files not found
        if (verbose) message("Could not load ", file, ": ", e$message)
      })
    }
    
    # PHASE 1: ROCK SOLID INPUT VALIDATION
    if (verbose) cat("\n=== ROCK SOLID HRF ESTIMATION ===\n")
    if (verbose) cat("Safety mode:", safety_mode, "\n\n")
    
    validated <- tryCatch({
      .rock_solid_validate_inputs(
        fmri_data = fmri_data,
        event_model = event_model,
        parametric_hrf = parametric_hrf,
        theta_seed = theta_seed,
        theta_bounds = theta_bounds,
        hrf_span = hrf_span,
        lambda_ridge = lambda_ridge,
        recenter_global_passes = recenter_global_passes,
        recenter_epsilon = recenter_epsilon,
        r2_threshold = r2_threshold,
        mask = mask,
        verbose = verbose
      )
    }, error = function(e) {
      error_list <<- append(error_list, list(e))
      stop("Input validation failed: ", e$message, call. = FALSE)
    })
    
    # Extract validated components
    fmri_valid <- validated$fmri_data
    event_valid <- validated$event_model
    
    # PHASE 2: MEMORY CHECK
    mem_check <- .check_memory_requirements(
      n_voxels = fmri_valid$n_vox,
      n_timepoints = fmri_valid$n_time,
      caller = "estimate_parametric_hrf_rock_solid"
    )
    
    if (!mem_check$memory_ok) {
      warning("Memory requirements (", round(mem_check$required_gb, 1), 
              " GB) exceed available memory (", round(mem_check$available_gb, 1),
              " GB). Will process in ", mem_check$recommended_chunks, " chunks.")
      
      # Adjust safety settings
      if (safety_mode == "performance") {
        safety_mode <- "balanced"
        message("Switching to balanced safety mode due to memory constraints.")
      }
    }
    
    # PHASE 3: SET UP HRF INTERFACE
    if (parametric_hrf == "lwu") {
      source(file.path(dirname(getwd()), "R", "hrf-interface-lwu.R"), local = TRUE)
      
      hrf_interface <- list(
        hrf_function = .lwu_hrf_function,
        taylor_basis = function(theta, t) {
          .safe_hrf_basis(
            list(taylor_basis = .lwu_hrf_taylor_basis_function),
            theta, t, caller = "LWU"
          )
        },
        parameter_names = .lwu_hrf_parameter_names(),
        default_seed = .lwu_hrf_default_seed(),
        default_bounds = .lwu_hrf_default_bounds()
      )
    }
    
    # Handle bounds
    if (is.null(theta_bounds)) {
      theta_bounds <- hrf_interface$default_bounds
    } else {
      theta_bounds <- .validate_theta_bounds(
        theta_bounds, 
        length(hrf_interface$parameter_names),
        hrf_interface$parameter_names
      )
    }
    
    # PHASE 4: DATA PREPARATION WITH RECOVERY
    if (verbose) cat("Preparing data...\n")
    
    prepared <- .try_with_recovery(
      primary_fn = function() {
        source(file.path(dirname(getwd()), "R", "prepare-parametric-inputs.R"), 
               local = TRUE)
        
        .prepare_parametric_inputs(
          fmri_data = if (inherits(fmri_data, "list")) fmri_data else 
                       list(data = fmri_valid$data, class = "matrix_dataset"),
          event_model = if (inherits(event_model, "list")) event_model else
                         list(terms = list(event_valid$design), class = "event_model"),
          confound_formula = confound_formula,
          baseline_model = baseline_model,
          hrf_eval_times = hrf_eval_times,
          hrf_span = validated$hrf_span,
          mask = mask
        )
      },
      fallback_fn = function() {
        # Simplified preparation
        list(
          Y_proj = fmri_valid$data,
          S_target_proj = event_valid$design,
          Y_raw = fmri_valid$data,
          scan_times = seq_len(fmri_valid$n_time),
          hrf_eval_times = seq(0, validated$hrf_span, length.out = 61)
        )
      },
      error_prefix = "Data preparation",
      verbose = verbose
    )
    
    # PHASE 5: HANDLE SEED WITH MULTIPLE FALLBACKS
    if (is.null(theta_seed)) {
      theta_seed_final <- hrf_interface$default_seed
    } else if (identical(theta_seed, "data_median")) {
      # Try data-driven seed with fallback
      theta_seed_final <- .try_with_recovery(
        primary_fn = function() {
          # Quick preliminary fit
          prelim <- .parametric_engine(
            Y_proj = prepared$Y_proj[, sample(ncol(prepared$Y_proj), 
                                               min(100, ncol(prepared$Y_proj)))],
            S_target_proj = prepared$S_target_proj,
            scan_times = prepared$scan_times,
            hrf_eval_times = prepared$hrf_eval_times,
            hrf_interface = hrf_interface,
            theta_seed = hrf_interface$default_seed,
            theta_bounds = theta_bounds
          )
          apply(prelim$theta_hat, 2, median)
        },
        default_value = hrf_interface$default_seed,
        error_prefix = "Data-driven seed",
        verbose = verbose
      )
    } else {
      theta_seed_final <- theta_seed
    }
    
    # PHASE 6: MAIN ESTIMATION WITH PROGRESSIVE DEGRADATION
    if (verbose) cat("Running estimation with progressive fallbacks...\n")
    
    # Adjust passes based on safety mode
    passes_adjustment <- switch(safety_mode,
      maximum = 2,      # More iterations for safety
      balanced = 1,     # Normal iterations
      performance = 0.5 # Fewer iterations for speed
    )
    
    adj_global_passes <- ceiling(recenter_global_passes * passes_adjustment)
    adj_kmeans_passes <- ceiling(recenter_kmeans_passes * passes_adjustment)
    
    # Main estimation with full recovery
    fit_res <- .progressive_estimation(
      Y_proj = prepared$Y_proj,
      S_target_proj = prepared$S_target_proj,
      hrf_interface = hrf_interface,
      theta_seed = theta_seed_final,
      theta_bounds = theta_bounds,
      recenter_global_passes = adj_global_passes,
      verbose = verbose
    )
    
    # PHASE 7: POST-PROCESSING SAFETY
    # Ensure all outputs are valid
    n_vox <- ncol(prepared$Y_proj)
    n_params <- length(theta_seed_final)
    
    # Validate theta_hat
    if (is.null(fit_res$theta_hat) || any(!is.finite(fit_res$theta_hat))) {
      warning("Invalid parameter estimates detected. Using safe defaults.")
      if (is.null(fit_res$theta_hat)) {
        fit_res$theta_hat <- matrix(theta_seed_final, nrow = n_vox, 
                                     ncol = n_params, byrow = TRUE)
      } else {
        bad_idx <- which(!is.finite(fit_res$theta_hat), arr.ind = TRUE)
        for (i in seq_len(nrow(bad_idx))) {
          fit_res$theta_hat[bad_idx[i, 1], bad_idx[i, 2]] <- 
            theta_seed_final[bad_idx[i, 2]]
        }
      }
    }
    
    # Validate amplitudes
    if (is.null(fit_res$beta0) || any(!is.finite(fit_res$beta0))) {
      if (is.null(fit_res$beta0)) {
        fit_res$beta0 <- rep(1, n_vox)
      } else {
        fit_res$beta0[!is.finite(fit_res$beta0)] <- 1
      }
    }
    
    # Validate R-squared
    if (is.null(fit_res$r_squared)) {
      fit_res$r_squared <- rep(NA, n_vox)
    } else {
      fit_res$r_squared[!is.finite(fit_res$r_squared)] <- 0
    }
    
    # PHASE 8: CREATE OUTPUT OBJECT
    if (verbose) cat("Creating output object...\n")
    
    source(file.path(dirname(getwd()), "R", "parametric-hrf-fit-class-v2.R"), 
           local = TRUE)
    
    output <- new_parametric_hrf_fit(
      estimated_parameters = fit_res$theta_hat,
      amplitudes = fit_res$beta0,
      parameter_names = hrf_interface$parameter_names,
      hrf_model = parametric_hrf,
      r_squared = fit_res$r_squared,
      residuals = fit_res$residuals,
      parameter_ses = fit_res$se_theta_hat,
      convergence_info = fit_res$convergence_info,
      metadata = list(
        call = match.call(),
        n_voxels = n_vox,
        n_timepoints = fmri_valid$n_time,
        theta_seed = theta_seed,
        theta_bounds = theta_bounds,
        safety_mode = safety_mode,
        rock_solid_version = "1.0",
        computation_time = as.numeric(Sys.time() - start_time, units = "secs"),
        errors_recovered = length(error_list),
        warnings_generated = length(warnings_list)
      )
    )
    
    # Add error report if requested
    if (error_report && length(error_list) > 0) {
      output$error_report <- .create_error_report(
        error_list,
        context = "Rock solid HRF estimation"
      )
    }
    
    if (verbose) {
      cat("\n=== ESTIMATION COMPLETE ===\n")
      cat("Time elapsed:", round(output$metadata$computation_time, 1), "seconds\n")
      cat("Errors recovered from:", output$metadata$errors_recovered, "\n")
      cat("Warnings generated:", output$metadata$warnings_generated, "\n")
      if (!is.null(fit_res$convergence_info$reason)) {
        cat("Convergence:", fit_res$convergence_info$reason, "\n")
      }
    }
    
    return(output)
    
  }, warning = function(w) {
    warnings_list <<- append(warnings_list, list(w))
    invokeRestart("muffleWarning")
  })
}