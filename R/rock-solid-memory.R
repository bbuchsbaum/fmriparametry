#' Rock Solid Memory Management Functions
#'
#' Memory protection functions that prevent out-of-memory errors and
#' optimize memory usage for large-scale analyses.

#' Check available memory and estimate requirements
#' @keywords internal
.check_memory_requirements <- function(n_voxels, n_timepoints, n_params = 3,
                                       safety_factor = 2.0, caller = "engine") {
  # Get system memory info
  sys_info <- Sys.info()
  
  # Platform-specific memory detection
  total_ram <- tryCatch({
    if (sys_info["sysname"] == "Windows") {
      # Windows: use memory.size
      memory.size(max = TRUE) * 1024 * 1024  # Convert MB to bytes
    } else if (sys_info["sysname"] == "Darwin") {
      # macOS: use system profiler
      cmd <- "sysctl -n hw.memsize"
      as.numeric(system(cmd, intern = TRUE))
    } else {
      # Linux: read from /proc/meminfo
      meminfo <- readLines("/proc/meminfo", n = 1)
      as.numeric(gsub("[^0-9]", "", meminfo)) * 1024  # KB to bytes
    }
  }, error = function(e) {
    warning(caller, ": Could not detect system memory. Assuming 8GB.")
    8 * 1024^3  # 8 GB default
  })
  
  # Estimate memory requirements
  # Main data matrices
  bytes_per_double <- 8
  
  memory_breakdown <- list(
    Y_matrix = n_voxels * n_timepoints * bytes_per_double,
    design_matrix = n_timepoints * (n_params + 1) * bytes_per_double,
    parameters = n_voxels * n_params * bytes_per_double,
    coefficients = n_voxels * (n_params + 1) * bytes_per_double,
    residuals = n_voxels * n_timepoints * bytes_per_double,
    working_memory = n_timepoints * n_timepoints * bytes_per_double  # QR decomp
  )
  
  total_required <- sum(unlist(memory_breakdown)) * safety_factor
  
  # Get current memory usage
  current_usage <- tryCatch({
    gc()  # Force garbage collection
    sum(gc()[, "used"]) * 1024 * 1024  # Convert MB to bytes
  }, error = function(e) 0)
  
  available_memory <- total_ram - current_usage
  
  # Memory status
  memory_ok <- total_required < available_memory * 0.8  # Use max 80% of available
  
  if (!memory_ok) {
    # Calculate recommended chunk size
    bytes_per_voxel <- total_required / n_voxels
    safe_n_voxels <- floor(available_memory * 0.8 / bytes_per_voxel)
    recommended_chunks <- ceiling(n_voxels / safe_n_voxels)
  } else {
    recommended_chunks <- 1
  }
  
  list(
    total_ram_gb = total_ram / 1024^3,
    available_gb = available_memory / 1024^3,
    required_gb = total_required / 1024^3,
    memory_ok = memory_ok,
    recommended_chunks = recommended_chunks,
    breakdown = lapply(memory_breakdown, function(x) x / 1024^3)  # Convert to GB
  )
}

#' Process data in memory-efficient chunks
#' @keywords internal
.chunked_processing <- function(data, chunk_size, process_fun, 
                                combine_fun = rbind, 
                                progress = TRUE, ...) {
  n_total <- nrow(data)
  n_chunks <- ceiling(n_total / chunk_size)
  
  if (progress) {
    cat("Processing", n_total, "items in", n_chunks, "chunks...\n")
  }
  
  results <- vector("list", n_chunks)
  
  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_total)
    
    if (progress && i %% 10 == 0) {
      cat("  Chunk", i, "/", n_chunks, 
          "(", round(100 * i / n_chunks), "% complete)\n")
    }
    
    # Process chunk
    chunk_data <- data[start_idx:end_idx, , drop = FALSE]
    
    results[[i]] <- tryCatch({
      process_fun(chunk_data, ...)
    }, error = function(e) {
      warning("Chunk ", i, " failed: ", e$message)
      NULL
    })
    
    # Aggressive garbage collection
    if (i %% 5 == 0) {
      gc()
    }
  }
  
  # Combine results
  combined <- tryCatch({
    do.call(combine_fun, results[!sapply(results, is.null)])
  }, error = function(e) {
    warning("Result combination failed: ", e$message)
    results
  })
  
  if (progress) {
    cat("Chunked processing complete.\n")
  }
  
  combined
}

#' Memory-aware matrix operations
#' @keywords internal
.memory_safe_matmul <- function(A, B, max_elements = 1e8) {
  dims_A <- dim(A)
  dims_B <- dim(B)
  
  # Check if standard multiplication is safe
  output_elements <- dims_A[1] * dims_B[2]
  
  if (output_elements < max_elements) {
    # Safe to use standard multiplication
    return(A %*% B)
  }
  
  # Use blocked multiplication
  block_size <- floor(sqrt(max_elements))
  n_blocks_i <- ceiling(dims_A[1] / block_size)
  n_blocks_j <- ceiling(dims_B[2] / block_size)
  
  # Pre-allocate result
  result <- matrix(0, nrow = dims_A[1], ncol = dims_B[2])
  
  for (i in seq_len(n_blocks_i)) {
    i_start <- (i - 1) * block_size + 1
    i_end <- min(i * block_size, dims_A[1])
    
    for (j in seq_len(n_blocks_j)) {
      j_start <- (j - 1) * block_size + 1
      j_end <- min(j * block_size, dims_B[2])
      
      # Compute block
      result[i_start:i_end, j_start:j_end] <- 
        A[i_start:i_end, ] %*% B[, j_start:j_end]
    }
    
    # Periodic garbage collection
    if (i %% 5 == 0) gc()
  }
  
  result
}

#' Safe pre-allocation with memory checks
#' @keywords internal
.safe_allocate <- function(dims, init_value = 0, type = "numeric", 
                           caller = "allocation") {
  if (length(dims) == 1) dims <- c(dims, 1)
  
  n_elements <- prod(dims)
  bytes_per_element <- switch(type,
                              numeric = 8,
                              integer = 4,
                              logical = 4,
                              character = 64,  # Rough estimate
                              8)  # Default
  
  required_bytes <- n_elements * bytes_per_element
  
  # Check if allocation is reasonable
  if (required_bytes > 4 * 1024^3) {  # > 4GB
    warning(caller, ": Large allocation requested (", 
            round(required_bytes / 1024^3, 1), " GB). ",
            "Consider chunked processing.")
  }
  
  # Try allocation with fallback
  allocated <- tryCatch({
    if (type == "numeric") {
      matrix(init_value, nrow = dims[1], ncol = dims[2])
    } else if (type == "integer") {
      matrix(as.integer(init_value), nrow = dims[1], ncol = dims[2])
    } else if (type == "logical") {
      matrix(as.logical(init_value), nrow = dims[1], ncol = dims[2])
    } else {
      matrix(init_value, nrow = dims[1], ncol = dims[2])
    }
  }, error = function(e) {
    stop(caller, ": Memory allocation failed for ", 
         dims[1], " x ", dims[2], " ", type, " matrix. ",
         "Error: ", e$message, call. = FALSE)
  })
  
  allocated
}

#' Clean up temporary objects and force garbage collection
#' @keywords internal
.cleanup_workspace <- function(objects_to_remove = NULL, 
                               gc_runs = 2, 
                               verbose = FALSE) {
  if (!is.null(objects_to_remove)) {
    # Remove specified objects from parent environment
    parent_env <- parent.frame()
    existing <- objects_to_remove[objects_to_remove %in% ls(parent_env)]
    if (length(existing) > 0) {
      rm(list = existing, envir = parent_env)
      if (verbose) {
        cat("Removed", length(existing), "temporary objects.\n")
      }
    }
  }
  
  # Multiple garbage collection runs
  mem_before <- gc()[, "used"]
  
  for (i in seq_len(gc_runs)) {
    gc()
  }
  
  mem_after <- gc()[, "used"]
  mem_freed <- sum(mem_before - mem_after)
  
  if (verbose && mem_freed > 0) {
    cat("Freed", round(mem_freed, 1), "MB of memory.\n")
  }
  
  invisible(mem_freed)
}

#' Monitor memory usage during operations
#' @keywords internal
.with_memory_tracking <- function(expr, label = "Operation", 
                                  warn_threshold_gb = 1, 
                                  caller = "memory_track") {
  # Initial memory state
  gc()
  mem_start <- sum(gc()[, "used"])
  time_start <- Sys.time()
  
  # Execute expression
  result <- tryCatch({
    expr
  }, error = function(e) {
    # Check if memory-related error
    if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
      stop(caller, ": Out of memory during ", label, ". ",
           "Consider using smaller chunks or increasing memory limits.",
           call. = FALSE)
    } else {
      stop(e)
    }
  })
  
  # Final memory state
  gc()
  mem_end <- sum(gc()[, "used"])
  time_end <- Sys.time()
  
  # Memory statistics
  mem_used_mb <- mem_end - mem_start
  mem_peak_mb <- max(gc()[, "max used"]) - mem_start
  time_elapsed <- as.numeric(time_end - time_start, units = "secs")
  
  if (mem_peak_mb > warn_threshold_gb * 1024) {
    warning(caller, ": ", label, " used ", round(mem_peak_mb / 1024, 1), 
            " GB of memory (peak). Consider optimization if this is excessive.")
  }
  
  attr(result, "memory_stats") <- list(
    used_mb = mem_used_mb,
    peak_mb = mem_peak_mb,
    time_seconds = time_elapsed,
    label = label
  )
  
  result
}