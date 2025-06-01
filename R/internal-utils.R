# Internal utility helpers
# Consolidated versions of .try_with_context(), .get_available_memory(),
# and .check_memory_available() stored in a dedicated environment to
# avoid accidental masking.

.fmriparametric_internal <- new.env(parent = emptyenv())

.fmriparametric_internal$try_with_context <- function(expr, context = "", fallback = NULL) {
  tryCatch(
    expr,
    error = function(e) {
      enhanced_msg <- sprintf(
        "Error in %s:\n  %s\n  Call stack: %s",
        context,
        e$message,
        paste(deparse(sys.calls()), collapse = " > ")
      )

      if (!is.null(fallback)) {
        warning(enhanced_msg, "\n  Attempting fallback...")
        tryCatch(
          fallback,
          error = function(e2) {
            stop(sprintf(
              "%s\n  Fallback also failed: %s",
              enhanced_msg, e2$message
            ), call. = FALSE)
          }
        )
      } else {
        stop(enhanced_msg, call. = FALSE)
      }
    }
  )
}

.fmriparametric_internal$get_available_memory <- function() {
  if (.Platform$OS.type == "windows") {
    memory.limit() * 1024^2
  } else if (Sys.info()["sysname"] == "Darwin") {
    # macOS
    tryCatch({
      vm_stat <- system("vm_stat", intern = TRUE)
      free_line <- grep("Pages free:", vm_stat, value = TRUE)
      inactive_line <- grep("Pages inactive:", vm_stat, value = TRUE)
      
      free_pages <- as.numeric(gsub(".*:\\s*(\\d+).*", "\\1", free_line))
      inactive_pages <- as.numeric(gsub(".*:\\s*(\\d+).*", "\\1", inactive_line))
      
      # Page size is typically 4096 bytes on macOS
      (free_pages + inactive_pages) * 4096
    }, error = function(e) {
      4e9  # Default 4GB
    })
  } else {
    # Linux
    tryCatch({
      as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE)) * 1024
    }, error = function(e) {
      4e9  # Default 4GB
    })
  }
}

.fmriparametric_internal$check_memory_available <- function(required_bytes, operation = "") {
  available <- .fmriparametric_internal$get_available_memory()
  
  # Handle case where available memory cannot be determined
  if (is.na(available) || is.null(available) || !is.numeric(available)) {
    # Can't determine memory, assume it's OK
    return(TRUE)
  }

  if (required_bytes > available * 0.8) {
    warning(sprintf(
      "Operation '%s' requires %.1f GB but only %.1f GB available",
      operation,
      required_bytes / 1e9,
      available / 1e9
    ))

    message("Consider:")
    message("  - Using smaller chunks")
    message("  - Increasing memory limit: memory.limit(size = ...)")
    message("  - Closing other applications")

    return(FALSE)
  }

  TRUE
}

# Backwards compatible wrappers
.try_with_context <- function(...) .fmriparametric_internal$try_with_context(...)
.get_available_memory <- function(...) .fmriparametric_internal$get_available_memory(...)
.check_memory_available <- function(...) .fmriparametric_internal$check_memory_available(...)
