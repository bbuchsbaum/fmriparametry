#' Memory mapped parametric engine
#'
#' Processes large fMRI datasets stored on disk in chunks. The data file
#' should be an RDS file containing a numeric matrix (time x voxels).
#'
#' @param fmri_file Path to RDS file containing the fMRI matrix
#' @param event_model Design matrix for events
#' @param hrf_interface HRF interface object
#' @param chunk_size Number of voxels to process per chunk
#' @param ... Additional arguments passed to `.parametric_engine`
#' @return List with `theta_hat`, `beta0` and `r_squared`
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x
.memory_mapped_engine <- function(fmri_file, event_model, hrf_interface,
                                  chunk_size = 1000, ...) {
  Y_all <- readRDS(fmri_file)
  n_time <- nrow(Y_all)
  n_vox <- ncol(Y_all)

  theta_hat <- matrix(NA_real_, n_vox, length(hrf_interface$parameter_names))
  beta0 <- numeric(n_vox)
  r2 <- numeric(n_vox)

  chunk_starts <- seq(1, n_vox, by = chunk_size)
  for (start in chunk_starts) {
    end <- min(start + chunk_size - 1, n_vox)
    Y_chunk <- Y_all[, start:end, drop = FALSE]
    res <- .parametric_engine(
      Y_proj = Y_chunk,
      S_target_proj = event_model,
      scan_times = seq_len(n_time),
      hrf_eval_times = hrf_interface$default_eval_times %||% seq(0, 30, length.out = 61),
      hrf_interface = hrf_interface,
      theta_seed = hrf_interface$default_seed(),
      theta_bounds = hrf_interface$default_bounds(),
      ...
    )
    theta_hat[start:end, ] <- res$theta_hat
    beta0[start:end] <- res$beta0
    r2[start:end] <- res$r_squared
  }

  list(theta_hat = theta_hat, beta0 = beta0, r_squared = r2)
}
