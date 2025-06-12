#' Prepare data and design matrices for parametric fitting
#'
#' Internal helper that extracts BOLD data, constructs stimulus regressors
#' and projects out confounds.
#'
#' @param fmri_data Matrix or dataset-like object
#' @param event_model Matrix or object convertible to design matrix
#' @param confound_formula Optional formula for confound regressors
#' @param baseline_model Baseline model specification (unused placeholder)
#' @param hrf_eval_times Optional vector of HRF evaluation times
#' @param hrf_span Span for HRF evaluation grid when `hrf_eval_times` is NULL
#' @param mask Optional logical vector selecting voxels
#'
#' @return List with projected data, design matrices and timing vectors
#' @keywords internal
.prepare_parametric_inputs <- function(
  fmri_data,
  event_model,
  confound_formula = NULL,
  baseline_model = "intercept",
  hrf_eval_times = NULL,
  hrf_span = 30,
  mask = NULL
) {
  # Step 1: extract data matrix and scan times
  if (is.matrix(fmri_data)) {
    Y_raw <- fmri_data
    scan_times <- seq_len(nrow(Y_raw))
  } else if (requireNamespace("fmrireg", quietly = TRUE) &&
             inherits(fmri_data, c("fmri_dataset", "matrix_dataset"))) {
    Y_raw <- fmrireg::get_data_matrix(fmri_data, mask = mask)
    # Extract scan times - matrix_dataset doesn't have a sampling_frame
    # so we'll use the number of rows
    scan_times <- seq_len(nrow(Y_raw))
  } else if (is.list(fmri_data) && !is.null(fmri_data$data)) {
    Y_raw <- fmri_data$data
    scan_times <- fmri_data$scan_times
    if (is.null(scan_times)) {
      scan_times <- seq_len(nrow(Y_raw))
    }
    if (!is.null(mask) && is.logical(mask)) {
      Y_raw <- Y_raw[, mask, drop = FALSE]
    }
  } else {
    stop("Unsupported fmri_data type", call. = FALSE)
  }



  # Step 2: stimulus design matrix
  if (is.matrix(event_model)) {
    S_target <- event_model
  } else if (requireNamespace("fmrireg", quietly = TRUE) &&
             inherits(event_model, "event_model")) {
    S_target <- as.matrix(fmrireg::design_matrix(event_model))
  } else if (is.list(event_model) && !is.null(event_model$design_matrix)) {
    S_target <- event_model$design_matrix
  } else {
    stop("Unsupported event_model type", call. = FALSE)
  }
  if (nrow(S_target) != length(scan_times)) {
    stop("event_model design matrix has wrong number of rows", call. = FALSE)
  }
  
  # Check if design matrix has any events
  if (all(S_target == 0)) {
    warning("No events detected in event_model", call. = FALSE)
  }

  # Step 3: confound matrix
  if (is.null(confound_formula)) {
    Z <- NULL
  } else {
    df_tmp <- data.frame(scan = scan_times)
    Z <- stats::model.matrix(confound_formula, data = df_tmp)
  }

  # Step 4: project out confounds
  if (!is.null(Z)) {
    qr_Z <- qr(Z)
    Q <- qr.Q(qr_Z)
    Y_proj <- Y_raw - Q %*% (t(Q) %*% Y_raw)
    S_target_proj <- S_target - Q %*% (t(Q) %*% S_target)
  } else {
    Y_proj <- Y_raw
    S_target_proj <- S_target
  }

  # Step 5: HRF evaluation times
  if (is.null(hrf_eval_times)) {
    dt <- if (length(scan_times) > 1) median(diff(scan_times)) else 1
    hrf_eval_times <- seq(0, hrf_span, by = dt)
  }

  # Ensure outputs are double for C++ compatibility
  storage.mode(Y_proj) <- "double"
  storage.mode(S_target_proj) <- "double"
  storage.mode(S_target) <- "double"

  list(
    Y_raw = Y_raw,
    Y_proj = Y_proj,
    S_target = S_target,
    S_target_proj = S_target_proj,
    Z = Z,
    scan_times = scan_times,
    hrf_eval_times = hrf_eval_times,
    baseline_model = baseline_model
  )
}
