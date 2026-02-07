#' Build reusable design inputs for parametric HRF estimation
#'
#' Prepares projected BOLD data and projected stimulus regressors once so they
#' can be reused across repeated model fits.
#'
#' @param fmri_data Matrix or dataset-like object accepted by
#'   [estimate_parametric_hrf()].
#' @param event_model Matrix or model object accepted by
#'   [estimate_parametric_hrf()].
#' @param confound_formula Optional confound formula.
#' @param baseline_model Optional baseline model.
#' @param hrf_eval_times Optional HRF evaluation grid.
#' @param hrf_span Span used when `hrf_eval_times` is `NULL`.
#' @param mask Optional logical mask selecting voxels.
#'
#' @return A list of class `parametric_design` with projected inputs and
#'   metadata.
#' @export
create_parametric_design <- function(
  fmri_data,
  event_model,
  confound_formula = NULL,
  baseline_model = NULL,
  hrf_eval_times = NULL,
  hrf_span = 30,
  mask = NULL
) {
  inputs <- .prepare_parametric_inputs(
    fmri_data = fmri_data,
    event_model = event_model,
    confound_formula = confound_formula,
    baseline_model = baseline_model,
    hrf_eval_times = hrf_eval_times,
    hrf_span = hrf_span,
    mask = mask
  )

  out <- list(
    Y_proj = inputs$Y_proj,
    S_target_proj = inputs$S_target_proj,
    hrf_eval_times = inputs$hrf_eval_times,
    Y_raw = inputs$Y_raw,
    S_target = inputs$S_target,
    Z = inputs$Z,
    scan_times = inputs$scan_times,
    baseline_model = inputs$baseline_model,
    metadata = list(
      n_time = nrow(inputs$Y_proj),
      n_vox = ncol(inputs$Y_proj),
      n_cond = ncol(inputs$S_target_proj),
      projected = !is.null(inputs$Z)
    )
  )

  class(out) <- c("parametric_design", "list")
  out
}

# Validate that a design list has required projected inputs.
.validate_parametric_design <- function(design) {
  required_fields <- c("Y_proj", "S_target_proj", "hrf_eval_times")
  missing_fields <- setdiff(required_fields, names(design))
  if (length(missing_fields) > 0) {
    stop(
      "'design' must contain: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.matrix(design$Y_proj) || !is.numeric(design$Y_proj)) {
    stop("design$Y_proj must be a numeric matrix", call. = FALSE)
  }
  if (!is.matrix(design$S_target_proj) || !is.numeric(design$S_target_proj)) {
    stop("design$S_target_proj must be a numeric matrix", call. = FALSE)
  }
  if (nrow(design$Y_proj) != nrow(design$S_target_proj)) {
    stop("design$Y_proj and design$S_target_proj must have same number of rows", call. = FALSE)
  }
  if (!is.numeric(design$hrf_eval_times) || length(design$hrf_eval_times) < 2) {
    stop("design$hrf_eval_times must be a numeric vector with length >= 2", call. = FALSE)
  }

  invisible(TRUE)
}

#' Estimate parametric HRFs from precomputed design inputs
#'
#' Interface-aligned entry point mirroring design-first workflows used in other
#' fMRI modeling packages.
#'
#' @param y Optional projected BOLD matrix (time x voxels). If `NULL`,
#'   `design$Y_proj` is used.
#' @param design A `parametric_design` object from [create_parametric_design()],
#'   or a compatible list containing `Y_proj`, `S_target_proj`, and
#'   `hrf_eval_times`.
#' @param parametric_model Parametric HRF model name.
#' @param ... Additional arguments forwarded to [estimate_parametric_hrf()].
#'
#' @return A `parametric_hrf_fit` object.
#' @export
estimate_parametric_hrf_from_design <- function(
  y = NULL,
  design,
  parametric_model = "lwu",
  ...
) {
  if (missing(design) || is.null(design) || !is.list(design)) {
    stop("'design' must be a non-null list", call. = FALSE)
  }
  .validate_parametric_design(design)

  if (is.null(y)) {
    y <- design$Y_proj
  }

  if (!is.matrix(y) || !is.numeric(y)) {
    stop("'y' must be a numeric matrix", call. = FALSE)
  }
  if (nrow(y) != nrow(design$S_target_proj)) {
    stop("'y' and design$S_target_proj must have the same number of rows", call. = FALSE)
  }

  prepared_inputs <- list(
    Y_raw = if (is.null(design$Y_raw)) y else design$Y_raw,
    Y_proj = y,
    S_target = if (is.null(design$S_target)) design$S_target_proj else design$S_target,
    S_target_proj = design$S_target_proj,
    Z = design$Z,
    scan_times = if (is.null(design$scan_times)) seq_len(nrow(y)) else design$scan_times,
    hrf_eval_times = design$hrf_eval_times,
    baseline_model = design$baseline_model
  )

  dots <- list(...)
  blocked <- intersect(
    names(dots),
    c(
      "fmri_data",
      "event_model",
      "parametric_model",
      "hrf_eval_times",
      "confound_formula",
      "baseline_model",
      "mask",
      "prepared_inputs"
    )
  )
  if (length(blocked) > 0) {
    warning(
      "Ignoring from_design argument(s): ",
      paste(blocked, collapse = ", "),
      call. = FALSE
    )
    dots[blocked] <- NULL
  }

  if (length(dots) == 0) {
    fit <- .run_hrf_estimation_engine(
      fmri_data = y,
      event_model = design$S_target_proj,
      parametric_model = parametric_model,
      hrf_eval_times = design$hrf_eval_times,
      confound_formula = NULL,
      baseline_model = NULL,
      prepared_inputs = prepared_inputs
    )
  } else {
    args <- c(
      list(
        fmri_data = y,
        event_model = design$S_target_proj,
        parametric_model = parametric_model,
        hrf_eval_times = design$hrf_eval_times,
        confound_formula = NULL,
        baseline_model = NULL,
        prepared_inputs = prepared_inputs
      ),
      dots
    )
    fit <- do.call(.run_hrf_estimation_engine, args)
  }

  fit$metadata$design_interface <- list(
    source = "from_design",
    design_class = class(design),
    projected = TRUE
  )

  fit
}
