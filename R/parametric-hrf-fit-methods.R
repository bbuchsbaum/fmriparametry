#' Methods for parametric_hrf_fit objects
#'
#' These S3 methods provide a unified user interface for objects returned by
#' `estimate_parametric_hrf()` and related functions.
#'
#' @name parametric_hrf_fit-methods
NULL

#' Pretty printing for parametric_hrf_fit
#'
#' Provides a concise summary of the fit using `cli`/`pillar` when available.
#'
#' @param x A \code{parametric_hrf_fit} object.
#' @param ... Additional arguments (ignored).
#' @return The input object invisibly.
#' @export
print.parametric_hrf_fit <- function(x, ...) {
  have_cli <- requireNamespace("cli", quietly = TRUE) &&
    requireNamespace("pillar", quietly = TRUE)

  summ <- summary(x)

  if (have_cli) {
    cli::cli_h1("Parametric HRF Fit")
    cli::cli_text("Model: {.val {x$hrf_model}}")
    cli::cli_text("Voxels: {n_voxels(x)}")
    cli::cli_h2("Summary")
    print(summ)
  } else {
    cat("Parametric HRF Fit\n")
    cat("Model:", x$hrf_model, "\n")
    cat("Voxels:", n_voxels(x), "\n\n")
    print(summ)
  }
  invisible(x)
}

#' Coefficients from a parametric_hrf_fit
#'
#' @param object A \code{parametric_hrf_fit} object.
#' @param type Which coefficients to return: "parameters", "amplitude", or "se".
#' @param ... Additional arguments (ignored).
#' @return Numeric matrix or vector depending on \code{type}.
#' @export
coef.parametric_hrf_fit <- function(object,
                                   type = c("parameters", "amplitude", "se"),
                                   ...) {
  type <- match.arg(type)
  switch(type,
         parameters = object$estimated_parameters,
         amplitude = object$amplitudes,
         se = {
           if (is.null(object$parameter_ses)) {
             warning("Standard errors not available")
             NULL
           } else {
             object$parameter_ses
           }
         })
}

#' Summarize a parametric_hrf_fit
#'
#' Returns a tidy data.frame of summary statistics for each parameter,
#' amplitudes, and R-squared values when available.
#'
#' @param object A \code{parametric_hrf_fit} object.
#' @param ... Additional arguments (ignored).
#' @return A data.frame with summary statistics and attributes `hrf_model`
#'   and `n_voxels`.
#' @export
summary.parametric_hrf_fit <- function(object, ...) {
  stat_vec <- function(v) {
    qs <- stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    c(min = qs[1], q1 = qs[2], median = qs[3], mean = mean(v, na.rm = TRUE),
      q3 = qs[4], max = qs[5])
  }

  param_stats <- t(apply(object$estimated_parameters, 2, stat_vec))
  df <- data.frame(parameter = rownames(param_stats), param_stats,
                   row.names = NULL, check.names = FALSE)

  amp_stats <- stat_vec(object$amplitudes)
  df <- rbind(df, data.frame(parameter = "amplitude", t(amp_stats)))

  if (!is.null(object$r_squared)) {
    r2_stats <- stat_vec(object$r_squared)
    df <- rbind(df, data.frame(parameter = "r_squared", t(r2_stats)))
  }

  attr(df, "hrf_model") <- object$hrf_model
  attr(df, "n_voxels") <- n_voxels(object)
  class(df) <- c("summary_parametric_hrf_fit", "data.frame")
  df
}

#' @export
print.summary_parametric_hrf_fit <- function(x, ...) {
  have_cli <- requireNamespace("cli", quietly = TRUE) &&
    requireNamespace("pillar", quietly = TRUE)

  if (have_cli) {
    cli::cli_h1("Parametric HRF Summary")
    cli::cli_text("Model: {.val {attr(x, 'hrf_model')}}")
    cli::cli_text("Voxels: {attr(x, 'n_voxels')}")
    pillar::print_table(x)
  } else {
    cat("Parametric HRF Summary\n")
    cat("Model:", attr(x, "hrf_model"), "\n")
    cat("Voxels:", attr(x, "n_voxels"), "\n\n")
    print.data.frame(x, row.names = FALSE)
  }
  invisible(x)
}

#' @export
fitted.parametric_hrf_fit <- function(object, Y_proj = NULL, ...) {
  if (!is.null(object$residuals)) {
    if (is.null(Y_proj)) {
      stop("Y_proj required to compute fitted values from residuals")
    }
    return(Y_proj - object$residuals)
  }
  stop("Cannot compute fitted values without residuals")
}

#' @export
residuals.parametric_hrf_fit <- function(object, ...) {
  if (is.null(object$residuals)) {
    warning("Residuals not stored in fit object")
    return(NULL)
  }
  object$residuals
}
