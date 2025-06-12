#' Calculate comprehensive fit metrics
#'
#' Calculates R-squared and other regression fit metrics with proper handling
#' of centered/uncentered models and pre-residualized data.
#'
#' @param y_true Numeric matrix or vector of true values (n_time x n_voxels)
#' @param y_pred Numeric matrix or vector of predicted values (same dimensions as y_true)
#' @param n_predictors Integer number of predictors (excluding intercept if present)
#' @param has_intercept Logical indicating if model includes an intercept term
#' @param precomputed_tss Optional numeric vector of pre-computed total sum of squares
#'   per voxel. Useful when y_true represents residuals from a previous stage.
#' @param tolerance Numeric tolerance for zero comparisons (default: 1e-10)
#'
#' @return List containing:
#'   \item{r_squared}{Numeric vector of R-squared values per voxel}
#'   \item{r_squared_adj}{Adjusted R-squared accounting for number of predictors}
#'   \item{mse}{Mean squared error}
#'   \item{rmse}{Root mean squared error}
#'   \item{mae}{Mean absolute error}
#'   \item{rss}{Residual sum of squares}
#'   \item{tss}{Total sum of squares}
#'   \item{n_obs}{Number of observations}
#'
#' @details
#' This function consolidates R-squared calculation logic across the package.
#' It properly handles:
#' - Models with and without intercepts
#' - Pre-residualized data (via precomputed_tss)
#' - Edge cases like zero variance
#' - Numerical stability issues
#'
#' @examples
#' # Simple example with intercept
#' y <- matrix(rnorm(100), ncol = 2)
#' y_pred <- y + rnorm(100, sd = 0.1)
#' metrics <- .calculate_fit_metrics(y, y_pred, n_predictors = 3, has_intercept = TRUE)
#' 
#' @keywords internal
.calculate_fit_metrics <- function(
  y_true,
  y_pred,
  n_predictors,
  has_intercept = TRUE,
  precomputed_tss = NULL,
  tolerance = 1e-10
) {
  # Convert to matrices if needed
  if (is.null(dim(y_true))) {
    y_true <- matrix(y_true, ncol = 1)
  }
  if (is.null(dim(y_pred))) {
    y_pred <- matrix(y_pred, ncol = 1)
  }
  
  # Input validation
  if (nrow(y_true) == 0 || ncol(y_true) == 0) {
    stop("Input arrays cannot be empty", call. = FALSE)
  }
  if (!identical(dim(y_true), dim(y_pred))) {
    stop(sprintf("Shape mismatch: y_true [%s] vs y_pred [%s]",
                 paste(dim(y_true), collapse = "x"),
                 paste(dim(y_pred), collapse = "x")), call. = FALSE)
  }
  if (!is.numeric(n_predictors) || n_predictors < 0) {
    stop("n_predictors must be a non-negative integer", call. = FALSE)
  }
  
  n_obs <- nrow(y_true)
  n_vox <- ncol(y_true)
  
  # Calculate residuals and RSS
  residuals <- y_true - y_pred
  rss <- colSums(residuals^2)
  
  # Calculate TSS
  if (!is.null(precomputed_tss)) {
    # Use pre-computed TSS (e.g., from original data before residualization)
    if (length(precomputed_tss) != n_vox) {
      stop(sprintf("precomputed_tss length (%d) must match number of voxels (%d)",
                   length(precomputed_tss), n_vox), call. = FALSE)
    }
    tss <- precomputed_tss
  } else {
    if (has_intercept) {
      # Standard R²: variance around the mean
      y_mean <- colMeans(y_true)
      tss <- colSums(sweep(y_true, 2, y_mean)^2)
    } else {
      # Uncentered R²: variance around zero
      tss <- colSums(y_true^2)
    }
  }
  
  # Calculate R-squared with edge case handling
  r_squared <- numeric(n_vox)
  for (v in seq_len(n_vox)) {
    if (abs(tss[v]) < tolerance) {
      # Zero variance in y_true
      if (abs(rss[v]) < tolerance) {
        # Perfect prediction of constant
        r_squared[v] <- 1.0
      } else {
        # Failed to predict constant
        r_squared[v] <- 0.0  # More conservative than -Inf for R
      }
    } else {
      r_squared[v] <- 1.0 - (rss[v] / tss[v])
      # Clamp to [0, 1] to handle numerical issues
      r_squared[v] <- max(0, min(1, r_squared[v]))
    }
  }
  
  # Calculate adjusted R-squared
  # Account for intercept in degrees of freedom if present
  p_total <- n_predictors + ifelse(has_intercept, 1, 0)
  adj_r2_denom <- n_obs - p_total
  
  if (adj_r2_denom <= 0) {
    r_squared_adj <- rep(NA_real_, n_vox)
  } else {
    r_squared_adj <- 1 - (1 - r_squared) * (n_obs - 1) / adj_r2_denom
  }
  
  # Calculate other metrics
  mse <- rss / n_obs
  rmse <- sqrt(mse)
  mae <- colMeans(abs(residuals))
  
  # Return comprehensive metrics
  list(
    r_squared = r_squared,
    r_squared_adj = r_squared_adj,
    mse = mse,
    rmse = rmse,
    mae = mae,
    rss = rss,
    tss = tss,
    n_obs = n_obs,
    n_predictors = n_predictors,
    has_intercept = has_intercept
  )
}

#' Simplified R-squared calculation
#'
#' Wrapper around .calculate_fit_metrics for backward compatibility
#' and simple use cases where only R-squared is needed.
#'
#' @param Y Numeric matrix of true values
#' @param Y_pred Numeric matrix of predicted values
#' @param has_intercept Logical, does the model include an intercept?
#'
#' @return Numeric vector of R-squared values
#' @keywords internal
.compute_r_squared <- function(Y, Y_pred, has_intercept = TRUE) {
  metrics <- .calculate_fit_metrics(
    y_true = Y,
    y_pred = Y_pred,
    n_predictors = 0,  # Not used for basic R² calculation
    has_intercept = has_intercept
  )
  metrics$r_squared
}