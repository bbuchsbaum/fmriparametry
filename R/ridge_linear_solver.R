#' Fast ridge regression solver using RcppEigen and OpenMP
#'
#' @param X Design matrix (time x parameters)
#' @param Y Response matrix (time x voxels)
#' @param lambda_ridge Non-negative numeric ridge penalty (scalar)
#' @keywords internal
.ridge_linear_solve <- function(X, Y, lambda_ridge = 0) {
  # Ensure inputs are matrices
  if (!is.matrix(X)) {
    stop("X must be a matrix in .ridge_linear_solve")
  }
  if (!is.matrix(Y)) {
    stop("Y must be a matrix in .ridge_linear_solve")
  }
  
  # Check for NULL or empty inputs
  if (is.null(X) || length(X) == 0) {
    stop("X is NULL or empty in .ridge_linear_solve")
  }
  if (is.null(Y) || length(Y) == 0) {
    stop("Y is NULL or empty in .ridge_linear_solve")
  }
  

  # Ensure lambda is numeric
  if (!is.numeric(lambda_ridge) || length(lambda_ridge) != 1) {
    stop("lambda_ridge must be a single numeric value")
  }

  # Dimension checks
  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows")
  }

  # Non-finite value checks
  if (any(!is.finite(X))) stop("X contains non-finite values")
  if (any(!is.finite(Y))) stop("Y contains non-finite values")


  # Validate ridge penalty
  lambda_ridge <- .validate_numeric_param(
    lambda_ridge, "lambda_ridge",
    min_val = 0, allow_null = FALSE,
    caller = ".ridge_linear_solve"
  )
  

  ridge_linear_solve_cpp(X, Y, lambda_ridge)
}
