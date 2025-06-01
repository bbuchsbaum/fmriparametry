#' Fast ridge regression solver using RcppEigen and OpenMP
#'
#' @param X Design matrix (time x parameters)
#' @param Y Response matrix (time x voxels)
#' @param lambda_ridge Ridge penalty
#' @keywords internal
.ridge_linear_solve <- function(X, Y, lambda_ridge = 0) {
  ridge_linear_solve_cpp(X, Y, lambda_ridge)
}
