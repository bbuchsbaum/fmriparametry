#' @section Bounds Enforcement in Parametric HRF Fitting:
#'
#' Parameter bounds serve three critical purposes in the fitting algorithm:
#'
#' \subsection{1. Physiological Plausibility}{
#' The HRF parameters must be physiologically meaningful:
#' \itemize{
#'   \item \code{tau} (time-to-peak): Must be positive, typically 0-20 seconds
#'   \item \code{sigma} (width): Must be positive, typically 0.1-10 seconds
#'   \item \code{rho} (undershoot ratio): Typically 0-1.5
#' }
#' Bounds prevent the optimizer from finding mathematically valid but
#' physiologically nonsensical solutions.
#' }
#'
#' \subsection{2. Numerical Stability}{
#' Certain parameter values can cause numerical errors in HRF computation:
#' \itemize{
#'   \item \code{sigma} must be > 0.05 to avoid singularities in the underlying
#'     gamma functions
#'   \item Extreme parameter values can cause overflow/underflow in convolution
#'   \item Taylor expansion requires parameters to be sufficiently far from
#'     bounds for accurate numerical derivatives
#' }
#' }
#'
#' \subsection{3. Regularization and Identifiability}{
#' Bounds act as a strong form of regularization:
#' \itemize{
#'   \item Constrains the solution space, improving convergence
#'   \item Prevents equivalent fits in unrealistic parameter regions
#'   \item Helps with identifiability when data is noisy or sparse
#' }
#' }
#'
#' \subsection{Implementation Details}{
#' Bounds are enforced using a simple projection (clamping) method at multiple
#' stages:
#' \itemize{
#'   \item Initial parameter seeds are clamped to bounds
#'   \item After each Taylor approximation update
#'   \item During line search in Gauss-Newton refinement
#'   \item Before HRF evaluation to ensure numerical stability
#' }
#'
#' The clamping operation is: \code{theta_bounded = pmax(lower, pmin(theta, upper))}
#'
#' Note: This is a form of projected gradient descent. While simple and robust,
#' it can cause the algorithm to "stick" at boundaries. More sophisticated
#' constrained optimization methods (e.g., interior point methods) could be
#' used for better boundary behavior but at increased computational cost.
#' }