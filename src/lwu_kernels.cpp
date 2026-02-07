#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>

using namespace Rcpp;

namespace {

inline double clamp_scalar(double x, double lo, double hi) {
  return std::min(std::max(x, lo), hi);
}

inline bool finite_vector(const NumericVector& x) {
  const int n = x.size();
  for (int i = 0; i < n; ++i) {
    if (!R_finite(x[i])) {
      return false;
    }
  }
  return true;
}

inline NumericVector lwu_formula_impl(
  const NumericVector& t,
  const double tau,
  const double sigma,
  const double rho,
  const std::string& normalize
) {
  const int n = t.size();
  NumericVector response(n);

  const double sigma2 = sigma * sigma;
  const double sigma_u = 1.6 * sigma;
  const double sigma_u2 = sigma_u * sigma_u;

  const double c1 = 2.0 * sigma2;
  const double c2 = 2.0 * sigma_u2;

  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    const double term1 = std::exp(-((ti - tau) * (ti - tau)) / c1);
    const double shifted = tau + 2.0 * sigma;
    const double term2 = rho * std::exp(-((ti - shifted) * (ti - shifted)) / c2);
    response[i] = term1 - term2;
  }

  if (normalize == "height") {
    double max_abs_val = 0.0;
    for (int i = 0; i < n; ++i) {
      const double abs_val = std::abs(response[i]);
      if (abs_val > max_abs_val) {
        max_abs_val = abs_val;
      }
    }
    if (max_abs_val > 1e-10) {
      for (int i = 0; i < n; ++i) {
        response[i] /= max_abs_val;
      }
    }
  }

  return response;
}

inline NumericVector convolve_signal_single(
  const NumericVector& signal,
  const NumericVector& kernel,
  const int output_length
) {
  const int signal_length = signal.size();
  const int kernel_length = kernel.size();
  NumericVector result(output_length);

  for (int i = 0; i < output_length; ++i) {
    double sum = 0.0;
    const int k_max = std::min(kernel_length - 1, i);
    for (int k = 0; k <= k_max; ++k) {
      const int idx = i - k;
      if (idx < signal_length) {
        sum += signal[idx] * kernel[k];
      }
    }
    result[i] = sum;
  }

  return result;
}

inline NumericMatrix convolve_signal_matrix(
  const NumericVector& signal,
  const NumericMatrix& kernels,
  const int output_length
) {
  const int signal_length = signal.size();
  const int kernel_length = kernels.nrow();
  const int n_kernels = kernels.ncol();

  NumericMatrix result(output_length, n_kernels);

  for (int j = 0; j < n_kernels; ++j) {
    for (int i = 0; i < output_length; ++i) {
      double sum = 0.0;
      const int k_max = std::min(kernel_length - 1, i);
      for (int k = 0; k <= k_max; ++k) {
        const int idx = i - k;
        if (idx < signal_length) {
          sum += signal[idx] * kernels(k, j);
        }
      }
      result(i, j) = sum;
    }
  }

  return result;
}

inline NumericVector clamp_lwu_theta(
  const NumericVector& params_vector0,
  const NumericVector& lower,
  const NumericVector& upper
) {
  NumericVector theta0(3);
  for (int k = 0; k < 3; ++k) {
    theta0[k] = clamp_scalar(params_vector0[k], lower[k] + 1e-6, upper[k] - 1e-6);
  }
  theta0[1] = std::max(theta0[1], 0.051);
  return theta0;
}

inline NumericMatrix lwu_taylor_basis_impl(
  const NumericVector& params_vector0,
  const NumericVector& t_hrf_eval,
  const NumericVector& lower,
  const NumericVector& upper,
  const double rel_step,
  const double min_step,
  const double bound_eps
) {
  NumericVector theta0 = clamp_lwu_theta(params_vector0, lower, upper);
  NumericVector base_hrf = lwu_formula_impl(t_hrf_eval, theta0[0], theta0[1], theta0[2], "none");

  const int n_time = t_hrf_eval.size();
  const int n_params = 3;
  NumericMatrix basis(n_time, n_params + 1);

  for (int i = 0; i < n_time; ++i) {
    basis(i, 0) = base_hrf[i];
  }

  const double step_floor[3] = {1.0, 0.2, 0.1};

  for (int k = 0; k < n_params; ++k) {
    NumericVector theta_plus = clone(theta0);
    NumericVector theta_minus = clone(theta0);

    const double step_scale = std::max(std::abs(theta0[k]), step_floor[k]);
    const double step_size = std::max(rel_step * step_scale, min_step);

    theta_plus[k] = std::min(theta0[k] + step_size, upper[k] - bound_eps);
    theta_minus[k] = std::max(theta0[k] - step_size, lower[k] + bound_eps);

    if (k == 1) {
      theta_plus[1] = std::max(theta_plus[1], 0.051);
      theta_minus[1] = std::max(theta_minus[1], 0.051);
    }

    const double denom = theta_plus[k] - theta_minus[k];

    if (denom > 1e-10) {
      NumericVector h_plus = lwu_formula_impl(t_hrf_eval, theta_plus[0], theta_plus[1], theta_plus[2], "none");
      NumericVector h_minus = lwu_formula_impl(t_hrf_eval, theta_minus[0], theta_minus[1], theta_minus[2], "none");
      for (int i = 0; i < n_time; ++i) {
        basis(i, k + 1) = (h_plus[i] - h_minus[i]) / denom;
      }
    } else if (theta_plus[k] > theta0[k] + 1e-10) {
      const double fwd_denom = theta_plus[k] - theta0[k];
      NumericVector h_plus = lwu_formula_impl(t_hrf_eval, theta_plus[0], theta_plus[1], theta_plus[2], "none");
      for (int i = 0; i < n_time; ++i) {
        basis(i, k + 1) = (h_plus[i] - base_hrf[i]) / fwd_denom;
      }
    } else if (theta_minus[k] < theta0[k] - 1e-10) {
      const double bwd_denom = theta0[k] - theta_minus[k];
      NumericVector h_minus = lwu_formula_impl(t_hrf_eval, theta_minus[0], theta_minus[1], theta_minus[2], "none");
      for (int i = 0; i < n_time; ++i) {
        basis(i, k + 1) = (base_hrf[i] - h_minus[i]) / bwd_denom;
      }
    } else {
      for (int i = 0; i < n_time; ++i) {
        basis(i, k + 1) = 0.0;
      }
    }
  }

  return basis;
}

} // anonymous namespace

// [[Rcpp::export]]
NumericVector lwu_hrf_formula_cpp(
  const NumericVector& t,
  const double tau,
  const double sigma,
  const double rho,
  std::string normalize = "none"
) {
  if (sigma <= 0.0 || !R_finite(sigma)) {
    stop("sigma must be positive and finite");
  }
  if (normalize != "none" && normalize != "height" && normalize != "area") {
    stop("normalize must be one of 'none', 'height', or 'area'");
  }

  // Keep parity with R wrapper behavior: treat 'area' as 'none'.
  if (normalize == "area") {
    normalize = "none";
  }

  return lwu_formula_impl(t, tau, sigma, rho, normalize);
}

// [[Rcpp::export]]
NumericMatrix lwu_taylor_basis_fd_cpp(
  const NumericVector& params_vector0,
  const NumericVector& t_hrf_eval,
  const NumericVector& lower,
  const NumericVector& upper,
  const double rel_step = 1e-4,
  const double min_step = 1e-5,
  const double bound_eps = 1e-8
) {
  if (params_vector0.size() != 3) {
    stop("params_vector0 must have length 3");
  }
  if (lower.size() != 3 || upper.size() != 3) {
    stop("lower and upper must have length 3");
  }

  return lwu_taylor_basis_impl(
    params_vector0,
    t_hrf_eval,
    lower,
    upper,
    rel_step,
    min_step,
    bound_eps
  );
}

// [[Rcpp::export]]
double lwu_gn_objective_cpp(
  const NumericVector& theta,
  const NumericVector& y,
  const NumericVector& signal,
  const NumericVector& t_hrf_eval,
  const NumericVector& lower,
  const NumericVector& upper
) {
  if (theta.size() != 3 || lower.size() != 3 || upper.size() != 3) {
    return R_PosInf;
  }
  if (y.size() == 0 || signal.size() == 0 || t_hrf_eval.size() == 0) {
    return R_PosInf;
  }
  if (!finite_vector(theta) || !finite_vector(y) || !finite_vector(signal) ||
      !finite_vector(t_hrf_eval) || !finite_vector(lower) || !finite_vector(upper)) {
    return R_PosInf;
  }

  NumericVector theta0 = clamp_lwu_theta(theta, lower, upper);
  NumericVector hrf_vals = lwu_formula_impl(t_hrf_eval, theta0[0], theta0[1], theta0[2], "none");
  NumericVector x_pred_raw = convolve_signal_single(signal, hrf_vals, y.size());

  double denom = 0.0;
  const int n_time = y.size();
  for (int i = 0; i < n_time; ++i) {
    denom += x_pred_raw[i] * x_pred_raw[i];
  }
  if (!R_finite(denom) || denom < 1e-8) {
    return R_PosInf;
  }

  double numer = 0.0;
  for (int i = 0; i < n_time; ++i) {
    numer += x_pred_raw[i] * y[i];
  }
  const double beta = numer / denom;
  if (!R_finite(beta)) {
    return R_PosInf;
  }

  double objective = 0.0;
  for (int i = 0; i < n_time; ++i) {
    const double residual = y[i] - beta * x_pred_raw[i];
    objective += residual * residual;
  }

  if (!R_finite(objective)) {
    return R_PosInf;
  }
  return objective;
}

// [[Rcpp::export]]
List lwu_gn_jacobian_cpp(
  const NumericVector& theta,
  const NumericVector& y,
  const NumericVector& signal,
  const NumericVector& t_hrf_eval,
  const NumericVector& lower,
  const NumericVector& upper,
  const double rel_step = 1e-4,
  const double min_step = 1e-5,
  const double bound_eps = 1e-8
) {
  if (theta.size() != 3 || lower.size() != 3 || upper.size() != 3) {
    return List::create(Named("ok") = false);
  }
  if (y.size() == 0 || signal.size() == 0 || t_hrf_eval.size() == 0) {
    return List::create(Named("ok") = false);
  }
  if (!finite_vector(theta) || !finite_vector(y) || !finite_vector(signal) ||
      !finite_vector(t_hrf_eval) || !finite_vector(lower) || !finite_vector(upper)) {
    return List::create(Named("ok") = false);
  }

  NumericMatrix basis = lwu_taylor_basis_impl(
    theta,
    t_hrf_eval,
    lower,
    upper,
    rel_step,
    min_step,
    bound_eps
  );
  NumericMatrix X_conv = convolve_signal_matrix(signal, basis, y.size());

  const int n_time = y.size();
  const int n_params = 3;

  NumericVector x_hrf(n_time);
  for (int i = 0; i < n_time; ++i) {
    x_hrf[i] = X_conv(i, 0);
  }

  double denom = 0.0;
  for (int i = 0; i < n_time; ++i) {
    denom += x_hrf[i] * x_hrf[i];
  }
  if (!R_finite(denom) || denom < 1e-8) {
    return List::create(Named("ok") = false);
  }

  double numer = 0.0;
  for (int i = 0; i < n_time; ++i) {
    numer += x_hrf[i] * y[i];
  }
  const double beta = numer / denom;
  if (!R_finite(beta)) {
    return List::create(Named("ok") = false);
  }

  NumericVector residuals(n_time);
  double objective = 0.0;
  for (int i = 0; i < n_time; ++i) {
    const double r = y[i] - beta * x_hrf[i];
    if (!R_finite(r)) {
      return List::create(Named("ok") = false);
    }
    residuals[i] = r;
    objective += r * r;
  }
  if (!R_finite(objective)) {
    return List::create(Named("ok") = false);
  }

  double dot_dy[3] = {0.0, 0.0, 0.0};
  double dot_xd[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < n_params; ++k) {
    for (int i = 0; i < n_time; ++i) {
      const double dx = X_conv(i, k + 1);
      dot_dy[k] += dx * y[i];
      dot_xd[k] += x_hrf[i] * dx;
    }
  }

  NumericVector dbeta(n_params);
  for (int k = 0; k < n_params; ++k) {
    dbeta[k] = (dot_dy[k] - 2.0 * beta * dot_xd[k]) / denom;
    if (!R_finite(dbeta[k])) {
      return List::create(Named("ok") = false);
    }
  }

  NumericMatrix jacobian(n_time, n_params);
  for (int k = 0; k < n_params; ++k) {
    for (int i = 0; i < n_time; ++i) {
      const double dx = X_conv(i, k + 1);
      const double jval = -beta * dx - x_hrf[i] * dbeta[k];
      if (!R_finite(jval)) {
        return List::create(Named("ok") = false);
      }
      jacobian(i, k) = jval;
    }
  }

  return List::create(
    Named("ok") = true,
    Named("objective") = objective,
    Named("jacobian") = jacobian,
    Named("residuals") = residuals,
    Named("amplitude") = beta
  );
}
