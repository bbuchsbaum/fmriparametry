#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
Rcpp::List calculate_fit_metrics_cpp(
  const NumericMatrix& y_true,
  const NumericMatrix& y_pred,
  const bool has_intercept = true,
  Nullable<NumericVector> precomputed_tss = R_NilValue,
  const double tolerance = 1e-10
) {
  const int n_obs = y_true.nrow();
  const int n_vox = y_true.ncol();
  if (n_obs <= 0 || n_vox <= 0) {
    Rcpp::stop("y_true must be non-empty");
  }
  if (y_pred.nrow() != n_obs || y_pred.ncol() != n_vox) {
    Rcpp::stop("y_true and y_pred must have identical dimensions");
  }

  std::vector<double> r_squared_vals(static_cast<size_t>(n_vox), NA_REAL);
  std::vector<double> r_squared_raw_vals(static_cast<size_t>(n_vox), NA_REAL);
  std::vector<double> rss_vals(static_cast<size_t>(n_vox), NA_REAL);
  std::vector<double> tss_vals(static_cast<size_t>(n_vox), NA_REAL);
  std::vector<double> mse_vals(static_cast<size_t>(n_vox), NA_REAL);
  std::vector<double> rmse_vals(static_cast<size_t>(n_vox), NA_REAL);
  std::vector<double> mae_vals(static_cast<size_t>(n_vox), NA_REAL);

  bool use_precomputed_tss = precomputed_tss.isNotNull();
  std::vector<double> tss_input;
  if (use_precomputed_tss) {
    NumericVector tss_r = NumericVector(precomputed_tss);
    if (tss_r.size() != n_vox) {
      Rcpp::stop("precomputed_tss length must match number of voxels");
    }
    tss_input.assign(tss_r.begin(), tss_r.end());
  }

  std::atomic<int> bad_voxel(-1);

#ifdef _OPENMP
#pragma omp parallel for if(n_vox >= 32)
#endif
  for (int v = 0; v < n_vox; ++v) {
    if (bad_voxel.load() >= 0) {
      continue;
    }

    double sum_y = 0.0;
    for (int i = 0; i < n_obs; ++i) {
      if (!R_finite(y_true(i, v)) || !R_finite(y_pred(i, v))) {
        int expected = -1;
        bad_voxel.compare_exchange_strong(expected, v);
        break;
      }
      sum_y += y_true(i, v);
    }
    if (bad_voxel.load() >= 0) {
      continue;
    }
    const double mean_y = sum_y / static_cast<double>(n_obs);

    double rss_v = 0.0;
    double mae_v = 0.0;
    double centered_tss_v = 0.0;
    double uncentered_tss_v = 0.0;

    for (int i = 0; i < n_obs; ++i) {
      const double y = y_true(i, v);
      const double yhat = y_pred(i, v);
      const double resid = y - yhat;
      const double centered = y - mean_y;

      rss_v += resid * resid;
      mae_v += std::abs(resid);
      centered_tss_v += centered * centered;
      uncentered_tss_v += y * y;
    }

    double tss_v = 0.0;
    if (use_precomputed_tss) {
      tss_v = tss_input[v];
    } else if (has_intercept) {
      tss_v = centered_tss_v;
    } else {
      tss_v = uncentered_tss_v;
      if (centered_tss_v < tolerance) {
        tss_v = centered_tss_v;
      }
    }

    double r2_raw_v = 0.0;
    if (std::abs(tss_v) < tolerance) {
      if (std::abs(rss_v) < tolerance) {
        r2_raw_v = 1.0;
      } else {
        r2_raw_v = 0.0;
      }
    } else {
      r2_raw_v = 1.0 - (rss_v / tss_v);
    }
    const double r2_v = std::min(1.0, std::max(0.0, r2_raw_v));

    rss_vals[static_cast<size_t>(v)] = rss_v;
    tss_vals[static_cast<size_t>(v)] = tss_v;
    r_squared_vals[static_cast<size_t>(v)] = r2_v;
    r_squared_raw_vals[static_cast<size_t>(v)] = r2_raw_v;
    mse_vals[static_cast<size_t>(v)] = rss_v / static_cast<double>(n_obs);
    rmse_vals[static_cast<size_t>(v)] = std::sqrt(mse_vals[static_cast<size_t>(v)]);
    mae_vals[static_cast<size_t>(v)] = mae_v / static_cast<double>(n_obs);
  }

  if (bad_voxel.load() >= 0) {
    Rcpp::stop("y_true and y_pred must contain only finite values");
  }

  NumericVector r_squared(n_vox);
  NumericVector r_squared_raw(n_vox);
  NumericVector rss(n_vox);
  NumericVector tss(n_vox);
  NumericVector mse(n_vox);
  NumericVector rmse(n_vox);
  NumericVector mae(n_vox);
  for (int v = 0; v < n_vox; ++v) {
    r_squared[v] = r_squared_vals[static_cast<size_t>(v)];
    r_squared_raw[v] = r_squared_raw_vals[static_cast<size_t>(v)];
    rss[v] = rss_vals[static_cast<size_t>(v)];
    tss[v] = tss_vals[static_cast<size_t>(v)];
    mse[v] = mse_vals[static_cast<size_t>(v)];
    rmse[v] = rmse_vals[static_cast<size_t>(v)];
    mae[v] = mae_vals[static_cast<size_t>(v)];
  }

  return List::create(
    _["r_squared"] = r_squared,
    _["r_squared_raw"] = r_squared_raw,
    _["rss"] = rss,
    _["tss"] = tss,
    _["mse"] = mse,
    _["rmse"] = rmse,
    _["mae"] = mae,
    _["n_obs"] = n_obs
  );
}
