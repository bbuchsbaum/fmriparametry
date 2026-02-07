#include <Rcpp.h>

using namespace Rcpp;

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

  NumericVector r_squared(n_vox);
  NumericVector rss(n_vox);
  NumericVector tss(n_vox);
  NumericVector mse(n_vox);
  NumericVector rmse(n_vox);
  NumericVector mae(n_vox);

  bool use_precomputed_tss = precomputed_tss.isNotNull();
  NumericVector tss_input;
  if (use_precomputed_tss) {
    tss_input = NumericVector(precomputed_tss);
  }

  for (int v = 0; v < n_vox; ++v) {
    double sum_y = 0.0;
    for (int i = 0; i < n_obs; ++i) {
      sum_y += y_true(i, v);
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

    double r2_v = 0.0;
    if (std::abs(tss_v) < tolerance) {
      if (std::abs(rss_v) < tolerance) {
        r2_v = 1.0;
      } else {
        r2_v = 0.0;
      }
    } else {
      r2_v = 1.0 - (rss_v / tss_v);
      if (r2_v < 0.0) r2_v = 0.0;
      if (r2_v > 1.0) r2_v = 1.0;
    }

    rss[v] = rss_v;
    tss[v] = tss_v;
    r_squared[v] = r2_v;
    mse[v] = rss_v / static_cast<double>(n_obs);
    rmse[v] = std::sqrt(mse[v]);
    mae[v] = mae_v / static_cast<double>(n_obs);
  }

  return List::create(
    _["r_squared"] = r_squared,
    _["rss"] = rss,
    _["tss"] = tss,
    _["mse"] = mse,
    _["rmse"] = rmse,
    _["mae"] = mae,
    _["n_obs"] = n_obs
  );
}
