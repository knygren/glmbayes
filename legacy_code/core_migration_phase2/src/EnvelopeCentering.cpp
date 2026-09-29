// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
#include "Envelopefuncs.h"
#include <cmath>

using namespace Rcpp;

namespace glmbayes {

namespace env {

List EnvelopeCentering(
    NumericVector y,
    NumericMatrix x,
    NumericVector mu,
    NumericMatrix P,
    NumericVector offset,
    NumericVector wt,
    double shape,
    double rate,
    int Gridtype,
    bool verbose
) {
  (void)verbose;
  (void)Gridtype;
  const int n_rss_iter = 10;
  Rcpp::Function lm_wfit("lm.wfit");

  int n_obs = y.size();
  NumericVector ystar(n_obs);
  for (int i = 0; i < n_obs; i++) {
    ystar[i] = y[i] - offset[i];
  }

  double n_w = 0.0;
  for (int i = 0; i < wt.size(); ++i) n_w += wt[i];

  Rcpp::List fit = lm_wfit(
    Rcpp::_["x"] = x,
    Rcpp::_["y"] = ystar,
    Rcpp::_["w"] = wt
  );

  NumericVector res = fit["residuals"];
  double RSS = 0.0;
  for (int i = 0; i < res.size(); i++) {
    RSS += res[i] * res[i];
  }
  int p = Rcpp::as<int>(fit["rank"]);
  double dispersion2 = RSS / (n_obs - p);

  // Precompute (fixed inputs): chol(P), sqrt(w)-scaled design/response, XtWX.
  const int l1 = x.ncol();
  const int l2 = x.nrow();
  const arma::mat X = as<arma::mat>(x);
  const arma::vec Y = as<arma::vec>(y);
  const arma::vec off = as<arma::vec>(offset);
  const arma::vec wv = as<arma::vec>(wt);
  const arma::vec mu_vec = as<arma::vec>(mu);
  arma::mat P_arma(P.begin(), P.nrow(), P.ncol(), false);

  arma::mat Xw(l2, l1);
  arma::vec yw(l2);
  for (int i = 0; i < l2; i++) {
    const double sw = std::sqrt(wv[i]);
    Xw.row(i) = sw * X.row(i);
    yw[i] = (Y[i] - off[i]) * sw;
  }

  const arma::mat RA = arma::chol(P_arma);
  const arma::vec z_bot = RA * mu_vec;
  const arma::mat XtWX = X.t() * (arma::diagmat(wv) * X);

  arma::mat W(l2 + l1, l1);
  arma::vec z(l2 + l1);

  double RSS_post_expected = NA_REAL;

  for (int j = 0; j < n_rss_iter; ++j) {
    const double s = 1.0 / std::sqrt(dispersion2);
    W.rows(0, l2 - 1) = s * Xw;
    W.rows(l2, l2 + l1 - 1) = RA;
    z.rows(0, l2 - 1) = s * yw;
    z.rows(l2, l2 + l1 - 1) = z_bot;

    const arma::mat WtW = W.t() * W;
    const arma::mat IR = arma::inv(arma::trimatu(arma::chol(WtW)));
    const arma::mat Sigma = IR * arma::trans(IR);
    const arma::vec b2_fast = Sigma * (W.t() * z);

    const arma::vec r_fast = Y - X * b2_fast - off;
    const double rss_at_mean_fast = arma::dot(wv, r_fast % r_fast);
    const double trace_term_fast = arma::trace(XtWX * Sigma);
    const double RSS_precomputed = rss_at_mean_fast + trace_term_fast;

    RSS_post_expected = RSS_precomputed;

    double shape2 = shape + n_w / 2.0;
    double rate2 = rate + RSS_precomputed / 2.0;
    dispersion2 = rate2 / (shape2 - 1.0);
  }

  return List::create(
    Named("dispersion") = dispersion2,
    Named("RSS_post") = RSS_post_expected
  );
}

}  // namespace env

}  // namespace glmbayes
