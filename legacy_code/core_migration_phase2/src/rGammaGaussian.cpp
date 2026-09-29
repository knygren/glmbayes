#include "RcppArmadillo.h"
#include "R_interface.h"
using namespace Rcpp;


namespace glmbayes {
namespace sim {

Rcpp::List rGammaGaussian(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector beta,
    Rcpp::NumericVector wt,
    Rcpp::NumericVector alpha,
    double shape,
    double rate,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    bool verbose
) {
  // --------------------------------------------------------------
  // R helper for truncated gamma
  // --------------------------------------------------------------
  Rcpp::Function rgamma_ct = glmbayes_R::r_rgamma_ct();
  
  int n_obs = y.size();
  int p     = x.ncol();
  
  if (x.nrow() != n_obs)
    Rcpp::stop("Number of rows in x must match length of y.");
  if (beta.size() != p)
    Rcpp::stop("Length of beta must match number of columns in x.");
  if (wt.size() != n_obs)
    Rcpp::stop("Length of wt must match length of y.");
  if (alpha.size() != n_obs)
    Rcpp::stop("Length of alpha must match length of y.");
  
  // --------------------------------------------------------------
  // y1 = y - alpha
  // --------------------------------------------------------------
  Rcpp::NumericVector y1(n_obs);
  for (int i = 0; i < n_obs; ++i)
    y1[i] = y[i] - alpha[i];
  
  // --------------------------------------------------------------
  // xb = X * beta
  // --------------------------------------------------------------
  Rcpp::NumericVector xb(n_obs);
  for (int i = 0; i < n_obs; ++i) {
    double acc = 0.0;
    for (int j = 0; j < p; ++j)
      acc += x(i, j) * beta[j];
    xb[i] = acc;
  }
  
  // --------------------------------------------------------------
  // residuals and weighted sum of squares
  // --------------------------------------------------------------
  double sum_wt = 0.0;
  double sum_wt_ss = 0.0;
  
  for (int i = 0; i < n_obs; ++i) {
    double r = y1[i] - xb[i];
    double ss = r * r;
    double wi = wt[i];
    sum_wt    += wi;
    sum_wt_ss += wi * ss;
  }
  
  // --------------------------------------------------------------
  // posterior parameters
  // --------------------------------------------------------------
  double a1 = shape + sum_wt / 2.0;
  double b1 = rate  + sum_wt_ss / 2.0;
  
  // --------------------------------------------------------------
  // output containers
  // --------------------------------------------------------------
  Rcpp::NumericVector dispersion(n);
  Rcpp::IntegerVector draws(n);
  
  bool lower_null = disp_lower.isNull();
  bool upper_null = disp_upper.isNull();
  
  // --------------------------------------------------------------
  // Case 1: no bounds → plain inverse-gamma
  // --------------------------------------------------------------
  if (lower_null && upper_null) {
    for (int i = 0; i < n; ++i) {
      // R::rgamma(shape, scale) with scale = 1/rate
      double v = R::rgamma(a1, 1.0 / b1);
      dispersion[i] = 1.0 / v;
      draws[i]      = 1;
    }
  }
  
  // --------------------------------------------------------------
  // Case 2: truncated via rgamma_ct
  // --------------------------------------------------------------
  else {
    double lower_prec = NA_REAL;
    double upper_prec = NA_REAL;
    
    if (!upper_null) {
      double du = Rcpp::as<double>(disp_upper);
      if (du <= 0.0)
        Rcpp::stop("disp_upper must be > 0.");
      lower_prec = 1.0 / du;
    }
    
    if (!lower_null) {
      double dl = Rcpp::as<double>(disp_lower);
      if (dl <= 0.0)
        Rcpp::stop("disp_lower must be > 0.");
      upper_prec = 1.0 / dl;
    }
    
    for (int i = 0; i < n; ++i) {
      Rcpp::NumericVector v_vec = rgamma_ct(
        Rcpp::_["n"]          = 1,
        Rcpp::_["shape"]      = a1,
        Rcpp::_["rate"]       = b1,
        Rcpp::_["lower_prec"] = Rcpp::NumericVector::create(lower_prec),
        Rcpp::_["upper_prec"] = Rcpp::NumericVector::create(upper_prec)
      );
      
      double v = v_vec[0];
      dispersion[i] = 1.0 / v;
      draws[i]      = 1;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("dispersion") = dispersion,
    Rcpp::Named("draws")      = draws
  );
}

} // namespace sim
} // namespace glmbayes