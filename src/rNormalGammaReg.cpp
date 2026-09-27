#include "RcppArmadillo.h"
#include "R_interface.h"
using namespace Rcpp;


namespace glmbayes {
namespace sim {


Rcpp::List rNormalGammaReg(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector offset,
    Rcpp::NumericVector wt,
    double shape,
    double rate,
    Rcpp::Nullable<double> max_disp_perc,   // unused (kept for future extensions)
    Rcpp::Nullable<double> disp_lower,      // unused
    Rcpp::Nullable<double> disp_upper,      // unused
    bool verbose
) {
  // --------------------------------------------------------------
  // R helpers
  // --------------------------------------------------------------
  Rcpp::Function rNormal_reg_wfit = glmbayes_R::r_rNormal_reg_wfit();
  Rcpp::Function glmbfamfunc = glmbayes_R::r_glmbfamfunc();
  Rcpp::Function gaussian("gaussian");
  
  // famfunc = glmbfamfunc(gaussian())
  Rcpp::List famfunc = glmbfamfunc( gaussian() );
  
  // --------------------------------------------------------------
  // Call rNormal_reg.wfit(x, y, P, mu, w = wt, offset = offset, ...)
  // --------------------------------------------------------------
  Rcpp::List fit = rNormal_reg_wfit(
    Rcpp::_["x"]          = x,
    Rcpp::_["y"]          = y,
    Rcpp::_["P"]          = P,
    Rcpp::_["mu"]         = mu,
    Rcpp::_["w"]          = wt,
    Rcpp::_["offset"]     = offset,
    Rcpp::_["method"]     = "qr",
    Rcpp::_["tol"]        = 1e-7,
    Rcpp::_["singular.ok"]= true
  );
  
  Rcpp::NumericVector Btilde = fit["Btilde"];
  Rcpp::NumericMatrix IR     = fit["IR"];
  double S                   = Rcpp::as<double>(fit["S"]);
  int k                      = Rcpp::as<int>(fit["k"]);
  
  // --------------------------------------------------------------
  // Effective sample size
  // --------------------------------------------------------------
  double sum_wt = 0.0;
  for (int i = 0; i < wt.size(); ++i) {
    sum_wt += wt[i];
  }
  
  // --------------------------------------------------------------
  // Posterior Gamma parameters for precision tau = 1/phi
  // --------------------------------------------------------------
  double a_post = shape + sum_wt / 2.0;
  double b_post = rate  + 0.5 * S;
  
  if (a_post <= 0.0 || b_post <= 0.0) {
    Rcpp::stop("Invalid posterior Gamma parameters: check shape, rate, and S.");
  }
  
  // --------------------------------------------------------------
  // Allocate outputs
  // --------------------------------------------------------------
  Rcpp::NumericMatrix coef(n, k);
  Rcpp::NumericVector dispersion(n);
  Rcpp::IntegerVector draws(n);
  
  Rcpp::NumericVector z(k);
  Rcpp::NumericVector IRz(k);
  
  // --------------------------------------------------------------
  // Draw dispersion and coefficients
  // --------------------------------------------------------------
  for (int i = 0; i < n; ++i) {
    double tau = R::rgamma(a_post, 1.0 / b_post);
    double phi = 1.0 / tau;
    
    dispersion[i] = phi;
    draws[i]      = 1;
    
    double sd = std::sqrt(phi);
    
    for (int j = 0; j < k; ++j) {
      z[j] = R::rnorm(0.0, 1.0);
    }
    
    for (int r = 0; r < k; ++r) {
      double acc = 0.0;
      for (int c = 0; c < k; ++c) {
        acc += IR(r, c) * z[c];
      }
      IRz[r] = acc;
    }
    
    for (int j = 0; j < k; ++j) {
      coef(i, j) = Btilde[j] + IRz[j] * sd;
    }
  }
  
  // --------------------------------------------------------------
  // Prior list
  // --------------------------------------------------------------
  Rcpp::List Prior = Rcpp::List::create(
    Rcpp::Named("mean")      = Rcpp::as<Rcpp::NumericVector>(mu),
    Rcpp::Named("Precision") = P
  );
  
  Rcpp::List Envelope = Rcpp::List::create();
  
  // --------------------------------------------------------------
  // Final return
  // --------------------------------------------------------------
  return Rcpp::List::create(
    Rcpp::Named("coefficients")  = coef,
    Rcpp::Named("coef.mode")     = Btilde,
    Rcpp::Named("dispersion")    = dispersion,
    Rcpp::Named("offset")        = offset,
    Rcpp::Named("Prior")         = Prior,
    Rcpp::Named("prior.weights") = wt,
    Rcpp::Named("y")             = y,
    Rcpp::Named("x")             = x,
    Rcpp::Named("fit")           = fit,
    Rcpp::Named("famfunc")       = famfunc,
    Rcpp::Named("iters")         = draws,
    Rcpp::Named("Envelope")      = Envelope
  );
}


} // namespace sim
} // namespace glmbayes