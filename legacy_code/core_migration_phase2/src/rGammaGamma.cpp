#include <RcppArmadillo.h>
#include "R_interface.h"
using namespace Rcpp;

namespace glmbayes {
namespace sim {


Rcpp::List rGammaGamma(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector beta,
    Rcpp::NumericVector wt,
    Rcpp::NumericVector alpha,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    bool verbose
) {
  // --------------------------------------------------------------
  // Stage 1: R functions we will use (for now)
  // --------------------------------------------------------------
  Rcpp::Function qgamma("qgamma");
  Rcpp::Function rgamma_ct = glmbayes_R::r_rgamma_ct();
  Rcpp::Function runif("runif");
  
  // Preserve whether user supplied bounds
  bool user_lower_null = disp_lower.isNull();
  bool user_upper_null = disp_upper.isNull();
  
  // --------------------------------------------------------------
  // Stage 2: Basic setup and Armadillo views
  // --------------------------------------------------------------
  int n_obs = y.size();
  
  arma::vec y_vec  = Rcpp::as<arma::vec>(y);
  arma::vec wt_vec = Rcpp::as<arma::vec>(wt);
  arma::vec a_vec  = Rcpp::as<arma::vec>(alpha);
  arma::vec b_vec  = Rcpp::as<arma::vec>(beta);
  arma::mat X      = Rcpp::as<arma::mat>(x);   // n_obs × p
  
  // mu = exp(alpha + X * beta)
  arma::vec linpred = a_vec + X * b_vec;
  arma::vec mu_vec  = arma::exp(linpred);
  
  double W_sum = arma::sum(wt_vec);
  
  // --------------------------------------------------------------
  // Stage 3: Log-likelihood and related functions
  // --------------------------------------------------------------
  auto loglik = [&](double v) -> double {
    double acc = 0.0;
    for (int i = 0; i < n_obs; ++i) {
      double yi  = y_vec[i];
      double mui = mu_vec[i];
      double wi  = wt_vec[i];
      
      acc += wi * (
        v * std::log(v) -
          v * std::log(mui) +
          (v - 1.0) * std::log(yi) -
          std::lgamma(v) -
          v * yi / mui
      );
    }
    return acc;
  };
  
  auto loglik_star = [&](double v) -> double {
    return loglik(v) - 0.5 * W_sum * std::log(v);
  };
  
  // rate1 = rate + sum_i w_i [ y_i / mu_i - log(y_i / mu_i) - 1 ]
  double rate1 = rate;
  for (int i = 0; i < n_obs; ++i) {
    double yi    = y_vec[i];
    double mui   = mu_vec[i];
    double wi    = wt_vec[i];
    double ratio = yi / mui;
    rate1 += wi * (ratio - std::log(ratio) - 1.0);
  }
  
  // --------------------------------------------------------------
  // Stage 4: Posterior mode via fixed-point iteration
  // --------------------------------------------------------------
  double shape2_mode = shape;
  double vstar1      = shape2_mode / rate1;
  
  auto vout_mode = [&](double v) -> double {
    double sum_term = 0.0;
    for (int i = 0; i < n_obs; ++i) {
      double wi = wt_vec[i];
      double wv = wi * v;
      sum_term += wi * (R::digamma(wv) - std::log(wv));
    }
    return vstar1 - (v / rate1) * sum_term;
  };
  
  double vstar = vstar1;
  for (int j = 0; j < 20; ++j) {
    vstar = vout_mode(vstar);
  }
  
  // --------------------------------------------------------------
  // Stage 5: Exact second derivative f''(v) of log-posterior
  // --------------------------------------------------------------
  auto f2_exact = [&](double v) -> double {
    double term1 = W_sum * (1.0 / v);
    double term2 = -W_sum * R::trigamma(v);
    double term3 = -(shape - 1.0) / (v * v);
    return term1 + term2 + term3;
  };
  
  double f2_vstar = f2_exact(vstar);
  
  // --------------------------------------------------------------
  // Stage 6: Gamma surrogate posterior for v
  // --------------------------------------------------------------
  double alpha_bar = 1.0 - f2_vstar * vstar * vstar;
  double beta_bar  = -f2_vstar * vstar;
  
  double v_mean = alpha_bar / beta_bar;
  
  // Surrogate-based bounds for v using qgamma
  double v_min_sur, v_max_sur;
  
  {
    Rcpp::NumericVector vmin_vec = qgamma(
      Rcpp::_["p"]          = 1.0 - max_disp_perc,
      Rcpp::_["shape"]      = alpha_bar,
      Rcpp::_["rate"]       = beta_bar,
      Rcpp::_["lower.tail"] = true,
      Rcpp::_["log.p"]      = false
    );
    Rcpp::NumericVector vmax_vec = qgamma(
      Rcpp::_["p"]          = max_disp_perc,
      Rcpp::_["shape"]      = alpha_bar,
      Rcpp::_["rate"]       = beta_bar,
      Rcpp::_["lower.tail"] = true,
      Rcpp::_["log.p"]      = false
    );
    v_min_sur = vmin_vec[0];
    v_max_sur = vmax_vec[0];
  }
  
  // --------------------------------------------------------------
  // Stage 7: Dispersion bounds and conversion to precision bounds
  // --------------------------------------------------------------
  double disp_lower_val, disp_upper_val;
  
  if (disp_lower.isNull()) {
    disp_lower_val = 1.0 / v_max_sur;
  } else {
    disp_lower_val = Rcpp::as<double>(disp_lower);
  }
  
  if (disp_upper.isNull()) {
    disp_upper_val = 1.0 / v_min_sur;
  } else {
    disp_upper_val = Rcpp::as<double>(disp_upper);
  }
  
  if (disp_lower_val <= 0.0 || disp_upper_val <= 0.0) {
    Rcpp::stop("disp_lower and disp_upper must be strictly positive when supplied.");
  }
  
  if (disp_lower_val >= disp_upper_val) {
    Rcpp::stop("Final dispersion bounds invalid: disp_lower must be < disp_upper.");
  }
  
  // Convert dispersion bounds to precision bounds
  double v_min = 1.0 / disp_upper_val;
  double v_max = 1.0 / disp_lower_val;
  
  // --------------------------------------------------------------
  // Stage 8: Tangency point and gradient cbar
  // --------------------------------------------------------------
  double v_tangent  = v_mean;
  double ll_star_v0 = loglik_star(v_tangent);
  
  auto ell_prime_at = [&](double v) -> double {
    double acc = 0.0;
    for (int i = 0; i < n_obs; ++i) {
      double yi  = y_vec[i];
      double mui = mu_vec[i];
      double wi  = wt_vec[i];
      
      acc += wi * (
        std::log(v) + 1.0 -
          std::log(mui) +
          std::log(yi) -
          R::digamma(v) -
          yi / mui
      );
    }
    return acc;
  };
  
  auto ell_star_prime_at = [&](double v) -> double {
    // adjustment term omitted as in R code
    return ell_prime_at(v);
  };
  
  double cbar = -ell_star_prime_at(v_tangent);
  
  // --------------------------------------------------------------
  // Stage 9: Proposal Gamma(shape_prop, rate_prop)
  // --------------------------------------------------------------
  double shape_prop = shape + 0.5 * W_sum;
  double rate_prop  = rate - cbar;
  
  if (rate_prop <= 0.0) {
    Rcpp::stop("Gamma proposal rate_prop <= 0; check v_tangent and curvature diagnostics.");
  }
  
  // --------------------------------------------------------------
  // Stage 10: Bounding function log h(v)
  // --------------------------------------------------------------
  auto lg_h = [&](double v) -> double {
    return loglik_star(v) -
      (ll_star_v0 - cbar * (v - v_tangent)) -
      0.5 * W_sum * std::log(v / v_min);
  };
  
  // --------------------------------------------------------------
  // Stage 11: Rejection sampling for precision v
  // --------------------------------------------------------------
  Rcpp::NumericVector dispersion(n);
  Rcpp::IntegerVector draws(n);
  
  for (int i = 0; i < n; ++i) {
    int draw_count = 1;
    bool accepted  = false;
    double v_curr  = NA_REAL;
    
    while (!accepted) {
      // Draw candidate precision from truncated Gamma via R's rgamma_ct
      Rcpp::RObject lower_prec_arg = R_NilValue;
      Rcpp::RObject upper_prec_arg = R_NilValue;
      
      if (user_lower_null && user_upper_null) {
        // use v_min, v_max
        lower_prec_arg = Rcpp::wrap(v_min);
        upper_prec_arg = Rcpp::wrap(v_max);
      } else {
        if (!user_upper_null) {
          lower_prec_arg = Rcpp::wrap(1.0 / disp_upper_val);
        }
        if (!user_lower_null) {
          upper_prec_arg = Rcpp::wrap(1.0 / disp_lower_val);
        }
      }
      
      Rcpp::NumericVector cand = rgamma_ct(
        Rcpp::_["n"]          = 1,
        Rcpp::_["shape"]      = shape_prop,
        Rcpp::_["rate"]       = rate_prop,
        Rcpp::_["lower_prec"] = lower_prec_arg,
        Rcpp::_["upper_prec"] = upper_prec_arg
      );
      
      v_curr = cand[0];
      
      // Accept–reject test
      Rcpp::NumericVector u_vec = runif(Rcpp::_["n"] = 1);
      double u = u_vec[0];
      
      double test_val = lg_h(v_curr) - std::log(u);
      
      if (test_val > 0.0) {
        accepted = true;
      } else {
        draw_count += 1;
      }
    }
    
    dispersion[i] = 1.0 / v_curr;  // convert precision to dispersion
    draws[i]      = draw_count;
  }
  
  // --------------------------------------------------------------
  // Final return
  // --------------------------------------------------------------
  return Rcpp::List::create(
    Rcpp::Named("dispersion") = dispersion,
    Rcpp::Named("draws")      = draws
  );
}

} // namespace sim
} // namespace glmbayes