// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include <RcppParallel.h>
#include "openclPort.h"
#include "progress_utils.h"
#include "rng_utils.h"  // for safe_runif()
#include <vector>
#include <algorithm>
#include <cmath>


using namespace Rcpp;
using namespace openclPort;
using namespace glmbayes::fam;
using namespace glmbayes::env;
using namespace glmbayes::rng;



NumericVector RSS(NumericVector y, NumericMatrix x,NumericMatrix b,NumericVector alpha,NumericVector wt)
{
  // Step 1: Set up dimensions
  
  int l1 = x.nrow(), l2 = x.ncol(); // Dimensions of x matrix (dims for y,alpha, and wt needs to be consistent) 
  int m1 = b.ncol();                // Number of columns for which output is needed
  
  // Step 2: Initialize b2temp and other Rcpp and arma objects used in calculations
  
  Rcpp::NumericMatrix b2temp(l2,1);
  Rcpp::NumericMatrix restemp(1,1);
  arma::mat y2(y.begin(), l1, 1, false);
  arma::mat x2(x.begin(), l1, l2, false); 
  arma::mat alpha2(alpha.begin(), l1, 1, false); 
  
  Rcpp::NumericVector xb(l1);
  arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
  
  NumericVector sqrt_wt=sqrt(wt);
  arma::mat sqrt_wt2(sqrt_wt.begin(), l1, 1, false); 
  
  //  NumericVector invwt=1/sqrt(wt);
  
  // Moving Loop inside the function is key for speed
  
  NumericVector yy(l1);
  NumericVector res(m1);
  arma::colvec res2(res.begin(),m1,false); // Reuse memory - update both below
  
  for(int i=0;i<m1;i++){
    
    // Grab one column at a time from b and one row at a time from res
    
    b2temp=b(Range(0,l2-1),Range(i,i));
    
    // Point b2 to memory for that column
    
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat restemp(res.begin()+i, 1, 1, false); 
    
    // calculate weighted residuals (element by element multiplication with weights)
    
    xb2=(y2-alpha2- x2 * b2)%sqrt_wt2;
    
    // This is where RSS should be calculated
    // Not sure if this will complain about type differences
    
    restemp=trans(xb2)*xb2;
    
  }
  
  return res;      
  
}




double rss_face_at_disp(double dispersion,
                               Rcpp::List cache,
                               Rcpp::NumericVector cbars_j,
                               Rcpp::NumericVector y,
                               Rcpp::NumericMatrix x,
                               Rcpp::NumericVector alpha,
                               Rcpp::NumericVector wt) {
  // Build 1×l1 matrix, then transpose to l1×1 for Inv_f3_with_disp
  int l1 = cbars_j.size();
  Rcpp::NumericMatrix cbars_small(1, l1);
  for (int k = 0; k < l1; ++k) cbars_small(0, k) = cbars_j[k];
  
  arma::mat theta_row = Inv_f3_with_disp(cache, dispersion, Rcpp::transpose(cbars_small));
  arma::vec beta = theta_row.t(); // 1×l1 -> l1×1
  
  arma::vec y2(y.begin(), y.size(), false);
  arma::vec a2(alpha.begin(), alpha.size(), false);
  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
  arma::vec w(wt.begin(), wt.size(), false);
  
  arma::vec resid = (y2 - a2 - X * beta) % arma::sqrt(w);
  return arma::as_scalar(resid.t() * resid);
}



double UB2(double dispersion,
           Rcpp::List cache,
           Rcpp::NumericVector cbars_j,
           Rcpp::NumericVector y,
           Rcpp::NumericMatrix x,
           Rcpp::NumericVector alpha,
           Rcpp::NumericVector wt,
           double rss_min_global) {
  
  // Call the existing RSS function
  double rss_val = rss_face_at_disp(dispersion, cache, cbars_j, y, x, alpha, wt);
  
  // Compute UB2
  double UB2_val = (0.5 / dispersion) * (rss_val - rss_min_global);
  
  return UB2_val;
}


// Utility: safe max for NumericVector
static inline double max_vec(const NumericVector& v) {
  double m = R_NegInf;
  for (int i = 0; i < v.size(); ++i) if (v[i] > m) m = v[i];
  return m;
}


NumericVector EnvBuildLinBound_cpp(NumericMatrix thetabars,
                                   NumericMatrix cbars,
                                   NumericVector y,
                                   NumericMatrix x,
                                   NumericMatrix P,
                                   NumericVector alpha,
                                   double dispstar) {
  // Convert to Armadillo
  arma::mat thetabarsA = as<arma::mat>(thetabars);
  arma::mat cbarsA     = as<arma::mat>(cbars);
  arma::vec yA         = as<arma::vec>(y);
  arma::mat xA         = as<arma::mat>(x);
  arma::mat PA         = as<arma::mat>(P);
  arma::vec alphaA     = as<arma::vec>(alpha);
  
  int gs = cbarsA.n_rows;
  
  arma::mat XtX   = xA.t() * xA;
  arma::vec rhs   = xA.t() * (yA - alphaA);
  arma::mat M     = XtX + dispstar * PA;
  arma::mat Minv  = arma::inv(M);           // match R's solve(M)
  arma::mat H1    = -Minv * PA * Minv;
  
  arma::mat V = -thetabarsA * PA + cbarsA;                 // gs x p
  arma::mat Minv_cbars = cbarsA * Minv.t();                // gs x p
  arma::vec term1 = arma::sum(V % Minv_cbars, 1);
  
  arma::mat rhs_mat = arma::repmat(rhs, 1, gs);            // p x gs
  arma::mat H1_rhs  = (H1 * (rhs_mat + dispstar * cbarsA.t())).t(); // gs x p
  arma::vec term2 = arma::sum(V % H1_rhs, 1);
  
  arma::vec result = term1 + term2;
  
  // Return explicitly as NumericVector
  NumericVector out(gs);
  std::copy(result.begin(), result.end(), out.begin());
  return out;
}


NumericVector thetabar_const_cpp(NumericMatrix P,
                                 NumericMatrix cbars,
                                 NumericMatrix thetabars) {
  arma::mat PA         = as<arma::mat>(P);
  arma::mat cbarsA     = as<arma::mat>(cbars);
  arma::mat thetabarsA = as<arma::mat>(thetabars);
  
  int gs = cbarsA.n_rows;
  arma::vec thetaconst(gs);
  
  for (int j = 0; j < gs; ++j) {
    arma::vec theta_temp = thetabarsA.row(j).t();
    arma::vec cbars_temp = cbarsA.row(j).t();
    thetaconst[j] = -0.5 * arma::as_scalar(theta_temp.t() * PA * theta_temp)
      + arma::as_scalar(cbars_temp.t() * theta_temp);
  }
  
  NumericVector out(gs);
  std::copy(thetaconst.begin(), thetaconst.end(), out.begin());
  return out;
}


// --- Internal helper: UB2 pilot timing block ---
// Not exported to R
Rcpp::List run_ub2_pilot_block(const Rcpp::Function& ub2_parallel_fn,
                               int gs, int l1,
                               double low, double upp,
                               const Rcpp::List& cache,
                               const Rcpp::NumericMatrix& cbars,
                               const Rcpp::NumericVector& y,
                               const Rcpp::NumericMatrix& x,
                               const Rcpp::NumericVector& alpha,
                               const Rcpp::NumericVector& wt,
                               double rss_min_global,
                               bool verbose) {
  double est_total = 0.0;
  // const int pilot_threshold = static_cast<int>(std::pow(3, 10)); // 59,049 faces
  
  // --- Warm-up pilot size ---
  int k1 = std::min(gs, 500);
  
  // Fractional pilots: ~0.5% and ~1.0% of total faces
  auto frac_round = [](double v) { return static_cast<int>(std::round(v)); };
  int k2_target = frac_round(0.005 * static_cast<double>(gs));
  int k3_target = frac_round(0.010 * static_cast<double>(gs));
  
  // Floors/caps
  int floor_k2 = 3000, floor_k3 = 6000;
  int cap_k2   = 50000, cap_k3 = 100000;
  
  int k2 = std::min(gs, std::max(floor_k2, std::min(k2_target, cap_k2)));
  int k3 = std::min(gs, std::max(floor_k3, std::min(k3_target, cap_k3)));
  if (k2 <= k1) k2 = std::min(gs, std::max(k1 + 1, floor_k2));
  if (k3 <= k2) k3 = std::min(gs, std::max(k2 + 1, floor_k3));
  
  auto make_slice = [&](int k) {
    Rcpp::NumericMatrix cbars_slice(k, l1);
    for (int i = 0; i < k; ++i)
      for (int j = 0; j < l1; ++j)
        cbars_slice(i, j) = cbars(i, j);
    return cbars_slice;
  };
  
  auto now_num = []() {
    return Rcpp::as<double>(
      Rcpp::Function("as.numeric")(Rcpp::Function("Sys.time")())
    );
  };
  
  // Pilot timings
  double t0 = now_num();
  ub2_parallel_fn(Rcpp::Named("par0")   = 0.5 * (low + upp),
                  Rcpp::Named("low")    = low,
                  Rcpp::Named("upp")    = upp,
                  Rcpp::Named("cache")  = cache,
                  Rcpp::Named("cbars")  = make_slice(k1),
                  Rcpp::Named("y")      = y,
                  Rcpp::Named("x")      = x,
                  Rcpp::Named("alpha")  = alpha,
                  Rcpp::Named("wt")     = wt,
                  Rcpp::Named("rss_min_global") = rss_min_global);
  double t1 = now_num();
  double elapsed1 = t1 - t0;
  
  double t2 = now_num();
  ub2_parallel_fn(Rcpp::Named("par0")   = 0.5 * (low + upp),
                  Rcpp::Named("low")    = low,
                  Rcpp::Named("upp")    = upp,
                  Rcpp::Named("cache")  = cache,
                  Rcpp::Named("cbars")  = make_slice(k2),
                  Rcpp::Named("y")      = y,
                  Rcpp::Named("x")      = x,
                  Rcpp::Named("alpha")  = alpha,
                  Rcpp::Named("wt")     = wt,
                  Rcpp::Named("rss_min_global") = rss_min_global);
  double t3 = now_num();
  double elapsed2 = t3 - t2;
  
  double t4 = now_num();
  ub2_parallel_fn(Rcpp::Named("par0")   = 0.5 * (low + upp),
                  Rcpp::Named("low")    = low,
                  Rcpp::Named("upp")    = upp,
                  Rcpp::Named("cache")  = cache,
                  Rcpp::Named("cbars")  = make_slice(k3),
                  Rcpp::Named("y")      = y,
                  Rcpp::Named("x")      = x,
                  Rcpp::Named("alpha")  = alpha,
                  Rcpp::Named("wt")     = wt,
                  Rcpp::Named("rss_min_global") = rss_min_global);
  double t5 = now_num();
  double elapsed3 = t5 - t4;
  
  // Estimate per-face slope
  double denom   = static_cast<double>(k3 - k2);
  double t_face  = (elapsed3 - elapsed2) / std::max(1.0, denom);
  double t_fixed = elapsed1;
  est_total      = t_fixed + static_cast<double>(gs) * t_face;
  
  auto fmt_hms = [](double seconds) {
    int s = static_cast<int>(std::round(seconds));
    int h = s / 3600; s %= 3600;
    int m = s / 60;   s %= 60;
    std::ostringstream oss;
    if (h) oss << h << "h ";
    if (h || m) oss << m << "m ";
    oss << s << "s";
    return oss.str();
  };
  
  Rcpp::Rcout << "[EnvelopeDispersionBuild:UB2:Pilot] k1=" << k1
              << " (" << (100.0 * k1 / gs) << "%) elapsed=" << elapsed1 << "s; "
              << "k2=" << k2 << " (" << (100.0 * k2 / gs) << "%) elapsed=" << elapsed2 << "s; "
              << "k3=" << k3 << " (" << (100.0 * k3 / gs) << "%) elapsed=" << elapsed3 << "s.\n";
  
  Rcpp::Rcout << "[EnvelopeDispersionBuild:UB2:Pilot] t_fixed=" << t_fixed
              << "s, t_face=" << t_face << "s/face.\n";
  
  Rcpp::Rcout << "[EnvelopeDispersionBuild:UB2:Pilot] Estimated full run = "
              << fmt_hms(est_total) << " (" << est_total << "s).\n";
  
  return Rcpp::List::create(Rcpp::Named("est_total") = est_total);
}



// ---------------------------------------------------------------------
// bound_rss_over_dispersion
// Closed-form lower bound for RSS_min per Appendix A07 (RSS ML Decomposition).
// RSS_j(d) = RSS_ML + quad_j(d); quad_j(d) >= C*m_j^2 with fixed a_j, b_j.
// LB = RSS_ML + C * m_j^2 (no d1_star).
// ---------------------------------------------------------------------
Rcpp::List bound_rss_over_dispersion(
    const Rcpp::List& cache,
    const Rcpp::List& Env,
    double RSS_ML,
    double shape2,
    double rate3,
    double low,
    double upp,
    bool verbose
) {
  using namespace Rcpp;
  using namespace arma;

  NumericMatrix cbars = Env["cbars"];
  int gs = cbars.nrow();
  int p  = cbars.ncol();

  arma::mat base_A  = cache["base_A"];
  arma::vec base_B0 = cache["base_B0"];
  arma::mat Pmat    = cache["Pmat"];
  arma::vec Pmu     = cache["Pmu"];
  Pmat = 0.5 * (Pmat + Pmat.t());

  bool dbg = verbose;
  if (dbg && !R_finite(RSS_ML)) Rcout << "[bound_rss:VALIDATE] RSS_ML non-finite\n";

  // Q = X' W X = base_A; beta_hat = ML
  arma::vec beta_hat = -arma::solve(base_A, base_B0);

  // M(low) and C = lambda_min(M(low)): use precomputed M_min from cache (Step 3B1.5 in
  // EnvelopeDispersionBuild) — same matrix as former internal block below, avoids
  // duplicate A_max / inv_sympd / M_min work.
  if (!cache.containsElementNamed("M_min")) {
    Rcpp::stop("bound_rss_over_dispersion: cache must contain M_min; call after Step 3B1.5 precompute.");
  }
  arma::mat M_min = Rcpp::as<arma::mat>(cache["M_min"]);
  M_min = 0.5 * (M_min + M_min.t());
  arma::vec evals_M = arma::eig_sym(M_min);
  double C = evals_M(0);

  // --- Internal duplicate of Step 3B1.5 (kept for reference; do not delete) ---
  // // A_max = P + base_A/low (needed for M-based C)
  // arma::mat A_max = Pmat + base_A / low;
  // A_max = 0.5 * (A_max + A_max.t());
  //
  // // Old C: C_old = lambda_min(Q) / lambda_max(A_max)^2 (commented out; M-based C used instead)
  // // arma::vec evals_Q = arma::eig_sym(base_A);
  // // double lambda_min_Q = evals_Q(0);
  // // if (dbg && (!R_finite(lambda_min_Q) || lambda_min_Q <= 0.0))
  // //   Rcout << "[bound_rss:VALIDATE] lambda_min(Q) non-finite or <=0\n";
  // // arma::vec evals_A_max = arma::eig_sym(A_max);
  // // double lambda_max_A = evals_A_max(evals_A_max.n_elem - 1);
  // // if (dbg && (!R_finite(lambda_max_A) || lambda_max_A <= 0.0))
  // //   Rcout << "[bound_rss:VALIDATE] lambda_max(A_max) non-finite or <=0\n";
  // // double C_old = lambda_min_Q / (lambda_max_A * lambda_max_A);
  //
  // // M-based bound: M(low) = A_max^{-1} Q A_max^{-1}, C = lambda_min(M(low)) bounds universally
  // arma::mat Ainv_max = arma::inv_sympd(A_max);
  // arma::mat M_min = Ainv_max.t() * base_A * Ainv_max;  // Q = base_A for Gaussian
  // M_min = 0.5 * (M_min + M_min.t());
  // arma::vec evals_M = arma::eig_sym(M_min);
  // double C = evals_M(0);

  // double t_min = 1.0 / upp;
  // double t_max = 1.0 / low;

  // Proposition 1 bound (two terms): RSS_j(d) >= RSS_ML + r*^T M(1/low) r*
  // M(t) decreasing in t => min at d=low. M(low) = A_max^{-1} Q A_max^{-1} = M_min
  // (A_max, M_min already computed above)

  NumericVector rss_bound_parallel(gs);
  NumericVector disp_min_parallel(gs);
  // NumericVector rss_bound_prop1_parallel(gs);  // old C*m^2 bound (kept for comparison)
  double rss_min_bound = R_PosInf;
  // double rss_min_bound_prop1 = R_PosInf;

  for (int j = 0; j < gs; ++j) {
    arma::vec cbar_j(p);
    for (int r = 0; r < p; ++r) cbar_j(r) = cbars(j, r);

    // r* = b_j for Prop1 (tangency at d -> low, t = 1/low)
    arma::vec b_j = cbar_j - Pmu - Pmat * beta_hat;

    // Old bound: a_j, t_star, t_tilde, r_min, m_j_sq, quad_lower_bound = C*m_j^2
    // arma::vec a_j = -(base_B0 + base_A * beta_hat);
    // double aa = arma::dot(a_j, a_j);
    // double ab = arma::dot(a_j, b_j);
    // double t_star = (aa > 0.0 ? -ab / aa : t_min);
    // double t_tilde = std::min(std::max(t_star, t_min), t_max);
    // arma::vec r_min = a_j * t_tilde + b_j;
    // double m_j_sq = arma::dot(r_min, r_min);
    // double quad_lower_bound = C * m_j_sq;
    // double lb_old = RSS_ML + quad_lower_bound;

    // Proposition 1 bound (two terms): RSS_ML + r*^T M(1/low) r* (output bound)
    double quad_prop1 = arma::as_scalar(b_j.t() * M_min * b_j);
    double lb_prop1_j = RSS_ML + quad_prop1;
    rss_bound_parallel(j) = lb_prop1_j;
    // rss_bound_prop1_parallel(j) = lb_old;
    disp_min_parallel(j) = low;                 // Prop1 bound achieved at d=low

    if (R_finite(lb_prop1_j) && lb_prop1_j < rss_min_bound) rss_min_bound = lb_prop1_j;
    // if (R_finite(lb_old) && lb_old < rss_min_bound_prop1) rss_min_bound_prop1 = lb_old;
  }

  // Optional: precompute UB2 curvature data for reuse (when gs <= 81)
  Rcpp::List out = List::create(
    Named("ok")                   = true,
    Named("rss_min_bound")        = rss_min_bound,
    Named("rss_bound_parallel")   = rss_bound_parallel,
    Named("disp_min_parallel")    = disp_min_parallel,
    Named("rss_min_global")       = rss_min_bound,
    Named("rss_min_parallel")     = rss_bound_parallel
  );
  if (gs <= 81) {
    // theta_ref = c^{-1}(cbar_j, d_ref): tangency at ref dispersion.
    // Residual r_j(d) = cbar_j - B0(d) - A(d)*theta_ref encodes d-dependence
    // of inverse; a_j, b_j give curvature f''(t) = 6A*t + 2B.
    double d_ref = 0.5 * (low + upp);
    double Delta0 = RSS_ML - rss_min_bound;
    NumericVector alpha_vec(gs), beta_vec(gs), gamma_vec(gs);
    NumericVector A_coef_vec(gs), B_coef_vec(gs), D_coef_vec(gs);
    double max_a_norm_sq = 0.0;
    double max_b0_plus_At = 0.0;
    for (int j = 0; j < gs; ++j) {
      arma::vec cbar_j(p);
      for (int r = 0; r < p; ++r) cbar_j(r) = cbars(j, r);
      NumericMatrix cbars_small(1, p);
      for (int r = 0; r < p; ++r) cbars_small(0, r) = cbar_j(r);
      arma::mat theta_mat = Inv_f3_with_disp(cache, d_ref, Rcpp::transpose(cbars_small));
      arma::vec theta_ref = theta_mat.row(0).t();
      arma::vec a_j = -(base_B0 + base_A * theta_ref);
      arma::vec b_j = cbar_j - Pmu - Pmat * theta_ref;
      double alpha_j = arma::dot(a_j, a_j);
      double beta_j  = 2.0 * arma::dot(a_j, b_j);
      double gamma_j = arma::dot(b_j, b_j);
      alpha_vec(j) = alpha_j;
      beta_vec(j)  = beta_j;
      gamma_vec(j) = gamma_j;
      A_coef_vec(j) = C * alpha_j / 2.0;
      B_coef_vec(j) = C * beta_j  / 2.0;
      D_coef_vec(j) = (Delta0 + C * gamma_j) / 2.0;
      if (alpha_j > max_a_norm_sq) max_a_norm_sq = alpha_j;
      double b0_plus_At_j = arma::norm(base_B0 + base_A * theta_ref, 2);
      if (b0_plus_At_j > max_b0_plus_At) max_b0_plus_At = b0_plus_At_j;
    }
    out["ub2_curv_ok"]       = true;
    out["ub2_curv_a_norm_sq"] = max_a_norm_sq;
    out["ub2_curv_b0_plus_Ab"] = max_b0_plus_At;
    out["ub2_curv_C"]        = C;
    out["ub2_curv_Delta0"]   = Delta0;
    out["ub2_curv_t_min"]    = 1.0 / upp;
    out["ub2_curv_t_max"]    = 1.0 / low;
    out["ub2_curv_alpha"]    = alpha_vec;
    out["ub2_curv_beta"]     = beta_vec;
    out["ub2_curv_gamma"]    = gamma_vec;
    out["ub2_curv_A_coef"]   = A_coef_vec;
    out["ub2_curv_B_coef"]   = B_coef_vec;
    out["ub2_curv_D_coef"]   = D_coef_vec;

  }
  return out;
}


// ---------------------------------------------------------------------
// rss_face_bound_from_cache_cpp — UNUSED on active path; entire definition
// commented out (do not delete). Only referenced from the large commented
// diagnostic block below (search rss_face_bound_from_cache_cpp). Uncomment
// this function if that block is re-enabled.
// ---------------------------------------------------------------------
/*
// Per-face diagnostic: compare closed-form RSS lower bound to actual
// minimized RSS. Uses same algebra as bound_rss_over_dispersion.
// Not exported; for internal/debug use only.
//
// Call sites: only inside a *commented-out* verbose diagnostic block in
// EnvelopeDispersionBuild (search for rss_face_bound_from_cache_cpp). It is not
// on the active code path unless that block is re-enabled.
//
// Using Q = base_A from cache (below) avoids recomputing X'WX from (X, wt); for
// Gaussian, Inv_f3_precompute_disp already stores Q as base_A — same matrix.
// ---------------------------------------------------------------------
Rcpp::List rss_face_bound_from_cache_cpp(
    Rcpp::List cache,
    Rcpp::NumericVector beta_hat_r,
    Rcpp::NumericVector cbars_j_r,
    double low,
    double upp,
    Rcpp::NumericVector y_r,
    Rcpp::NumericMatrix X_r,
    Rcpp::NumericVector alpha_r,
    Rcpp::NumericVector wt_r,
    double shape2,
    double rate3,
    double RSS_ML,
    double rss_min_j_actual,
    bool verbose = true,
    int face_j = -1
) {
  using namespace arma;

  arma::vec beta_hat(beta_hat_r.begin(), beta_hat_r.size(), false);
  arma::vec cbars_j(cbars_j_r.begin(), cbars_j_r.size(), false);
  arma::vec y(y_r.begin(), y_r.size(), false);
  arma::mat X(X_r.begin(), X_r.nrow(), X_r.ncol(), false);
  arma::vec alpha(alpha_r.begin(), alpha_r.size(), false);
  arma::vec wt(wt_r.begin(), wt_r.size(), false);

  arma::mat Pmat    = cache["Pmat"];
  arma::vec Pmu     = cache["Pmu"];
  arma::vec base_B0 = cache["base_B0"];
  arma::mat base_A  = cache["base_A"];
  Pmat = 0.5 * (Pmat + Pmat.t());

  // Q = X' W X — identical to cache's base_A from Inv_f3_precompute_disp (Gaussian).
  // Extract instead of recomputing from X, wt:
  //   arma::vec sqrtw = arma::sqrt(wt);
  //   arma::mat Xw = X.each_col() % sqrtw;
  //   arma::mat Q = Xw.t() * Xw;
  arma::mat Q = base_A;
  // Retain y, X, alpha, wt in signature for diagnostic callers; Q matches X'WX from cache.
  (void)y;
  (void)X;
  (void)alpha;
  (void)wt;

  arma::vec eig_Q = arma::eig_sym(Q);
  double lambda_min_Q = eig_Q.min();


  
  
  
  arma::mat A_max = Pmat + base_A / low;
  A_max = 0.5 * (A_max + A_max.t());
  arma::vec eig_Amax = arma::eig_sym(A_max);
  double lambda_max_Amax = eig_Amax(eig_Amax.n_elem - 1);
  double C = lambda_min_Q / (lambda_max_Amax * lambda_max_Amax);

  // M_min for eigenvalue diagnostic: prefer cache (same as bound_rss_over_dispersion) when present.
  arma::mat M_min_diag;
  if (cache.containsElementNamed("M_min")) {
    M_min_diag = Rcpp::as<arma::mat>(cache["M_min"]);
    M_min_diag = 0.5 * (M_min_diag + M_min_diag.t());
  } else {
    // Was: duplicate of Step 3B1.5
    // arma::mat A_max2 = Pmat + base_A / low;
    // A_max2 = 0.5 * (A_max2 + A_max2.t());
    // arma::mat Ainv_max = arma::inv_sympd(A_max2);
    // M_min_diag = Ainv_max.t() * base_A * Ainv_max;
    arma::mat A_max2 = Pmat + base_A / low;
    A_max2 = 0.5 * (A_max2 + A_max2.t());
    arma::mat Ainv_max = arma::inv_sympd(A_max2);
    M_min_diag = Ainv_max.t() * base_A * Ainv_max;
    M_min_diag = 0.5 * (M_min_diag + M_min_diag.t());
  }

  arma::vec evals_M_min = arma::eig_sym(M_min_diag);
  Rcout << "[rss_face_bound:eigenvalues of M] ";
  for (arma::uword i = 0; i < evals_M_min.n_elem; ++i)
    Rcout << evals_M_min(i) << (i + 1 < evals_M_min.n_elem ? " " : "\n");
  double C2 = evals_M_min(0);  // lambda_min(M_min)

  Rcout << "[Curvature bounds]  C_old=" << C 
        << "   C_new=" << C2 
        << "   (ratio=" << C2 / C << ")\n";
    
  // Fixed a_j, b_j: r_j(t) = a_j*t + b_j with t = 1/d
  arma::vec a_j = -(base_B0 + base_A * beta_hat);
  arma::vec b_j = cbars_j - Pmu - Pmat * beta_hat;

  double t_min = 1.0 / upp;
  double t_max = 1.0 / low;

  double aa = arma::dot(a_j, a_j);
  double ab = arma::dot(a_j, b_j);
  double t_star = (aa > 0.0 ? -ab / aa : t_min);
  double t_tilde = std::min(std::max(t_star, t_min), t_max);

  arma::vec r_min = a_j * t_tilde + b_j;
  double m_j2 = arma::dot(r_min, r_min);
  double quad_lower_bound = C * m_j2;
  double RSS_lb_j = RSS_ML + quad_lower_bound;

  if (verbose) {
    double diff_j = rss_min_j_actual - RSS_lb_j;
    if (face_j >= 0) {
      Rcpp::Rcout << "  " << face_j << " | bound=" << RSS_lb_j
                  << " actual=" << rss_min_j_actual
                  << " diff=" << diff_j << "\n";
    } else {
      Rcpp::Rcout << "Closed-form RSS lower bound: " << RSS_lb_j << "\n";
      Rcpp::Rcout << "Actual minimized RSS:        " << rss_min_j_actual << "\n";
      Rcpp::Rcout << "Difference (actual - bound): " << diff_j << "\n";
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("RSS_lower_bound") = RSS_lb_j,
    Rcpp::Named("RSS_actual")      = rss_min_j_actual,
    Rcpp::Named("C")               = C,
    Rcpp::Named("m_j2")            = m_j2,
    Rcpp::Named("t_tilde")         = t_tilde,
    Rcpp::Named("a_j")             = a_j,
    Rcpp::Named("b_j")             = b_j,
    Rcpp::Named("Q")               = Q
  );
}
*/


// ---------------------------------------------------------------------
// rss_face_quadratic_sum_internal
// Computes closed-form lower bound on quad: quad >= quad_lower_bound for all d in [low, upp].
// No dispersion, no printing, no RSS-minimization dependency. Q = base_A from cache.
// Caller does quad-at-d diagnostics externally (e.g. via rss_face_at_disp - RSS_ML).
// NOTE: no [[Rcpp::export]] — internal only
// ---------------------------------------------------------------------
Rcpp::List rss_face_quadratic_sum_internal(
    Rcpp::List cache,
    Rcpp::NumericVector cbars_j_r,
    Rcpp::NumericVector beta_hat_r,
    double low,
    double upp
) {
  using namespace arma;

  vec cbars_j(cbars_j_r.begin(), cbars_j_r.size(), false);
  vec beta_hat(beta_hat_r.begin(), beta_hat_r.size(), false);

  mat Pmat    = cache["Pmat"];
  vec Pmu     = cache["Pmu"];
  vec base_B0 = cache["base_B0"];
  mat base_A  = cache["base_A"];
  Pmat = 0.5 * (Pmat + Pmat.t());
  mat Q = base_A;

  vec a_j = -(base_B0 + base_A * beta_hat);
  vec b_j = cbars_j - Pmu - Pmat * beta_hat;

  double t_min = 1.0 / upp;
  double t_max = 1.0 / low;
  double aa = dot(a_j, a_j);
  double ab = dot(a_j, b_j);
  double t_star = (aa > 0.0 ? -ab / aa : t_min);
  double t_tilde = std::min(std::max(t_star, t_min), t_max);

  vec r_min = a_j * t_tilde + b_j;
  double m_j2 = dot(r_min, r_min);

  mat A_max = Pmat + base_A / low;
  A_max = 0.5 * (A_max + A_max.t());
  vec eig_Amax = eig_sym(A_max);
  double lambda_max_Amax = eig_Amax(eig_Amax.n_elem - 1);
  vec eig_Q = eig_sym(Q);
  double lambda_min_Q = eig_Q.min();
  double C = lambda_min_Q / (lambda_max_Amax * lambda_max_Amax);
  double quad_lower_bound = C * m_j2;

  return Rcpp::List::create(
    Rcpp::Named("quad_lower_bound") = quad_lower_bound,
    Rcpp::Named("C")                = C,
    Rcpp::Named("m_j2")             = m_j2,
    Rcpp::Named("t_tilde")          = t_tilde,
    Rcpp::Named("a_j")              = a_j,
    Rcpp::Named("b_j")              = b_j,
    Rcpp::Named("beta_hat")         = beta_hat,
    Rcpp::Named("Q")                = Q
  );
}



// ---------------------------------------------------------------------
// Exact (root-finding) minimization of UB2_j(d) over d in [low, upp] for
// anisotropic coefficient priors.
//
// Background (see vignettes/Chapter-A07.Rmd, Remark 5.5.4/5.5.7, and
// glmbayesCore's data-raw/ub2_root_finding_prototype.R /
// data-raw/README_ub2_rootfinding_fix.md): with t = 1/d,
// K = Q^{-1/2} P Q^{-1/2}, v_j = Q^{-1/2}*(cbars_j - P*mu - P*beta_hat), and
// w_i = (u_i^T v_j)^2 (u_i = eigenvectors of K, lambda_i its eigenvalues),
//
//   tilde{UB2}_j(t) = (t/2) * (g(t) - Delta),   g(t) = sum_i w_i/(lambda_i+t)^2.
//
// Claim 7 (Chapter A07) assumed the minimum over t always occurs at an
// endpoint, but the underlying proof (Remark 5.5.7) only guarantees any
// critical point t* satisfies t* < lambda_max(K), not t* < lambda_min(K).
// For anisotropic K this allows genuine interior minima, which the
// endpoint-only shortcut misses -- the mechanism behind observed
// "Sign violation: UB2 < 0" errors downstream (the endpoint estimate is too
// large, so an actually-evaluated dispersion near the true interior minimum
// undercuts it). This block finds the true minimum exactly, at negligible
// extra cost: t is always scalar regardless of the coefficient dimension p,
// and K depends only on P/Q (not on the face j), so its eigendecomposition
// is computed once per envelope build and reused across all faces.
// ---------------------------------------------------------------------

namespace ub2_exact_detail {

inline double g_of_t(const arma::vec& lambda, const arma::vec& w, double t) {
  double s = 0.0;
  for (arma::uword i = 0; i < lambda.n_elem; ++i) {
    double m = lambda(i) + t;
    s += w(i) / (m * m);
  }
  return s;
}

inline double hprime_of_t(const arma::vec& lambda, const arma::vec& w, double t) {
  double s = 0.0;
  for (arma::uword i = 0; i < lambda.n_elem; ++i) {
    double m = lambda(i) + t;
    s += w(i) * (lambda(i) - t) / (m * m * m);
  }
  return s;
}

inline double ub2_reduced(const arma::vec& lambda, const arma::vec& w, double Delta, double t) {
  return 0.5 * t * (g_of_t(lambda, w, t) - Delta);
}

// Robust bracketed bisection root finder for f(t) = 0 on [a, b], f(a) and
// f(b) assumed to have opposite signs (or zero). Not the fastest possible
// (Brent's method would converge faster) but simple and reliably correct;
// this runs a handful of times per face, once per envelope build, so
// performance is not a concern.
template <typename F>
inline double bisection_root(F f, double a, double b, double fa, double fb,
                              int max_iter = 100, double tol = 1e-12) {
  for (int it = 0; it < max_iter; ++it) {
    double mid = 0.5 * (a + b);
    double fm = f(mid);
    if (std::abs(fm) < tol || (b - a) < tol * std::max(1.0, std::abs(mid))) return mid;
    if ((fa > 0.0) == (fm > 0.0)) { a = mid; fa = fm; } else { b = mid; fb = fm; }
  }
  return 0.5 * (a + b);
}

struct ExactResult { double ub2_min; double t_star; };

// Finds the exact minimum of tilde{UB2}_j(t) over t in [t_lo, t_hi] by
// bracketing sign changes of h'(t) - Delta on a grid anchored at t_lo, t_hi,
// and min(t_hi, lambda_max(K)) (any critical point must lie strictly below
// lambda_max(K); see Remark 5.5.7), refined near each eigenvalue of K, then
// polishing each bracket via bisection. Always includes t_lo and t_hi among
// the evaluated candidates, so the result is never worse than the old
// endpoint-only estimate.
inline ExactResult ub2_min_exact_1d(const arma::vec& lambda, const arma::vec& w,
                                     double Delta, double t_lo, double t_hi,
                                     int grid_mult = 40) {
  double lam_max = lambda.max();
  double hi_search = std::min(t_hi, lam_max * (1.0 - 1e-9));

  std::vector<double> cands;
  cands.push_back(t_lo);
  cands.push_back(t_hi);

  if (hi_search > t_lo) {
    std::vector<double> anchors;
    anchors.push_back(t_lo);
    anchors.push_back(hi_search);
    for (arma::uword i = 0; i < lambda.n_elem; ++i) {
      if (lambda(i) > t_lo && lambda(i) < hi_search) anchors.push_back(lambda(i));
    }
    std::sort(anchors.begin(), anchors.end());
    anchors.erase(std::unique(anchors.begin(), anchors.end()), anchors.end());

    std::vector<double> grid;
    for (size_t i = 0; i + 1 < anchors.size(); ++i) {
      double lo_i = anchors[i], hi_i = anchors[i + 1];
      if (hi_i <= lo_i) continue;
      for (int k = 0; k < grid_mult; ++k) {
        grid.push_back(lo_i + (hi_i - lo_i) * static_cast<double>(k) / (grid_mult - 1));
      }
    }
    if (grid.size() < 2) {
      for (int k = 0; k < grid_mult; ++k) {
        grid.push_back(t_lo + (hi_search - t_lo) * static_cast<double>(k) / (grid_mult - 1));
      }
    }
    std::sort(grid.begin(), grid.end());
    grid.erase(std::unique(grid.begin(), grid.end()), grid.end());

    std::vector<double> fvals(grid.size());
    for (size_t i = 0; i < grid.size(); ++i) fvals[i] = hprime_of_t(lambda, w, grid[i]) - Delta;

    for (size_t i = 0; i + 1 < grid.size(); ++i) {
      if ((fvals[i] > 0.0) != (fvals[i + 1] > 0.0)) {
        double root = bisection_root(
          [&](double t) { return hprime_of_t(lambda, w, t) - Delta; },
          grid[i], grid[i + 1], fvals[i], fvals[i + 1]
        );
        if (root > t_lo && root < t_hi) cands.push_back(root);
      }
    }
  }

  double best_val = ub2_reduced(lambda, w, Delta, cands[0]);
  double best_t = cands[0];
  for (size_t i = 1; i < cands.size(); ++i) {
    double v = ub2_reduced(lambda, w, Delta, cands[i]);
    if (v < best_val) { best_val = v; best_t = cands[i]; }
  }
  return ExactResult{ best_val, best_t };
}

}  // namespace ub2_exact_detail

// ---------------------------------------------------------------------
// bound_ub2_over_dispersion
// Evaluate UB2 at the dispersion endpoints [low, upp] for each face, plus
// (when the coefficient prior precision P/Q pair yields a well-defined K)
// the exact interior minimum via ub2_exact_detail::ub2_min_exact_1d above;
// return the per-face minimum and the matching dispersion.
// ---------------------------------------------------------------------
Rcpp::List bound_ub2_over_dispersion(
    int gs,
    double low,
    double upp,
    const Rcpp::List& cache,
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    double rss_min_global
    ) {
  using namespace Rcpp;
  using namespace arma;

  int p = static_cast<int>(cbars.ncol());

  NumericVector disp_min_ub2(gs);
  NumericVector ub2_min(gs);

  // Local selector for UB2 endpoint method:
  // 1 = UB2() via rss_face_at_disp (current default)
  // 2 = RSS-based quadratic form using M_min / M_max at low/upp
  int UB_Min_Method = 2;

  // Shared matrices (extracted once; only used when UB_Min_Method == 2)
  mat base_A  = cache["base_A"];
  vec base_B0 = cache["base_B0"];
  mat Pmat    = cache["Pmat"];
  vec Pmu     = cache["Pmu"];
  Pmat = 0.5 * (Pmat + Pmat.t());

  vec beta_hat;
  mat M_min, M_max;
  double RSS_ML_local = 0.0;

  // Exact (root-finding) minimization needs K = Q^{-1/2} P Q^{-1/2} and its
  // eigendecomposition. K depends only on P and Q = base_A, not on the face
  // j, so both are computed once here and reused for every face below.
  bool have_K = false;
  vec K_eigval;
  mat K_eigvec;
  mat Qinvhalf;   // Rq^{-1}, where base_A = Rq^T * Rq (arma::chol, upper-tri)
  double Delta = 0.0;
  double t_lo = 1.0 / upp;   // t = 1/d; d=upp -> t_lo (smallest t)
  double t_hi = 1.0 / low;   // d=low -> t_hi (largest t)

  // Near-isotropic fast path: Claim 7 part 3 (Chapter A07 vignette) proves
  // the endpoint-only minimum is *exact* -- not merely a heuristic -- once
  // kappa(K) = lambda_max(K)/lambda_min(K) <= 2 (any critical point lies in
  // t* < lambda_max(K) by part 1, while any inflection point needs
  // t* >= 2*lambda_min(K) by part 2; these two ranges cannot overlap when
  // lambda_max(K) <= 2*lambda_min(K), so no genuine interior local minimum
  // is possible). When this certificate holds we skip the per-face
  // root-finding search entirely (it would just re-derive the endpoint
  // values at extra cost) and fall straight through to the endpoint
  // comparison below. This is a single, cheap, once-per-envelope-build
  // check -- not a per-face search -- since kappa(K) does not depend on j.
  bool K_is_near_isotropic = false;
  double kappa_K = NA_REAL;

  if (UB_Min_Method == 2) {
    M_min = Rcpp::as<mat>(cache["M_min"]);
    M_max = Rcpp::as<mat>(cache["M_max"]);

    beta_hat      = -solve(base_A, base_B0);
    mat X         = Rcpp::as<mat>(x);
    vec yv        = Rcpp::as<vec>(y);
    vec alphav    = Rcpp::as<vec>(alpha);
    vec wv        = Rcpp::as<vec>(wt);
    vec resid_ml  = yv - X * beta_hat - alphav;
    RSS_ML_local  = as_scalar(resid_ml.t() * (wv % resid_ml));

    Delta = rss_min_global - RSS_ML_local;
    if (Delta < 0.0) Delta = 0.0;  // guard against tiny floating-point noise

    mat Rq;
    if (arma::chol(Rq, base_A)) {
      Qinvhalf = arma::inv(arma::trimatu(Rq));
      mat K = Qinvhalf.t() * Pmat * Qinvhalf;
      K = 0.5 * (K + K.t());
      if (arma::eig_sym(K_eigval, K_eigvec, K) && K_eigval.min() > 0.0) {
        have_K = true;
        kappa_K = K_eigval.max() / K_eigval.min();
        // Small numerical cushion above the exact "2" cutoff from Claim 7
        // part 3, since kappa_K is itself only known to eigendecomposition
        // precision; using the exact cutoff risks needlessly falling back
        // to the (harmless, still-correct) search for borderline cases.
        K_is_near_isotropic = kappa_K <= 2.0 * (1.0 + 1e-8);
      }
    }
  }

  for (int j = 0; j < gs; ++j) {
    NumericVector cbars_j(p);
    for (int r = 0; r < p; ++r) cbars_j[r] = cbars(j, r);

    double ub2_low = NA_REAL;
    double ub2_upp = NA_REAL;
    double best_ub2, best_disp;
    bool used_exact = false;

    if (UB_Min_Method == 1) {
      // Method 1: original UB2 helper (endpoint-only; retained for
      // compatibility if UB_Min_Method is manually switched back to 1).
      ub2_low = UB2(low, cache, cbars_j, y, x, alpha, wt, rss_min_global);
      ub2_upp = UB2(upp, cache, cbars_j, y, x, alpha, wt, rss_min_global);
    } else {
      // Method 2: RSS-based quadratic form with M_min / M_max
      vec cbar_j_vec(p);
      for (int r = 0; r < p; ++r) cbar_j_vec(r) = cbars(j, r);
      vec b_j = cbar_j_vec - Pmu - Pmat * beta_hat;

      double rss_low_approx = RSS_ML_local + as_scalar(b_j.t() * M_min * b_j);
      double rss_upp_approx = RSS_ML_local + as_scalar(b_j.t() * M_max * b_j);

      ub2_low = (0.5 / low) * (rss_low_approx - rss_min_global);
      ub2_upp = (0.5 / upp) * (rss_upp_approx - rss_min_global);

      if (have_K && !K_is_near_isotropic) {
        // Exact minimum over the whole [low, upp] interval, not just the
        // endpoints -- fixes the Claim 7 gap for anisotropic K. Always
        // evaluates at t_lo/t_hi too, so this is never worse than the
        // endpoint-only estimate above, and reduces to it automatically
        // when K is (numerically) isotropic.
        vec v_j = Qinvhalf.t() * b_j;
        vec w_coords = arma::square(K_eigvec.t() * v_j);
        ub2_exact_detail::ExactResult ex =
          ub2_exact_detail::ub2_min_exact_1d(K_eigval, w_coords, Delta, t_lo, t_hi);
        best_ub2   = ex.ub2_min;
        best_disp  = 1.0 / ex.t_star;
        used_exact = true;
      }
      // else: K_is_near_isotropic (kappa(K) <= 2) certifies -- once, for the
      // whole envelope build, per the note above bound_ub2_over_dispersion's
      // K/kappa_K computation -- that the endpoint-only comparison below is
      // already exact, so the per-face root-finding search is skipped.
    }

    if (!used_exact) {
      if (ub2_low <= ub2_upp) {
        best_ub2  = ub2_low;
        best_disp = low;
      } else {
        best_ub2  = ub2_upp;
        best_disp = upp;
      }
    }

    disp_min_ub2[j] = best_disp;
    ub2_min[j]      = best_ub2;
  }

  return Rcpp::List::create(
    Rcpp::Named("disp_min_ub2")        = disp_min_ub2,
    Rcpp::Named("ub2_min")             = ub2_min,
    Rcpp::Named("kappa_K")             = kappa_K,
    Rcpp::Named("K_is_near_isotropic") = K_is_near_isotropic
  );
}


// ---------------------------------------------------------------------
// Internal helper: minimize UB2 over dispersion for all faces
// Not exported. Only visible inside this .cpp file.
// ---------------------------------------------------------------------
Rcpp::List minimize_ub2_over_dispersion(
    int gs,
    int l1,
    double low,
    double upp,
    const Rcpp::List& cache,
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    double rss_min_global,
    const Rcpp::NumericVector& rss_min_parallel,
    double RSS_ML,
    int RSS_Min_Type,
    int UB2_Min_Type,
    bool verbose,
    const Rcpp::List& rss_bound_res,
    const Rcpp::NumericVector& disp_start
) {
  using namespace Rcpp;
  using namespace arma;

  NumericVector disp_min_ub2(gs);
  NumericVector ub2_min(gs);
  
  double ub2_min_global      = R_PosInf;
  double disp_min_global_ub2 = NA_REAL;
  int    j_best_ub2          = -1;
  
  // -------------------------------------------------------------------
  // Case 1: UB2 minimization is performed (UB2_Min_Type == 1)
  // -------------------------------------------------------------------
  if (UB2_Min_Type == 1) {
    // UB2 parallel minimization via the former R callback has been removed.
    // We keep `UB2_Min_Type == 1` as a compatibility mode, but compute UB2
    // minima using endpoint evaluation only.
    if (RSS_Min_Type == 1) {
      for (int j = 0; j < gs; ++j) {
        NumericVector cbars_j(cbars.ncol());
        for (int k = 0; k < cbars.ncol(); ++k) cbars_j[k] = cbars(j, k);
        
        double rss_at_low = rss_face_at_disp(low, cache, cbars_j, y, x, alpha, wt);
        double rss_at_upp = rss_face_at_disp(upp, cache, cbars_j, y, x, alpha, wt);
        
        double ub2_at_low = (0.5 / low) * (rss_at_low - rss_min_global);
        double ub2_at_upp = (0.5 / upp) * (rss_at_upp - rss_min_global);
        
        if (ub2_at_low <= ub2_at_upp) {
          ub2_min[j]      = ub2_at_low;
          disp_min_ub2[j] = low;
        } else {
          ub2_min[j]      = ub2_at_upp;
          disp_min_ub2[j] = upp;
        }
      }
      if (verbose) {
        Rcout << "[EnvelopeDispersionBuild] UB2 source = endpoint eval via rss_face_at_disp (no R callback)\n";
      }
    } else if (RSS_Min_Type == 2) {
      for (int j = 0; j < gs; ++j) {
        ub2_min[j]      = 0.0;
        disp_min_ub2[j] = upp;  // enforce upper bound anchor
      }
      if (verbose) {
        Rcout << "[EnvelopeDispersionBuild] UB2 source = Set to 0 (skip RSS_Min and UB2 Min)\n";
      }
    }

    // Find global UB2 minimum
    ub2_min_global      = R_PosInf;
    disp_min_global_ub2 = NA_REAL;
    j_best_ub2          = -1;
    for (int j = 0; j < gs; ++j) {
      if (ub2_min[j] < ub2_min_global) {
        ub2_min_global      = ub2_min[j];
        disp_min_global_ub2 = disp_min_ub2[j];
        j_best_ub2          = j;
      }
    }

    return List::create(
      Named("ub2_min")           = ub2_min,
      Named("disp_min_ub2")      = disp_min_ub2,
      Named("ub2_min_global")    = ub2_min_global,
      Named("disp_min_global_ub2") = disp_min_global_ub2,
      Named("j_best_ub2")        = j_best_ub2
    );
  } else { // UB2_Min_Type == 2
    
    if (RSS_Min_Type == 1) {
      // Evaluate UB2 at low and upp using rss_face_at_disp (same as UB2 minimizer)
      for (int j = 0; j < gs; ++j) {
        NumericVector cbars_j(cbars.ncol());
        for (int k = 0; k < cbars.ncol(); ++k) cbars_j[k] = cbars(j, k);
        double rss_at_low = rss_face_at_disp(low, cache, cbars_j, y, x, alpha, wt);
        double rss_at_upp = rss_face_at_disp(upp, cache, cbars_j, y, x, alpha, wt);
        double ub2_at_low = (0.5 / low) * (rss_at_low - rss_min_global);
        double ub2_at_upp = (0.5 / upp) * (rss_at_upp - rss_min_global);
        if (ub2_at_low <= ub2_at_upp) {
          ub2_min[j]      = ub2_at_low;
          disp_min_ub2[j] = low;
        } else {
          ub2_min[j]      = ub2_at_upp;
          disp_min_ub2[j] = upp;
        }
      }
      if (verbose) {
        Rcout << "[EnvelopeDispersionBuild] UB2 source = endpoint eval via rss_face_at_disp (skip UB2)\n";
      }
      
    } else if (RSS_Min_Type == 2) {
      // RSS not minimized, UB2 skipped: set UB2 to 0
      for (int j = 0; j < gs; ++j) {
        ub2_min[j]      = 0.0;
        disp_min_ub2[j] = upp;  // enforce upper bound anchor
      }
      if (verbose) {
        Rcout << "[EnvelopeDispersionBuild] UB2 source = Set to 0 (skip RSS_Min and UB2 Min)\n";
      }
    }
  }
  
  return List::create(
    Named("ub2_min")           = ub2_min,
    Named("disp_min_ub2")      = disp_min_ub2,
    Named("ub2_min_global")    = ub2_min_global,
    Named("disp_min_global_ub2") = disp_min_global_ub2,
    Named("j_best_ub2")        = j_best_ub2
  );
}


// ---------------------------------------------------------------------
// Internal helper: Envelope geometry construction
// Pure geometry: no mixture weights, no UB2, no packaging.
// ---------------------------------------------------------------------
Rcpp::List compute_envelope_geometry_cpp(
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericMatrix& thetabars,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& P,      // FIXED
    const Rcpp::NumericVector& alpha,
    double low,
    double upp,
    double shape2,
    double rate3,
    double n_w
) {
  using namespace Rcpp;
  
   int gs = cbars.nrow();
  
  // Step 4: Base face constants
  NumericVector thetabar_const_base =thetabar_const_cpp(P, cbars, thetabars);
  
  // Step 5: initial anchor (posterior mean)

  double d1_star = rate3 / (shape2 - 1.0);   // posterior-mode anchor
  
  // Rcpp::Rcout << "compute_envelope_geometry_cpp: dispstar = " << dispstar << "\n";
  
  // Step 6: Face slopes at dispstar
  NumericVector New_LL_Slope = EnvBuildLinBound_cpp(thetabars, cbars, y, x, P, alpha, d1_star);
  

  // Step 7: Linear extrapolation to bounds
  NumericVector thetabar_const_upp_apprx(gs), thetabar_const_low_apprx(gs);
  for (int j = 0; j < gs; ++j) {
    thetabar_const_upp_apprx[j] =
      thetabar_const_base[j] + (upp - d1_star) * New_LL_Slope[j];
    thetabar_const_low_apprx[j] =
      thetabar_const_base[j] + (low - d1_star) * New_LL_Slope[j];
  }
  
  // Step 8: Global upper line geometry
  double max_low = max_vec(thetabar_const_low_apprx);
  double max_upp = max_vec(thetabar_const_upp_apprx);
  
  // We later compute
  // 1) double pf_upp = thetabar_const_upp_apprx[j] - max_upp;
  // 2) double pf_low = thetabar_const_low_apprx[j] - max_low;
  // 3) prob_factor[j]  = (pf_upp > pf_low ? pf_upp : pf_low);
  
  //  And add the result to UB3A--->We subtract the most allowed to keep the global line above the face specific line
  
  //////////////////////////////////////////
  // ---------------------------------------------------------------------------
  // UB3A/UB3B COMPATIBLE SLOPE GEOMETRY
  //
  // UB3A introduces a global linear upper bound in d:
  //
  //      UB3A_line(d) = lmc1 + lmc2 * d
  //
  // where lmc2 is constructed as an "average" (or weighted average) of the
  // face-specific slopes New_LL_Slope[j].  Once lmc2 is chosen, lmc1 is the
  // smallest intercept such that UB3A_line(d) dominates all face lines on
  // [low, upp].  Thus lmc1 depends on lmc2, but the slope is the primary driver.
  //
  // UB3B must then wrap this same UB3A line inside a log-linear bounding term:
  //
  //      UB3B_bound(d) = C + lm_log1 + lm_log2 * log(d)
  //
  // and UB3B(d) = UB3B_bound(d) - UB3A_line(d) must be >= 0 on [low, upp].
  //
  // To ensure consistency between UB3A and UB3B, we force the two lines to match
  // at two anchor points d_a and d_b.  Subtracting the two matching equations
  // gives the exact slope relation:
  //
  //      lm_log2 = lmc2 * (d_b - d_a) / (log(d_b) - log(d_a))
  //
  // This is the "compatible UB3A slope" formula: the UB3B log-tilt coefficient
  // lm_log2 is determined entirely by the UB3A slope lmc2 and the chosen anchors.
  // The intercept lmc1 cancels out and plays no role in this slope relation.
  //
  // The Gamma proposal uses:
  //
  //      shape3 = shape2 - lm_log2
  //
  // so lm_log2 must satisfy lm_log2 <= shape2 to keep shape3 > 0.  We impose the
  // stronger condition lm_log2 <= n_w / 2 (the data contribution to shape2).
  //
  // Therefore, before constructing UB3A, we cap the global slope so that the
  // implied lm_log2 remains feasible:
  //
  //      new_slope <= (n_w / 2) * (log(d_b) - log(d_a)) / (d_b - d_a)
  //
  // In practice, we compute the mean (or weighted mean) of New_LL_Slope[j], then
  // replace it by the smaller of:
  //      (i) that mean, and
  //      (ii) the maximum slope compatible with lm_log2 <= n_w/2.
  //
  // This ensures UB3A and UB3B remain algebraically consistent and guarantees
  // shape3 > 0 for the Gamma proposal.
  // ---------------------------------------------------------------------------
  
  double lmc2_max = (n_w / 2.0) * (std::log(upp) - std::log(low)) / (upp - low);
  
  // Mean-slope correction (parity with original)
//  double m_New_LL_Slope = Rcpp::mean(New_LL_Slope);
  

  double mean_slope = static_cast<double>(Rcpp::mean(New_LL_Slope));

  // Issue a warning if the UB3A mean slope exceeds the UB3B-compatible maximum
  if (mean_slope > lmc2_max) {
    Rcpp::Rcout
    << "[WARNING] UB3A mean slope (" << mean_slope
    << ") exceeds UB3B-compatible maximum (" << lmc2_max << ").\n"
    << "          Capping global slope to preserve lm_log2 <= n_w/2 "
    << "and ensure shape3 > 0.\n";
  }
  
  // Update the formular to cap m_New_LL_slope
  double m_New_LL_Slope = std::min(mean_slope, lmc2_max);
  
  // Compute the three quantities
  double max_low_mean      = max_upp - m_New_LL_Slope * (upp - low);
  double max_low_mean_old  = max_upp - mean_slope      * (upp - low);
  
  max_low = max_low_mean;

  double new_slope = (max_upp - max_low) / (upp - low);
  double new_int   = max_low - new_slope * low;

  double d2_star = (upp - low) / (std::log(upp / low));   // log-tilt anchor
  
    
  // Return all geometry objects
  return List::create(
    Named("thetabar_const_low_apprx") = thetabar_const_low_apprx,
    Named("thetabar_const_upp_apprx") = thetabar_const_upp_apprx,
    Named("max_low")                  = max_low_mean_old,
    Named("max_upp")                  = max_upp,
    Named("new_slope")                = new_slope,
    Named("new_int")                  = new_int,
    Named("d2_star")                 = d2_star
  );
}



// ---------------------------------------------------------------------
// Internal helper: mixture weights, gamma tilt, UB_list, diagnostics
// Updates the existing Env by adding/overwriting PLSD.
// ---------------------------------------------------------------------



Rcpp::List compute_mixture_and_outputs_cpp(
    Rcpp::List Env,   // existing envelope (must contain "cbars")
    

    // Existing extrapolated constants
    const Rcpp::NumericVector& thetabar_const_low_apprx,
    const Rcpp::NumericVector& thetabar_const_upp_apprx,
    

    // UB2 minima
    const Rcpp::NumericVector& ub2_min,
    
    // Mixture log-weights
    const Rcpp::NumericVector& logP1,
    
    // Global UB3A geometry (still used for legacy UB3A)
    double max_low,
    double max_upp,
    double new_slope,
    double new_int,
    
    double d2_star,
    
    // Gamma proposal parameters
    double shape2,
    double Rate,
    
    // Dispersion bounds
    double low,
    double upp,
    
    // RSS bounds
    double RSS_ML,
    double rss_min_global,
    
    bool verbose
) {
  int gs = logP1.size();
  
  // cbars from Env (needed for 0.5 * ||c_j||^2 term)
  NumericMatrix cbars = Env["cbars"];
  int l1 = cbars.ncol();
  
  NumericVector New_logP2(gs);
  NumericVector prob_factor(gs);
  NumericVector prob_factor2(gs);

  // --- Step 9: Mixture weights per face (match original) ---
  for (int j = 0; j < gs; ++j) {
    Rcpp::checkUserInterrupt();
    
    double pf_upp = thetabar_const_upp_apprx[j] - max_upp;
    double pf_low = thetabar_const_low_apprx[j] - max_low;
    
    prob_factor[j]  = (pf_upp > pf_low ? pf_upp : pf_low);
    prob_factor2[j] = prob_factor[j] - ub2_min[j];

    
    
    /////////////////////////////////////////////////////////
        
    double norm2 = 0.0;
    for (int k = 0; k < l1; ++k) {
      double cjk = cbars(j, k);
      norm2 += cjk * cjk;
    }
    New_logP2[j] = logP1[j] + 0.5 * norm2;
  }
  
  
  NumericVector lg_prob_factor  = clone(prob_factor);

  // --- Stable PLSD computation (prob_factor_exp2 only) ---
  NumericVector prob_factor_exp2(gs);
  
  NumericVector logw2(gs);
  for (int j = 0; j < gs; ++j) {
    logw2[j] = New_logP2[j] + prob_factor2[j];
  }
  
  double max_logw2 = Rcpp::max(logw2);
  double sumP2 = 0.0;
  for (int j = 0; j < gs; ++j) {
    prob_factor_exp2[j] = std::exp(logw2[j] - max_logw2);
    sumP2 += prob_factor_exp2[j];
  }
  
  if (sumP2 <= 0.0 || !R_finite(sumP2)) {
    Rcpp::stop("PLSD normalization failed: sumP2 non-finite or non-positive");
  }

  
  for (int j = 0; j < gs; ++j) {
    prob_factor_exp2[j] /= sumP2;
      }
  
  Env["PLSD"] = prob_factor_exp2;
  
  double lm_log2 = new_slope * d2_star;
  double lm_log1 = new_int + new_slope * d2_star - new_slope * std::log(d2_star);
  double shape3  = shape2 - lm_log2;
  
  

  ////////////////////////////////////////////////////////////////////////
  
  // --- Sanity checks on the tilted gamma parameters (global) ---
  if (!R_finite(lm_log1) || !R_finite(lm_log2)) {
    Rcpp::stop("EnvelopeDispersionBuild: lm_log1/lm_log2 non-finite; envelope tilt is undefined.");
  }
   if (!R_finite(shape3) || shape3 <= 0.0) {
     Rcpp::stop("EnvelopeDispersionBuild: implied shape3 <= 0; tilted inverse-gamma is invalid.");
   }
  
  double rate2 = Rate + rss_min_global / 2.0;
  if (!R_finite(rate2) || rate2 <= 0.0) {
    Rcpp::stop("EnvelopeDispersionBuild: implied rate2 <= 0; tilted inverse-gamma is invalid.");
  }
  
  Rcpp::List gamma_list = Rcpp::List::create(
    // Global proposal parameters (unchanged)
    Rcpp::Named("shape3")       = shape3,
    Rcpp::Named("rate2")        = rate2,
    Rcpp::Named("disp_upper")   = upp,
    Rcpp::Named("disp_lower")   = low
  );
  

  List UB_list = List::create(
    Named("RSS_ML")          = RSS_ML,
    Named("RSS_Min")         = rss_min_global,
    Named("max_New_LL_UB")   = max_upp,
    Named("max_LL_log_disp") = lm_log1 + lm_log2 * std::log(upp),
    Named("lm_log1")         = lm_log1,
    Named("lm_log2")         = lm_log2,
    Named("lg_prob_factor")  = lg_prob_factor,
    Named("lmc1")            = new_int,
    Named("lmc2")            = new_slope,
    Named("UB2min")          = ub2_min

  );
  
  
  
  List diagnostics = List::create(
    Named("shape2")          = shape2,
    Named("rate3")           = Rate,
    Named("shape3")          = shape3,
    Named("max_low")         = max_low,
    Named("max_upp")         = max_upp,
    Named("new_slope")       = new_slope,
    Named("new_int")         = new_int,
    Named("UB2min")          = ub2_min
  );
  
  if (!R_finite(shape3) || shape3 <= 0.0) {
    Rcpp::stop("EnvelopeDispersionBuild: implied shape3 <= 0; tilted inverse-gamma is invalid.");
  }
  
    
  return List::create(
    Named("Env")         = Env,
    Named("gamma_list")  = gamma_list,
    Named("UB_list")     = UB_list,
    Named("diagnostics") = diagnostics
  );
}





namespace glmbayes {

namespace env {
List EnvelopeDispersionBuild(
    List Env,
    double Shape,
    double Rate,
    NumericMatrix P,
    NumericVector y,
    NumericMatrix x,
    NumericVector alpha,
    int n_obs,
    double RSS_post,
    double RSS_ML,
    NumericMatrix mu,         // ← new
    NumericVector wt,         // ← new
    double max_disp_perc ,
    Nullable<double> disp_lower ,
    Nullable<double> disp_upper ,
    bool verbose ,
    bool use_parallel    // ← add flag here
  
)
  {
  
  
  // int disp_grid_type=2;
  // 
  // if(use_parallel) disp_grid_type=2;
  
  
  
  // --- NEW: selector for RSS source ---
  // 1 = use minimization (default)
  // 2 = use RSS_ML (skip minimization)
  // int RSS_Min_Type = 1;  // change manually for testing (currently unused: RSS minimization disabled)
  // int UB2_Min_Type = 1;  // change manually for testing (currently unused: UB2 minimization disabled)
  
  double n_w = 0.0;
  for (int i = 0; i < wt.size(); ++i)    n_w += wt[i];
  
  // Step 1: Posterior Gamma parameters (precision prior)
  double shape2 = Shape + n_w / 2.0;
  double rate3  = Rate  + RSS_post / 2.0;
  
  // Step 2: Dispersion bounds (on sigma^2)
  double low, upp;
  if (disp_lower.isNull() || disp_upper.isNull()) {
    // Call R's qgamma for tail quantiles, then invert to get sigma^2 bounds
    Function qgamma("qgamma");
    NumericVector q_low = qgamma(
      Named("p")     = max_disp_perc,
      Named("shape") = shape2,
      Named("rate")  = rate3
    );
    NumericVector q_upp = qgamma(
      Named("p")     = 1.0 - max_disp_perc,
      Named("shape") = shape2,
      Named("rate")  = rate3
    );
    low = 1.0 / q_low[0];
    upp = 1.0 / q_upp[0];
  } else {
    low = as<double>(disp_lower);
    upp = as<double>(disp_upper);
    if (!R_finite(low) || !R_finite(upp))
      stop("disp_lower/disp_upper must be finite.");
    if (low <= 0.0 || upp <= 0.0)
      stop("disp_lower/disp_upper must be positive.");
    if (upp <= low)
      stop("disp_upper must be strictly greater than disp_lower.");
  }
  
  // Step 3: Extract envelope faces
  NumericMatrix cbars     = Env["cbars"];      // gs x l1
  NumericMatrix thetabars = Env["thetabars"];  // gs x l1 (grid of tangencies)
  NumericVector logP1     = Env["logP"];       // length gs
  int gs = cbars.nrow();
  // int l1 = cbars.ncol();  // unused after disabling UB2 minimization
  
  /// Step 3B: Precompute elements for finding inverse function for cbars
  
  
  if (verbose) {
    
    Rcpp::Rcout << "[EnvelopeDipsersionBuild:Inv_f3_precompute_disp] Entering: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  
  Rcpp::List cache = Inv_f3_precompute_disp(cbars, y, x, mu, P, alpha, wt);

  if (verbose) {
    
    Rcpp::Rcout << "[EnvelopeDipsersionBuild:Inv_f3_precompute_disp] Exiting: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  // Step 3B1: Compute RSS_ML when caller did not supply it (rIndepNormalGammaReg passes NA)
  if (!R_finite(RSS_ML)) {
    arma::mat base_A  = cache["base_A"];
    arma::vec base_B0 = cache["base_B0"];
    arma::vec beta_hat = -arma::solve(arma::mat(base_A), base_B0);
    arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
    arma::vec yv(y.begin(), y.size(), false);
    arma::vec alphav(alpha.begin(), alpha.size(), false);
    arma::vec wv(wt.begin(), wt.size(), false);
    arma::vec resid = yv - X * beta_hat - alphav;
    RSS_ML = arma::as_scalar(resid.t() * (wv % resid));
  }

  // Step 3B1.5: Precompute A_max, M_min, A_min, M_max once so they can be
  // shared by both RSS and UB2 bounding logic (even though the current
  // implementations still compute their own versions internally).
  {
    arma::mat base_A  = cache["base_A"];
    arma::mat Pmat    = cache["Pmat"];
    Pmat = 0.5 * (Pmat + Pmat.t());

    // A_max and M_min at d = low
    arma::mat A_max = Pmat + base_A / low;
    A_max = 0.5 * (A_max + A_max.t());
    arma::mat Ainv_max = arma::inv_sympd(A_max);
    arma::mat M_min = Ainv_max.t() * base_A * Ainv_max;
    M_min = 0.5 * (M_min + M_min.t());

    // A_min and M_max at d = upp
    arma::mat A_min = Pmat + base_A / upp;
    A_min = 0.5 * (A_min + A_min.t());
    arma::mat Ainv_min = arma::inv_sympd(A_min);
    arma::mat M_max = Ainv_min.t() * base_A * Ainv_min;
    M_max = 0.5 * (M_max + M_max.t());

    // Store in cache for potential reuse by bounding functions
    cache["A_max"] = A_max;
    cache["M_min"] = M_min;
    cache["A_min"] = A_min;
    cache["M_max"] = M_max;
  }

  // Step 3B2: Closed-form RSS_min bound (dispersion-aware ball)
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeDispersionBuild:bound_rss] Entering: "
                << glmbayes::progress::timestamp_cpp()
                << "\n";
  }
  Rcpp::List rss_bound_res = bound_rss_over_dispersion(
    cache, Env, RSS_ML, shape2, rate3, low, upp, verbose
  );
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeDispersionBuild:bound_rss] Exiting: "
                << glmbayes::progress::timestamp_cpp()
                << "\n";
  }
  Rcpp::NumericVector rss_bound_parallel;
  bool bound_ok = false;
  if (rss_bound_res.containsElementNamed("ok") && Rcpp::as<bool>(rss_bound_res["ok"])) {
    rss_bound_parallel = Rcpp::as<Rcpp::NumericVector>(rss_bound_res["rss_bound_parallel"]);
    bound_ok           = true;
  }

  double rss_min_global;
  NumericVector rss_min_parallel;
  NumericVector disp_min_parallel;
  if (!bound_ok || rss_bound_parallel.size() != static_cast<R_xlen_t>(gs)) {
    Rcpp::stop("bound_rss_over_dispersion failed or size mismatch; RSS must come from closed-form bound.");
  }
  rss_min_global    = Rcpp::as<double>(rss_bound_res["rss_min_global"]);
  rss_min_parallel  = Rcpp::as<Rcpp::NumericVector>(rss_bound_res["rss_min_parallel"]);
  disp_min_parallel = Rcpp::as<Rcpp::NumericVector>(rss_bound_res["disp_min_parallel"]);

  // [COMMENTED OUT] Temporary check + Comparison + quadratic test (require minimize_rss)
  /*
  // Temporary check: bound <= rss_min from minimization for all faces, and diff is small
  if (verbose && gs <= 81) {
    NumericVector rss_min_from_minimize = rss_res["rss_min_parallel"];
    const double tol = 1e-6 * std::max(1.0, rss_min_global);
    bool ok = true;
    for (int j = 0; j < gs; ++j) {
      double diff = rss_min_parallel[j] - rss_min_from_minimize[j];
      if (diff > tol) {
        Rcpp::Rcout << "[rss_bound:CHECK] face " << j << " bound > min: bound="
                    << rss_min_parallel[j] << " min=" << rss_min_from_minimize[j]
                    << " diff=" << diff << "\n";
        ok = false;
      }
    }
    if (!ok) {
      Rcpp::Rcout << "[rss_bound:CHECK] bound should be <= min for all faces; diff should be small.\n";
    }
  }
  // Comparison: bound (source) vs minimize (check), when gs <= 81 and verbose
  if (gs <= 81 && verbose) {
    NumericVector rss_min_from_minimize = rss_res["rss_min_parallel"];
    double rss_min_global_from_minimize = rss_res["rss_min_global"];
    double d1_star_comp = rate3 / (shape2 - 1.0);
    Rcpp::Rcout << "[EnvelopeDispersionBuild:RSS comparison] gs=" << gs
                << " | bound_global=" << rss_min_global
                << " min_global=" << rss_min_global_from_minimize
                << " d1_star=" << d1_star_comp << "\n";
    Rcpp::Rcout << "  face | bound | min | disp_min"
                << " | bound-RSS_ML | min-RSS_ML | bound-min\n";
    for (int j = 0; j < gs; ++j) {
      double diff_bound = rss_min_parallel[j] - RSS_ML;
      double diff_min   = rss_min_from_minimize[j] - RSS_ML;
      double diff_check = rss_min_parallel[j] - rss_min_from_minimize[j];
      Rcpp::Rcout << "  " << j << " | "
                  << rss_min_parallel[j] << " | "
                  << rss_min_from_minimize[j] << " | "
                  << disp_min_parallel[j] << " | "
                  << diff_bound << " | " << diff_min << " | " << diff_check << "\n";
    }
    // Per-face diagnostic: rss_face_bound_from_cache_cpp + rss_face_quadratic_sum_internal (quad vs quad_lb)
    if (gs > 0) {
      Rcpp::Rcout << "  face | bound | actual | diff\n";
      arma::mat X_mat(x.begin(), x.nrow(), x.ncol(), false);
      arma::vec yv(y.begin(), y.size(), false);
      arma::vec alphav(alpha.begin(), alpha.size(), false);
      arma::vec wv(wt.begin(), wt.size(), false);
      arma::mat Q = X_mat.t() * (X_mat.each_col() % wv);
      Q = 0.5 * (Q + Q.t());
      arma::vec rhs = X_mat.t() * (wv % (yv - alphav));
      arma::vec beta_hat_vec = arma::solve(Q, rhs);
      arma::vec resid_ml = (yv - alphav - X_mat * beta_hat_vec) % arma::sqrt(wv);
      double rss_at_beta_hat = arma::as_scalar(resid_ml.t() * resid_ml);
      double rel_err = std::abs(rss_at_beta_hat - RSS_ML) / std::max(1e-15, std::abs(RSS_ML));
      if (rel_err > 1e-8) {
        Rcpp::Rcout << "[rss_face_bound:CHECK] beta_hat does not reproduce RSS_ML: "
                    << "RSS(beta_hat)=" << rss_at_beta_hat << " RSS_ML=" << RSS_ML
                    << " rel_err=" << rel_err << "\n";
      }
      Rcpp::NumericVector beta_hat_r = Rcpp::wrap(beta_hat_vec);
      Rcpp::Rcout << "  [quad vs quad_lb at disp_min]\n";
      for (int j = 0; j < gs; ++j) {
        Rcpp::NumericVector cbars_j_r(cbars.ncol());
        for (int k = 0; k < cbars.ncol(); ++k) cbars_j_r[k] = cbars(j, k);
        rss_face_bound_from_cache_cpp(
          cache, beta_hat_r, cbars_j_r,
          low, upp, y, x, alpha, wt,
          shape2, rate3, RSS_ML, rss_min_parallel[j], true, j
        );
        double disp_j = disp_min_parallel[j];
        Rcpp::List quad_res = rss_face_quadratic_sum_internal(
          cache, cbars_j_r, beta_hat_r, low, upp
        );
        double quad_lb = Rcpp::as<double>(quad_res["quad_lower_bound"]);
        double rss_at_d = rss_face_at_disp(disp_j, cache, cbars_j_r, y, x, alpha, wt);
        double quad = rss_at_d - RSS_ML;
        Rcpp::Rcout << "  face " << j << " | quad=" << quad << " | quad_lb=" << quad_lb
                    << " | diff=" << (quad - quad_lb) << "\n";
      }
    }
  }
  */

  if (verbose) {
    Rcpp::Rcout << "[EnvelopeDispersionBuild:bound_ub2] Entering: "
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  Rcpp::List bound_res = bound_ub2_over_dispersion(
    gs, low, upp,
    cache, cbars,
    y, x, alpha, wt,
    rss_min_global
  );

  if (verbose) {
    Rcpp::Rcout << "[EnvelopeDispersionBuild:bound_ub2] Exiting: "
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
    if (bound_res.containsElementNamed("kappa_K")) {
      double kappa_K_report = Rcpp::as<double>(bound_res["kappa_K"]);
      bool near_iso_report  = Rcpp::as<bool>(bound_res["K_is_near_isotropic"]);
      Rcpp::Rcout << "[EnvelopeDispersionBuild:bound_ub2] kappa(K) = "
                  << kappa_K_report
                  << ", K_is_near_isotropic = " << (near_iso_report ? "TRUE" : "FALSE")
                  << "\n";
    }
  }

  NumericVector disp_min_ub2_bound = bound_res["disp_min_ub2"];

  // UB2 minimization over dispersion is currently disabled; the endpoint
  // evaluation in bound_ub2_over_dispersion is used as the UB2 source.

  // Results from bounding (used downstream)
  NumericVector ub2_min_bound = bound_res["ub2_min"];
  NumericVector disp_min_ub2  = disp_min_ub2_bound;
  NumericVector ub2_min       = ub2_min_bound;

    if (verbose) {
      
      Rcpp::Rcout << "[EnvelopeDipsersionBuild:compute_geometry] Entering: "
      //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << glmbayes::progress::timestamp_cpp()
                    << "\n";
    }
    

    Rcpp::List geom;

  geom = compute_envelope_geometry_cpp(
    cbars,
    thetabars,
    y,
    x,
    P,
    alpha,
    low,
    upp,
    shape2,
    rate3,
    n_w
  );

      
    if (verbose) {
      
      Rcpp::Rcout << "[EnvelopeDipsersionBuild:compute_geometry] Exiting: "
      //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << glmbayes::progress::timestamp_cpp()
                    << "\n";
    }
    
    
  NumericVector thetabar_const_low_apprx = geom["thetabar_const_low_apprx"];
  NumericVector thetabar_const_upp_apprx = geom["thetabar_const_upp_apprx"];

  double max_low  = geom["max_low"];
  double max_upp  = geom["max_upp"];
  double new_slope = geom["new_slope"];
  double new_int   = geom["new_int"];
  // NEW: extract both dispersion anchors
  double d2_star = geom["d2_star"];   // log‑tilt anchor (UB3B geometry)
  
  
  

  ////////////////////////////////////////////////////////////
  
  if (verbose) {
    
    Rcpp::Rcout << "[EnvelopeDipsersionBuild:compute_mixture_outputs] Entering: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  Rcpp::List mix;
  

  mix = compute_mixture_and_outputs_cpp(
    Env,                              // ← pass existing envelope
    thetabar_const_low_apprx,
    thetabar_const_upp_apprx,
    ub2_min,
    logP1,
    max_low,
    max_upp,
    new_slope,
    new_int,
    d2_star,
    shape2,
    Rate,
    low,
    upp,
    RSS_ML,
    rss_min_global,
    verbose
  );
  

  if (verbose) {
    
    Rcpp::Rcout << "[EnvelopeDipsersionBuild:compute_mixture_outputs] Exiting: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  
  Env         = mix["Env"];
  List gamma_list  = mix["gamma_list"];
  List UB_list     = mix["UB_list"];
  List diagnostics = mix["diagnostics"];

  double kappa_K_diag            = bound_res.containsElementNamed("kappa_K")
    ? Rcpp::as<double>(bound_res["kappa_K"]) : NA_REAL;
  bool   K_is_near_isotropic_diag = bound_res.containsElementNamed("K_is_near_isotropic")
    ? Rcpp::as<bool>(bound_res["K_is_near_isotropic"]) : false;
  diagnostics["kappa_K"]             = kappa_K_diag;
  diagnostics["K_is_near_isotropic"] = K_is_near_isotropic_diag;

  return List::create(
    Named("Env_out")    = Env,
    Named("gamma_list") = gamma_list,
    Named("UB_list")    = UB_list,
    Named("diagnostics")= diagnostics
  );
}


}

}


