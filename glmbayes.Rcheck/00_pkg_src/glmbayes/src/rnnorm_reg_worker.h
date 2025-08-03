

#pragma once

// Required headers
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/mutex.h>
#include <string>
#include "rng_utils.h"  // for safe_runif()

using namespace Rcpp;
using namespace RcppParallel;

//-----------------------------------------------------------------------------
// rnnorm_reg_worker: parallel sampler with envelope logic
//-----------------------------------------------------------------------------
struct rnnorm_reg_worker : public RcppParallel::Worker {
  // --- Inputs ---
  int n;
  
  RVector<double>       y_r;       // observed counts
  RMatrix<double>       x_r;       // design matrix
  RMatrix<double>       mu_r;      // mode matrix
  RMatrix<double>       P_r;       // precision matrix
  RVector<double>       alpha_r;   // predictor offset
  RVector<double>       wt_r;      // observation weights
  
  arma::vec             PLSD;      // slice density
  arma::vec             LLconst;   // envelope constants
  arma::mat             loglt;     // envelope lower bounds
  arma::mat             logrt;     // envelope upper bounds
  arma::mat             cbars;     // envelope centers
  
  CharacterVector       family;    // GLM family
  CharacterVector       link;      // link function
  int                   progbar;   // progress bar toggle
  
  // --- Outputs ---
  RMatrix<double>       out;       // accepted draws
  RVector<double>       draws;     // trial counts
  int                   ncol;      // dimensionality
  
  // --- Constructor ---
  rnnorm_reg_worker(
    int n_,
    const RVector<double>& y_r_,
    const RMatrix<double>& x_r_,
    const RMatrix<double>& mu_r_,
    const RMatrix<double>& P_r_,
    const RVector<double>& alpha_r_,
    const RVector<double>& wt_r_,
    const arma::vec& PLSD_,
    const arma::vec& LLconst_,
    const arma::mat& loglt_,
    const arma::mat& logrt_,
    const arma::mat& cbars_,
    const CharacterVector& family_,
    const CharacterVector& link_,
    int progbar_,
    RMatrix<double>& out_,
    RVector<double>& draws_
  )
    : n(n_),
      y_r(y_r_), x_r(x_r_), mu_r(mu_r_), P_r(P_r_),
      alpha_r(alpha_r_), wt_r(wt_r_),
      PLSD(PLSD_), LLconst(LLconst_),
      loglt(loglt_), logrt(logrt_), cbars(cbars_),
      family(family_), link(link_), progbar(progbar_),
      out(out_), draws(draws_), ncol(out_.ncol())
  {}
  
  // --- Parallel Loop ---
  void operator()(std::size_t begin, std::size_t end);
};
