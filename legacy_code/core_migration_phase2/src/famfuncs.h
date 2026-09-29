// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file famfuncs.h
 * @brief Family-specific likelihood and envelope components for glmbayes.
 *
 * @namespace glmbayes::fam
 * @brief Core log-likelihood, log-posterior, and corresponding log gradient components for each GLM family.
 *
 * @section ImplementedIn
 *   These declarations are implemented in:
 *     - famfuncs_binomial.cpp
 *     - famfuncs_Gamma.cpp
 *     - famfuncs_gaussian.cpp
 *     - famfuncs_poisson.cpp
 *
 * @section UsedBy
 *   These functions are consumed by:
 *     - Envelopefuncs.cpp
 *     - EnvelopeEval.cpp
 *     - EnvelopeBuild.cpp
 *     - rNormalGLM.cpp
 *     - rIndepNormalGammaReg.cpp
 *
 * @section Responsibilities
 *   Provides f1, f2, f3 components for all GLM families and link functions.
 *   Supplies both standard and parallel-safe (RMatrix/RVector) variants.
 *   
 *   f1 functions --> Negative log-likelihoods
 *   f2 functions --> Negative log-likelihood+ small prior component (used for models in standard form)
 *   f3 functions --> Gradients for negative log-likelihoods + small prior components (used for models in standard form)
 *   inv_f3       --> Functions used to compute inverse gradients (with respect to dispersion) for gaussian families
 *
 *   Functions are used during both envelope construction and simulation
 *   
 */


#ifndef GLMBAYES_FAM_H
#define GLMBAYES_FAM_H

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;




namespace glmbayes{

namespace fam {

NumericVector dbinom_glmb( NumericVector x, NumericVector N, NumericVector means, int lg);
NumericVector dpois_glmb( NumericVector x, NumericVector means, int lg);
NumericVector dgamma_glmb( NumericVector x, NumericVector shape, NumericVector scale, int lg);
NumericVector dnorm_glmb( NumericVector x, NumericVector means, NumericVector sds,int lg);

void neg_dpois_glmb_rmat(const RVector<double>& x,     // observed counts
                         const std::vector<double>& means, // Poisson rates
                         std::vector<double>& res,         // output buffer (preallocated)
                         const int lg);                     // log=TRUE?


//----------------- binomial_Logit -------------------------------------------------


NumericVector  f1_binomial_logit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
arma::vec f2_binomial_logit_rmat(
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    int progbar);
arma::mat  f3_binomial_logit(NumericMatrix b,
                             NumericVector y, 
                             NumericMatrix x,
                             NumericMatrix mu,
                             NumericMatrix P,
                             NumericVector alpha,
                             NumericVector wt,
                             int progbar);

// Combined f2/f3 evaluator (binomial–logit)
// Computes negative log‑posterior (f2) and gradient (f3) in one pass.
// Returns a list with:
//   $qf   : NumericVector length = m1   (same as f2_binomial_logit)
//   $grad : NumericMatrix m1 × l2       (same as f3_binomial_logit)
Rcpp::List f2_f3_binomial_logit(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);

//----------------- binomial_Probit -------------------------------------------------
  
NumericVector  f1_binomial_probit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
arma::vec f2_binomial_probit_rmat(
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    int progbar);
arma::mat  f3_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

Rcpp::List f2_f3_binomial_probit(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);

//----------------- binomial_cloglog -------------------------------------------------
  

NumericVector  f1_binomial_cloglog(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);

    arma::vec f2_binomial_cloglog_rmat(
      const RMatrix<double>& b,
      const RVector<double>& y,
      const RMatrix<double>& x,
      const RMatrix<double>& mu,
      const RMatrix<double>& P,
      const RVector<double>& alpha,
      const RVector<double>& wt,
      int progbar);
arma::mat  f3_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);


Rcpp::List f2_f3_binomial_cloglog(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);
  
//----------------- Poisson -------------------------------------------------

NumericVector  f1_poisson(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
arma::vec f2_poisson_rmat(const RMatrix<double>& b,       // candidate coefficients
                             const RVector<double>& y,       // observed counts
                             const RMatrix<double>& x,       // design matrix
                             const RMatrix<double>& mu,      // mode vector
                             const RMatrix<double>& P,       // precision matrix
                             const RVector<double>& alpha,   // predictor offset
                             const RVector<double>& wt,      // observation weights
                             int progbar );          // progress toggle

arma::mat  f3_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

Rcpp::List f2_f3_poisson(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);

//----------------- Gamma -------------------------------------------------
  
NumericVector  f1_gamma(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
arma::vec f2_gamma_rmat(const RMatrix<double>& b,       // candidate coefficients
                          const RVector<double>& y,       // observed counts
                          const RMatrix<double>& x,       // design matrix
                          const RMatrix<double>& mu,      // mode vector
                          const RMatrix<double>& P,       // precision matrix
                          const RVector<double>& alpha,   // predictor offset
                          const RVector<double>& wt,      // observation weights
                          int progbar );          // progress toggle

arma::mat  f3_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

Rcpp::List f2_f3_gamma(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);

//----------------- Gaussian -------------------------------------------------
  
NumericVector  f1_gaussian(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_gaussian(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
arma::vec f2_gaussian_rmat(const RcppParallel::RMatrix<double>& b,       // candidate coefficients
                           const RcppParallel::RVector<double>& y,       // observed counts
                           const RcppParallel::RMatrix<double>& x,       // design matrix
                           const RcppParallel::RMatrix<double>& mu,      // mode vector
                           const RcppParallel::RMatrix<double>& P,       // precision matrix
                           const RcppParallel::RVector<double>& alpha,   // predictor offset
                           const RcppParallel::RVector<double>& wt,      // observation weights
                           int progbar);                                 // progress toggle

// New parallel-safe variant (wt as RMatrix)
arma::vec f2_gaussian_rmat_mat(const RcppParallel::RMatrix<double>& b,   // candidate coefficients
                               const RcppParallel::RVector<double>& y,   // observed counts
                               const RcppParallel::RMatrix<double>& x,   // design matrix
                               const RcppParallel::RMatrix<double>& mu,  // mode vector
                               const RcppParallel::RMatrix<double>& P,   // precision matrix
                               const RcppParallel::RVector<double>& alpha, // predictor offset
                               const RcppParallel::RMatrix<double>& wt,  // observation weights (matrix view)
                               int progbar);                             // progress toggle

arma::mat  f3_gaussian(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);


Rcpp::List f2_f3_gaussian(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);

//----------------- Inv_f3 --Gaussian -------------------------------------------------
  
Rcpp::List Inv_f3_precompute_disp(NumericMatrix cbars,
                                  NumericVector y,
                                  NumericMatrix x,
                                  NumericMatrix mu,
                                  NumericMatrix P,
                                  NumericVector alpha,
                                  NumericVector wt);

// Dispersion-aware envelope solver
arma::mat Inv_f3_with_disp(Rcpp::List cache,
                           double dispersion,
                           Rcpp::NumericMatrix cbars_small);


arma::mat Inv_f3_with_disp_rmat(
    const RcppParallel::RMatrix<double>& Pmat_r,
    const RcppParallel::RMatrix<double>& Pmu_r,
    const RcppParallel::RVector<double>& base_B0_r,
    const RcppParallel::RMatrix<double>& base_A_r,
    double dispersion,
    const RcppParallel::RMatrix<double>& cbars_r // p × m
);


} // famfuncs

} //glmbayes

#endif