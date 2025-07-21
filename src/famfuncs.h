// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

NumericVector dbinom_glmb( NumericVector x, NumericVector N, NumericVector means, int lg);
NumericVector  f1_binomial_logit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
// arma::vec  f2_binomial_logit_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
arma::vec f2_binomial_logit_rmat(
    // NumericMatrix b, NumericVector y,
    //                              NumericMatrix x, NumericMatrix mu,
    //                              NumericMatrix P, NumericVector alpha,
    //                              NumericVector wt, int progbar = 0
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    int progbar);

arma::mat  f3_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

NumericVector  f1_binomial_probit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
// arma::vec  f2_binomial_probit_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
arma::vec f2_binomial_probit_rmat(
    // NumericMatrix b, NumericVector y,
    //                              NumericMatrix x, NumericMatrix mu,
    //                              NumericMatrix P, NumericVector alpha,
    //                              NumericVector wt, int progbar = 0
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    int progbar);

arma::mat  f3_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

NumericVector  f1_binomial_cloglog(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
// arma::vec  f2_binomial_cloglog_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
arma::vec f2_binomial_cloglog_rmat(
    // NumericMatrix b, NumericVector y,
    //                              NumericMatrix x, NumericMatrix mu,
    //                              NumericMatrix P, NumericVector alpha,
    //                              NumericVector wt, int progbar = 0
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    int progbar);


arma::mat  f3_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

void neg_dpois_glmb_rmat(const RVector<double>& x,     // observed counts
                         const std::vector<double>& means, // Poisson rates
                         std::vector<double>& res,         // output buffer (preallocated)
                         const int lg);                     // log=TRUE?
  
NumericVector dpois_glmb( NumericVector x, NumericVector means, int lg);
NumericVector  f1_poisson(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
// arma::vec   f2_poisson_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
// arma::vec f2_poisson_rmat_old(const NumericMatrix& b,
//                           const NumericVector& y,
//                           const NumericMatrix& x,
//                           const NumericMatrix& mu,
//                           const NumericMatrix& P,
//                           const NumericVector& alpha,
//                           const NumericVector& wt,
//                           int progbar);

arma::vec f2_poisson_rmat(const RMatrix<double>& b,       // candidate coefficients
                             const RVector<double>& y,       // observed counts
                             const RMatrix<double>& x,       // design matrix
                             const RMatrix<double>& mu,      // mode vector
                             const RMatrix<double>& P,       // precision matrix
                             const RVector<double>& alpha,   // predictor offset
                             const RVector<double>& wt,      // observation weights
                             int progbar );          // progress toggle

arma::mat  f3_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);


NumericVector dgamma_glmb( NumericVector x, NumericVector shape, NumericVector scale, int lg);
NumericVector  f1_gamma(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
// arma::vec  f2_gamma_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
arma::vec f2_gamma_rmat(const RMatrix<double>& b,       // candidate coefficients
                          const RVector<double>& y,       // observed counts
                          const RMatrix<double>& x,       // design matrix
                          const RMatrix<double>& mu,      // mode vector
                          const RMatrix<double>& P,       // precision matrix
                          const RVector<double>& alpha,   // predictor offset
                          const RVector<double>& wt,      // observation weights
                          int progbar );          // progress toggle

arma::mat  f3_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

NumericVector dnorm_glmb( NumericVector x, NumericVector means, NumericVector sds,int lg);
NumericVector  f1_gaussian(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_gaussian(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
arma::vec  f2_gaussian_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
arma::mat  f3_gaussian(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
NumericVector RSS(NumericVector y, NumericMatrix x,NumericMatrix b,NumericVector alpha,NumericVector wt);
arma::mat  Inv_f3_gaussian(NumericMatrix cbars,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
