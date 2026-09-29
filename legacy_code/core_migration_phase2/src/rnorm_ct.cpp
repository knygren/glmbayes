// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <math.h>
#include "rng_utils.h"


#include <Rmath.h>       // For R::qnorm


// Don't add mutex if EMSCRIPTEN

#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
#include <tbb/mutex.h>   // For thread locking
tbb::mutex qnorm_mutex;  // Local mutex for this file
#endif


using namespace Rcpp;
using namespace glmbayes::rng;


double safe_qnorm_logp(double logp, double mu, double sigma, bool lower_tail) {
#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
    tbb::mutex::scoped_lock lock(qnorm_mutex);
#endif  
  return R::qnorm(logp, mu, sigma, lower_tail, true);  // log.p = TRUE
}

namespace glmbayes {
namespace rng {

    double rnorm_ct(double lgrt, double lglt, double mu, double sigma) {
  double U = 0;
  double out = 0;
  
  if (lgrt >= lglt) {
    U = runif_safe();
    
    // u1 = 1 - exp(lgrt)  →  -expm1(lgrt)
    double u1   = -std::expm1(lgrt);
    double lgu1 = std::log(u1);
    
    // log(1 - exp(lgu1 - lglt)) → log(-expm1(lgu1 - lglt))
    double lgU2 = std::log(U) + lglt + std::log(-std::expm1(lgu1 - lglt));
    double lgU3 = lgU2 + std::log1p(std::exp(lgu1 - lgU2));
    
    out = safe_qnorm_logp(lgU3, mu, sigma, true);
    
  } else {
    U = runif_safe();
    
    // e1mu2 = 1 - exp(lglt) → -expm1(lglt)
    double e1mu2  = -std::expm1(lglt);
    double lg1mu2 = std::log(e1mu2);
    
    // log(1 - exp(lg1mu2 - lgrt)) → log(-expm1(lg1mu2 - lgrt))
    double lgU2 = std::log(U) + lgrt + std::log(-std::expm1(lg1mu2 - lgrt));
    double lgU3 = lgU2 + std::log1p(std::exp(lg1mu2 - lgU2));
    
    out = safe_qnorm_logp(lgU3, mu, sigma, false);
  }
  
  return out;
}

}
}