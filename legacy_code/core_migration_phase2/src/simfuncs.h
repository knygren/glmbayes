// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file simfuncs.h
 * @brief Simulation and posterior sampling routines for glmbayes.
 *
 * @namespace glmbayes::sim
 * @brief Core Normal, Normal–Gamma, and independent Normal–Gamma samplers.
 *
 * @section ImplementedIn
 *   These declarations are implemented in:
 *     - rNormalGLM.cpp
 *     - rNormalGLMBlocks.cpp
 *     - rIndepNormalGammaReg.cpp
 *     - rNormalGammaReg.cpp
 *     - rGammaGaussian.cpp
 *     - rGammaGamma.cpp
 *     - glmb_Standardize_Model.cpp
 *
 * @section UsedBy
 *   These functions are consumed by:
 *     - export_wrappers.cpp 
 *
 * @section Responsibilities
 *   Provides simulation kernels for:
 *     - Normal GLM posterior draws (standardized and unstandardized)
 *     - Normal–Gamma regression models
 *     - Independent Normal–Gamma regression models (standard and parallel variants)
 *     - Gamma–Gaussian and Gamma–Gamma models
 *
 *   All routines:
 *     - assume validated inputs from R wrappers,
 *     - use Rcpp::List for structured return objects,
 *   
 *   Some routines:
 *     - rely on envelope objects and f2/f3 functions for accept–reject sampling,
 *     - support optional parallelization and OpenCL where applicable.
 */

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef GLMBAYES_SIM_H
#define GLMBAYES_SIM_H


// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;


namespace glmbayes {

namespace sim {

Rcpp::List  rNormalGLM_std(int n,
                               NumericVector y,
                               NumericMatrix x,
                               NumericMatrix mu,
                               NumericMatrix P,
                               NumericVector alpha,
                               NumericVector wt,
                               Function f2,
                               Rcpp::List  Envelope,
                               Rcpp::CharacterVector   family,
                               Rcpp::CharacterVector   link, 
                               int progbar=1,
                               bool verbose = false                                 
                                 );

Rcpp::List rNormalGLM(
    int n,
    NumericVector y,
    NumericMatrix x,
    NumericVector mu,
    NumericMatrix P,
    NumericVector offset,
    NumericVector wt,
    double dispersion,
    Function f2,
    Function f3,
    NumericVector start,
    std::string family = "binomial",
    std::string link = "logit",
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
);

/// C++ counterpart to R \c rNormalGLM_reg_block(): loop over \c row_blocks,
/// call \c rNormalGLM() per block. Not exported to R until Phase 2.
Rcpp::List rNormalGLMBlocks(
    int n,
    NumericVector y,
    NumericMatrix x,
    NumericVector offset,
    NumericVector wt,
    NumericVector dispersion,
    NumericMatrix mu,
    List P_blocks,
    bool prior_by_block,
    List row_blocks,
    Function f2,
    Function f3,
    std::string family,
    std::string link,
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
);

Rcpp::List rNormalReg(int n,NumericVector y,NumericMatrix x, 
                         NumericVector mu,NumericMatrix P,
                         NumericVector offset,NumericVector wt,
                         double dispersion,
                         Function f2,Function f3,
                         NumericVector start,
                         std::string family="gaussian",
                         std::string link="identity",
                         int Gridtype=2      
);



Rcpp::List  rIndepNormalGammaReg_std(int n,NumericVector y,NumericMatrix x,
                                          NumericMatrix mu, /// This is typically standardized to be a zero vector
                                          NumericMatrix P, /// Part of prior precision shifted to the likelihood
                                          NumericVector alpha,NumericVector wt,
                                          Function f2,Rcpp::List  Envelope,
                                          Rcpp::List  gamma_list,
                                          Rcpp::List  UB_list,
                                          Rcpp::CharacterVector   family,Rcpp::CharacterVector   link,
                                          bool progbar=true,
                                          bool verbose=false);

Rcpp::List rIndepNormalGammaReg_std_parallel(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,   // typically standardized to be a zero vector
    Rcpp::NumericMatrix P,    // part of prior precision shifted to the likelihood
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    Rcpp::Function f2,        // kept for signature parity
    Rcpp::List Envelope,
    Rcpp::List gamma_list,
    Rcpp::List UB_list,
    Rcpp::CharacterVector family,
    Rcpp::CharacterVector link,
    bool progbar = true,
    bool verbose = false
);


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
    Rcpp::Nullable<double> max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    bool verbose
);


Rcpp::List rIndepNormalGammaReg(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector offset,
    Rcpp::NumericVector wt,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<Rcpp::NumericVector> disp_lower,
    Rcpp::Nullable<Rcpp::NumericVector> disp_upper,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    bool progbar
);

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
    bool verbose = false
);


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
    bool verbose = false
);

Rcpp::List glmb_Standardize_Model(
    NumericVector y, 
    NumericMatrix x,   // Original design matrix (to be adjusted)
    NumericMatrix P,   // Prior Precision Matrix (to be adjusted)
    NumericMatrix bstar, // Posterior Mode from optimization (to be adjusted)
    NumericMatrix A1  // Precision for Log-Posterior at posterior mode (to be adjusted)
);


}
}

#endif