// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file Envelopefuncs.h
 * @brief Core envelope–construction routines for glmbayes.
 *
 * @namespace glmbayes::env
 * @brief Algorithms for grid construction, envelope building, evaluation,
 *        sorting, and dispersion–aware refinement used in accept–reject
 *        sampling for GLM Bayesian models.
 *
 * @section ImplementedIn
 *   These declarations are implemented in:
 *     - EnvelopeSize.cpp
 *     - EnvelopeBuild.cpp
 *     - EnvelopeEval.cpp
 *     - EnvelopeBuild_Ind_Normal_Gamma.cpp
 *     - EnvelopeDispersionBuild.cpp
 *
 * @section UsedBy
 *   These functions are consumed by:
 *     - EnvelopeOrchestrator.cpp (high‑level orchestration)
 *     - rNormalGLM.cpp and rIndepNormalGammaReg.cpp (simulation routines)
 *     - export_wrappers wrappers that expose envelope construction to the user
 *
 * @section Responsibilities
 *   Provides the computational kernels for:
 *     - grid sizing and initialization (EnvelopeSize)
 *     - envelope construction for standard GLM families (EnvelopeBuild)
 *     - envelope evaluation at grid points (EnvelopeEval)
 *     - independent Normal–Gamma envelope variants (EnvelopeBuild_Ind_Normal_Gamma)
 *     - dispersion‑aware envelope refinement (EnvelopeDispersionBuild)
 *     - envelope sorting and UB‑list assembly (EnvelopeSort_cpp)
 *     - grid and log‑probability updates (EnvelopeSet_Grid, EnvelopeSet_LogP)
 *
 *   These routines:
 *     - assume validated inputs from R wrappers,
 *     - operate on Armadillo matrices and Rcpp vectors,
 *     - support optional OpenCL and parallel execution where applicable,
 *     - form the backbone of accept–reject sampling in glmbayes.
 */

#ifndef GLMBAYES_ENV_H
#define GLMBAYES_ENV_H


// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;


namespace glmbayes{

namespace env{

inline void check_disp_bounds_or_stop(
    double disp_lower,
    double disp_upper,
    const char* fn_name
) {
  if (!R_finite(disp_lower) || !R_finite(disp_upper) ||
      disp_lower <= 0.0 || disp_upper <= 0.0 || disp_upper <= disp_lower) {
    Rcpp::stop(
      "invalid disp_lower or disp_upper in function %s. disp_lower=%g, disp_upper=%g",
      fn_name, disp_lower, disp_upper
    );
  }
}

inline void check_disp_bounds_or_stop(
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    const char* fn_name
) {
  if (disp_lower.isNull() || disp_upper.isNull()) return;
  check_disp_bounds_or_stop(
    Rcpp::as<double>(disp_lower), Rcpp::as<double>(disp_upper), fn_name
  );
}

inline void check_disp_bounds_or_stop(
    Rcpp::Nullable<Rcpp::NumericVector> disp_lower,
    Rcpp::Nullable<Rcpp::NumericVector> disp_upper,
    const char* fn_name
) {
  if (disp_lower.isNull() || disp_upper.isNull()) return;
  Rcpp::NumericVector lo_v = Rcpp::as<Rcpp::NumericVector>(disp_lower);
  Rcpp::NumericVector up_v = Rcpp::as<Rcpp::NumericVector>(disp_upper);
  if (lo_v.size() < 1 || up_v.size() < 1) {
    Rcpp::stop(
      "invalid disp_lower or disp_upper in function %s. disp_lower=<empty>, disp_upper=<empty>",
      fn_name
    );
  }
  check_disp_bounds_or_stop(lo_v[0], up_v[0], fn_name);
}

Rcpp::List EnvelopeSize(const arma::vec& a,
                        const Rcpp::NumericMatrix& G1,
                        int Gridtype   = 2,
                        int n          = 1000,
                        int n_envopt   = -1,
                        bool use_opencl = false,
                        bool verbose    = false);



List EnvelopeBuild(NumericVector bStar,
                       NumericMatrix A,
                       NumericVector y,
                       NumericMatrix x,
                       NumericMatrix mu,
                       NumericMatrix P,
                       NumericVector alpha,
                       NumericVector wt,
                       std::string family = "binomial",
                       std::string link   = "logit",
                       int Gridtype       = 2,
                       int n              = 1,
                       int n_envopt       = -1,   // NEW: effective sample size for EnvelopeOpt (defaults to n if -1)
                       bool sortgrid      = false,
                       bool use_opencl    = false, // Enables OpenCL acceleration during envelope construction
                       bool verbose       = false  // Enables diagnostic output
);


Rcpp::List EnvelopeEval(const Rcpp::NumericMatrix& G4,   // grid (parameters × grid points)
                        const Rcpp::NumericVector& y,
                        const Rcpp::NumericMatrix& x,
                        const Rcpp::NumericMatrix& mu,
                        const Rcpp::NumericMatrix& P,
                        const Rcpp::NumericVector& alpha,
                        const Rcpp::NumericVector& wt,
                        const std::string& family,
                        const std::string& link,
                        bool use_opencl = false,
                        bool verbose = false);

List EnvelopeBuild_Ind_Normal_Gamma(NumericVector bStar,
                                    NumericMatrix A,
                                    NumericVector y, 
                                    NumericMatrix x,
                                    NumericMatrix mu,
                                    NumericMatrix P,
                                    NumericVector alpha,
                                    NumericVector wt,
                                    std::string family="binomial",
                                    std::string link="logit",
                                    int Gridtype=2, 
                                    int n=1,
                                    int n_envopt=-1,
                                    bool sortgrid=false,
                                    bool use_opencl    = false,
                                    bool verbose       = false);


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
    
);


/// EnvelopeCentering: dispersion anchoring for Normal-Gamma; RSS_post is closed-form
/// E[RSS] each iteration and drives the Gamma update.
List EnvelopeCentering(
    NumericVector y,
    NumericMatrix x,
    NumericVector mu,
    NumericMatrix P,
    NumericVector offset,
    NumericVector wt,
    double shape,
    double rate,
    int Gridtype = 2,
    bool verbose = false
);

Rcpp::List EnvelopeOrchestrator(
    NumericVector bstar2,
    NumericMatrix A,
    NumericVector y,
    NumericMatrix x2,
    NumericMatrix mu2,
    NumericMatrix P2,
    NumericVector alpha,
    NumericVector wt,
    
    int n,
    int Gridtype,
    Nullable<int> n_envopt,
    
    double shape,
    double rate,
    double RSS_Post2,
    double RSS_ML,
    
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    
    bool use_parallel,
    bool use_opencl,
    bool verbose
);

Rcpp::List EnvelopeSort_cpp(
    int l1,
    int l2,
    const Rcpp::NumericMatrix& GIndex,
    const Rcpp::NumericMatrix& G3,
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericMatrix& logU,   // l2 x l1 (or l2 x ncol); R also accepts vector
    const Rcpp::NumericMatrix& logrt,
    const Rcpp::NumericMatrix& loglt,
    const Rcpp::NumericMatrix& logP,   // l2 x 2 (or l2 x ncol)
    const Rcpp::NumericMatrix& LLconst,  // l2 x 1 (or l2 x ncol)
    const Rcpp::NumericVector& PLSD,
    const Rcpp::NumericVector& a1,
    double E_draws,
    const Rcpp::Nullable<Rcpp::NumericVector>& lg_prob_factor = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector>& UB2min        = R_NilValue
);



List EnvelopeSet_Grid(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint);
void EnvelopeSet_Grid_C2(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint,Rcpp::NumericMatrix Down,Rcpp::NumericMatrix Up,Rcpp::NumericMatrix lglt,Rcpp::NumericMatrix lgrt,Rcpp::NumericMatrix lgct,Rcpp::NumericMatrix logU,Rcpp::NumericMatrix logP);
void EnvelopeSet_Grid_C2_pointwise(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint,Rcpp::NumericMatrix Down,Rcpp::NumericMatrix Up,Rcpp::NumericMatrix lglt,Rcpp::NumericMatrix lgrt,Rcpp::NumericMatrix lgct,Rcpp::NumericMatrix logU,Rcpp::NumericMatrix logP);


List   EnvelopeSet_LogP(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3);


} //env

}  //glmbayes




NumericVector RSS(NumericVector y, NumericMatrix x,NumericMatrix b,NumericVector alpha,NumericVector wt);




// Rcpp::List Set_Grid_C(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint,Rcpp::NumericMatrix Down,Rcpp::NumericMatrix Up,Rcpp::NumericMatrix lglt,Rcpp::NumericMatrix lgrt,Rcpp::NumericMatrix lgct,Rcpp::NumericMatrix logU,Rcpp::NumericMatrix logP);


// Rcpp::List   setlogP_C(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst);
void setlogP_C2(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst);




double rss_face_at_disp(double dispersion,
                        Rcpp::List cache,
                        Rcpp::NumericVector cbars_j,
                        Rcpp::NumericVector y,
                        Rcpp::NumericMatrix x,
                        Rcpp::NumericVector alpha,
                        Rcpp::NumericVector wt);

double UB2(double dispersion,
           Rcpp::List cache,
           Rcpp::NumericVector cbars_j,
           Rcpp::NumericVector y,
           Rcpp::NumericMatrix x,
           Rcpp::NumericVector alpha,
           Rcpp::NumericVector wt,
           double rss_min_global);



#endif