// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppParallel)]]
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "RcppParallel.h"

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include "simfuncs.h"
#include "progress_utils.h"


#include <math.h>
#include <cmath>         // for std::log or std::exp if used


#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
#include <tbb/mutex.h>
static tbb::mutex f2_mutex;
#endif

#include <thread>
#include "rng_utils.h"

#include <atomic>
#include <memory>
#include <string>


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;
using namespace glmbayes::fam;
using namespace glmbayes::env;
using namespace glmbayes::sim;
using namespace glmbayes::rng;
using namespace glmbayes::progress;


namespace glmbayes {

namespace sim {

Rcpp::List glmb_Standardize_Model(
    NumericVector y, 
    NumericMatrix x,   // Original design matrix (to be adjusted)
    NumericMatrix P,   // Prior Precision Matrix (to be adjusted)
    NumericMatrix bstar, // Posterior Mode from optimization (to be adjusted)
    NumericMatrix A1  // Precision for Log-Posterior at posterior mode (to be adjusted)
  ) {
  
  // Get dimensions for x matrix
  
  int l1=x.ncol();
  int l2=x.nrow();

  // Arma matrix versions of most inputs
    
  arma::mat x2(x.begin(), l2, l1, false);
  arma::mat P2(P.begin(), P.nrow(), P.ncol(), false);
  arma::mat b2(bstar.begin(), bstar.nrow(), bstar.ncol(), false);
  arma::mat A1_b(A1.begin(), l1, l1, false); 
  
  // For now keep these name (legacy from older code inside glmbsim_NGauss
  
  NumericMatrix A4_1(l1, l1);  
  NumericMatrix b4_1(l1,1);  // Maybe this can be set as vector to avoid the below conversi
  NumericMatrix x4_1(l2, l1);
  NumericMatrix mu4_1(l1,1);  
  NumericMatrix P5_1(l1, l1);
  NumericMatrix P6Temp_1(l1, l1);
  NumericMatrix L2Inv_1(l1, l1); 
  NumericMatrix L3Inv_1(l1, l1); 

  double scale=1;
  int check=0;
  double eigval_temp;
  
  Rcpp::Function asVec("as.vector");

  arma::mat L2Inv(L2Inv_1.begin(), L2Inv_1.nrow(), L2Inv_1.ncol(), false);
  arma::mat L3Inv(L3Inv_1.begin(), L3Inv_1.nrow(), L3Inv_1.ncol(), false);
  arma::mat b4(b4_1.begin(), b4_1.nrow(), b4_1.ncol(), false);
//  arma::mat mu4(mu4_1.begin(), mu4_1.nrow(), mu4_1.ncol(), false);
  arma::mat x4(x4_1.begin(), x4_1.nrow(), x4_1.ncol(), false);
  arma::mat A4(A4_1.begin(), A4_1.nrow(), A4_1.ncol(), false);
  arma::mat P5(P5_1.begin(), P5_1.nrow(), P5_1.ncol(), false);
  arma::mat P6Temp(P6Temp_1.begin(), P6Temp_1.nrow(), P6Temp_1.ncol(), false);
  arma::vec eigval_1;
  arma::mat eigvec_1;
  arma::vec eigval_2;
  arma::mat eigvec_2;
  arma::mat ident=arma::mat (l1,l1,arma::fill::eye);
  
// Standardize Model to Have Diagonal Variance-Covariance Matrix at Posterior Mode

      eig_sym(eigval_1, eigvec_1, A1_b);
      
//////////////  Adding warnings if ill-conditioned /////  
      
      double lambda_min = eigval_1.min();
      double lambda_max = eigval_1.max();
      double kappa_H    = lambda_max / lambda_min;
      
      // --- Numerical stability warnings based on κ(H) ---
      
      if (!R_finite(kappa_H)) {
        Rcpp::Rcout <<
          "[glmb_Standardize_Model][WARNING] Posterior Hessian is not finite.\n"
          "  kappa(H) is NaN or Inf.\n"
          "  Standardization is likely to be numerically unstable.\n";
      }
      else if (kappa_H > 1e8) {
        Rcpp::Rcout <<
          "[glmb_Standardize_Model][WARNING] Posterior Hessian is effectively singular.\n"
          "  kappa(H) = " << kappa_H << "\n"
          "  Standardization may be unreliable; curvature is dominated by roundoff.\n";
      }
      else if (kappa_H > 1e6) {
        Rcpp::Rcout <<
          "[glmb_Standardize_Model][WARNING] Posterior Hessian is numerically dangerous.\n"
          "  kappa(H) = " << kappa_H << "\n"
          "  Curvature is extremely uneven; standardization may be unstable.\n";
      }
      else if (kappa_H > 1e5) {
        Rcpp::Rcout <<
          "[glmb_Standardize_Model][WARNING] Posterior Hessian is severely ill-conditioned.\n"
          "  kappa(H) = " << kappa_H << "\n"
          "  Expect sensitivity to rounding and potential instability.\n";
      }
      else if (kappa_H > 1e4) {
        Rcpp::Rcout <<
          "[glmb_Standardize_Model][NOTE] Posterior Hessian is moderately ill-conditioned.\n"
          "  kappa(H) = " << kappa_H << "\n";
      }
      
      
////////////////////////////////////////////////////////////      
      
      arma::mat D1=arma::diagmat(eigval_1);
      arma::mat L2= arma::sqrt(D1)*trans(eigvec_1);
      L2Inv=eigvec_1*sqrt(inv_sympd(D1));  // Also used to undo normalization later

      // output variables used in latter step

      arma::mat b3=L2*b2;   
//      arma::mat mu3=L2*mu2; // These are needed but will not be used to pass 
      //  arma::mat mu3=L2*mu1b; // Corrected - mu1b is zero vector

      arma::mat x3=x2*L2Inv;
      arma::mat P3=trans(L2Inv)*P2*L2Inv;

      // Seup for loop that sets epsilon
    
      arma::mat P3Diag=arma::diagmat(arma::diagvec(P3));// diagonal part of P3
      arma::mat epsilon=P3Diag;
      arma::mat P4=P3Diag;   

    // Find scaled matrix epsilon


    while(check==0){
      epsilon=scale*P3Diag;  // scaled version of diagonal matrix
  
      // Checks if difference between Prior precision and diagonal matrix
      // is positive definite
      // is positive definite - to be added to likelihood 
      
      
  
      P4=P3-epsilon;				
      eig_sym(eigval_2, eigvec_2, P4);
      eigval_temp=arma::min(eigval_2);
      if(eigval_temp>0){check=1;
    
      //      Rcpp::Rcout << "scale after step1 " << std::flush << scale << std::endl;
      //      P4.print("P4 after step 1");  
      //      epsilon.print("epsilon after step 1");  
      }
      else{scale=scale/2;}
    }


    
      // Setup prior to Eigenvalue decomposition

      // NOTE: Epsilon is a diagonal part of the prior
      // At this point, the Posterior is the identity matrix
      // We keep only a diagnonal subset (epsilon) of the Prior in the prior and shift the rest P3-Epsilon to the likelihood

      arma::mat A3=ident-epsilon;	// This should be a diagonal matrix and represents "data" precision in transformed model

    //   Put into Standard form where prior is identity matrix

    // Perform second eigenvalue decomposition

    eig_sym(eigval_2, eigvec_2, epsilon);
    arma::mat D2=arma::diagmat(eigval_2);
    
    // Apply second transformation

      
    
    
    arma::mat L3= arma::sqrt(D2)*trans(eigvec_2);
    L3Inv=eigvec_2*sqrt(inv_sympd(D2));
    b4=L3*b3; 
//    mu4=L3*mu3; 
    x4=x3*L3Inv;
    A4=trans(L3Inv)*A3*L3Inv;  // Should be transformed data precision matrix
    P5=trans(L3Inv)*P4*L3Inv;  // Should be precision matrix without epsilon
    P6Temp=P5+ident;           // Should be precision matrix for posterior (may not be used)

        // "Oddball extra steps due to legacy code 
    
    NumericVector b5=asVec(b4_1); // Maybe this causes error?
    
    NumericMatrix mu5_1=0*mu4_1; // Does this modify mu4_1? Set mu5_1 to 0 more efficiently

    // Interpretation at end of standardization (Nygren & Nygren 2006, standard form):
    // I     --> Prior for standardized model (identity).
    // P5    --> Part of prior shifted to log-likelihood (modified log-likelihood).
    // I+P5  --> Prior for original full model in standardized coords (equiv. to P passed in).
    // A4    --> Data precision for modified model.
    // A4-P5 --> Data precision for original log-likelihood (equiv. to X'WX in standardized coords).
    // I+A4  --> Posterior precision in standardized space (equiv. to A1 passed in).
    //
    // Special case: When P3 is diagonal and only one "divide by 2" is needed (scale=0.5), we have
    // epsilon = (1/2)*P3 and P4 = P3 - epsilon = (1/2)*P3, so epsilon = P4. The second transform
    // makes both the prior and P4 into I, hence P5 = I. With the default (Zellner-like) prior
    // from Prior_Setup, P3 is diagonal, so P5 equals I and downstream use of P5 vs I cannot be
    // distinguished unless a test prior that yields non-diagonal P3 (or scale < 0.5) is used.

    return Rcpp::List::create(
      Rcpp::Named("bstar2")=b5,       // Transformed posterior mode (untransposed also used)
      Rcpp::Named("A")=A4_1,                 // Precision for Standardized data precision 
      Rcpp::Named("x2")=x4_1,                // Transformed Design matrix
      Rcpp::Named("mu2")=mu5_1,               // Transformed prior mean (should really always be 0)
      Rcpp::Named("P2")=P5_1,               // Precision matrix for Normal component shifted to Log-Likelihood
      Rcpp::Named("L2Inv")=L2Inv,               // Precision matrix for Normal component shifted to Log-Likelihood
      Rcpp::Named("L3Inv")=L3Inv               // Precision matrix for Normal component shifted to Log-Likelihood
    );
  
}

}
}




//-----------------------------------------------------------------------------
// rNormalGLM_worker: parallel sampler with envelope logic
//-----------------------------------------------------------------------------
struct rNormalGLM_worker : public RcppParallel::Worker {
  // --- Inputs ---
  int n;
  
  RVector<double>       y_r;       // observed counts
  RMatrix<double>       x_r;       // design matrix
  RMatrix<double>       mu_r;      // mode matrix
  RMatrix<double>       P_r;       // precision matrix
  RVector<double>       alpha_r;   // predictor offset
  RVector<double>       wt_r;      // observation weights
  
  // Envelope components as thread-safe handles (no copies)
  RVector<double> PLSD_r;
  RVector<double> LLconst_r;
  RMatrix<double> loglt_r;
  RMatrix<double> logrt_r;
  RMatrix<double> cbars_r;
  
  
  //   arma::vec             PLSD;      // slice density
  // arma::vec             LLconst;   // envelope constants
  // arma::mat             loglt;     // envelope lower bounds
  // arma::mat             logrt;     // envelope upper bounds
  // arma::mat             cbars;     // envelope centers
  
  CharacterVector       family;    // GLM family
  CharacterVector       link;      // link function
  int                   progbar;   // progress bar toggle
  
  // --- Outputs ---
  RMatrix<double>       out;       // accepted draws
  RVector<double>       draws;     // trial counts
  int                   ncol;      // dimensionality
  
  // --- Optional test controls ---
  // shared atomic flag: set to 1 by any thread if it hits the cap
  std::shared_ptr<std::atomic<int>> any_maxdraw_flag; // default nullptr (no reporting)
  int                   max_draws;                   // -1 => no per-index cap
  
  // --- Constructor ---
  rNormalGLM_worker(
    int n_,
    const RVector<double>& y_r_,
    const RMatrix<double>& x_r_,
    const RMatrix<double>& mu_r_,
    const RMatrix<double>& P_r_,
    const RVector<double>& alpha_r_,
    const RVector<double>& wt_r_,
    
    const RcppParallel::RVector<double>& PLSD_r_,
    const RcppParallel::RVector<double>& LLconst_r_,
    const RcppParallel::RMatrix<double>& loglt_r_,
    const RcppParallel::RMatrix<double>& logrt_r_,
    const RcppParallel::RMatrix<double>& cbars_r_,
    
    // const arma::vec& PLSD_,
    // const arma::vec& LLconst_,
    // const arma::mat& loglt_,
    // const arma::mat& logrt_,
    // const arma::mat& cbars_,
    const CharacterVector& family_,
    const CharacterVector& link_,
    int progbar_,
    RMatrix<double>& out_,
    RVector<double>& draws_,
    std::shared_ptr<std::atomic<int>> any_maxdraw_flag_ = nullptr, // optional shared flag
    int max_draws_ = -1                                              // optional per-index cap
  )
    : n(n_),
      y_r(y_r_), x_r(x_r_), mu_r(mu_r_), P_r(P_r_),
      alpha_r(alpha_r_), wt_r(wt_r_),
      PLSD_r(PLSD_r_), LLconst_r(LLconst_r_),
      loglt_r(loglt_r_), logrt_r(logrt_r_), cbars_r(cbars_r_),
      // PLSD(PLSD_), LLconst(LLconst_),
      // loglt(loglt_), logrt(logrt_), cbars(cbars_),
      family(family_), link(link_), progbar(progbar_),
      out(out_), draws(draws_), ncol(out_.ncol())
    , any_maxdraw_flag(any_maxdraw_flag_),
    max_draws(max_draws_)
    
    
  {}
  
  // --- Parallel Loop ---
  void operator()(std::size_t begin, std::size_t end);
};
  


// operator() implements the parallel loop
void rNormalGLM_worker::operator()(std::size_t begin, std::size_t end) {
  // Per-call lightweight views; no allocations of K×p
  arma::vec  PLSD   (PLSD_r.begin(),    PLSD_r.length(),               false, true);
  arma::vec  LLconst(LLconst_r.begin(), LLconst_r.length(),            false, true);
  arma::mat  loglt  (loglt_r.begin(),   loglt_r.nrow(), loglt_r.ncol(), false, true);
  arma::mat  logrt  (logrt_r.begin(),   logrt_r.nrow(), logrt_r.ncol(), false, true);
  arma::mat  cbars  (cbars_r.begin(),   cbars_r.nrow(), cbars_r.ncol(), false, true);



  // Create Armadillo views directly from RMatrix/RVector memory
  arma::vec y2(y_r.begin(), y_r.length(), false);
  arma::vec alpha2(alpha_r.begin(), alpha_r.length(), false);
  arma::vec wt2(wt_r.begin(), wt_r.length(), false);

  arma::mat x2(x_r.begin(), x_r.nrow(), x_r.ncol(), false);
  arma::mat mu2(mu_r.begin(), mu_r.nrow(), mu_r.ncol(), false);
  arma::mat P2(P_r.begin(), P_r.nrow(), P_r.ncol(), false);


  // Precompute dimensions and envelope pieces
  int l1 = mu_r.nrow();


  // Convert family/link once per thread
  std::string fam2 = as<std::string>(family);
  std::string lnk2 = as<std::string>(link);


  // Thread‐local buffers and views
  std::vector<double> outtemp_buf(l1), cbartemp_buf(l1);
  arma::rowvec        outtemp2(outtemp_buf.data(),   l1, false);
  arma::rowvec        cbartemp2(cbartemp_buf.data(), l1, false);


  std::vector<double> btemp_buf(l1);
  arma::mat btemp2(btemp_buf.data(), l1, 1, false);
  RcppParallel::RMatrix<double> btemp_r(btemp_buf.data(), l1, 1); // optional: only if still needed


  arma::mat testtemp2(1, 1);  // Allocated directly on the heap
  arma::vec testll2(1, arma::fill::none);  // Uninitialized vector of size m1

  /////////////////////////////////////////////////////////


  //    Rcpp::Rcout << "1.0 Launching Worker: " << begin << std::endl;

  // Main loop over indices
  for (std::size_t i = begin; i < end; ++i) {

    draws[i] = 1.0;


    //      Rcpp::Rcout << "i=" << i  << "\n";

    double a1 = 0.0;

    while (a1 == 0.0) {


      // 1) slice selection
      //double U  = R::runif(0.0, 1.0)
      double U = runif_safe();
      double a2 = 0.0;
      int    J  = 0;
      while (a2 == 0.0) {
        if (U <= PLSD[J]) {
          //if (U <= PLSD2[J]) {
          a2 = 1.0;
        } else {
          U -= PLSD[J];
          ++J;
        }
      }

      // 2) draw truncated‐normal candidates
      for (int j = 0; j < l1; ++j) {
        out(i, j) = rnorm_ct(logrt(J, j),loglt(J, j),-cbars(J, j), 1.0 );
        //  out(i, j) = rnorm_ct(logrt2(J, j),loglt2(J, j),-cbars2(J, j), 1.0 );
      }

      // 3) prepare for test
      for (int j = 0; j < l1; ++j) {
        outtemp_buf[j]  = out(i, j);
        cbartemp_buf[j] = cbars(J, j);
        //cbartemp_buf[j] = cbars2(J, j);


      }
      testtemp2 = outtemp2 * trans(cbartemp2);
      //      double U2 = R::runif(0.0, 1.0);

      double U2 = runif_safe();


      btemp2   = trans(outtemp2);

      // declare test here so it’s in scope below
      //double test;



      // 4) compute log‐lik and print test under lock
      {

#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
        tbb::mutex::scoped_lock lock(f2_mutex);
#endif

        //          std::lock_guard<std::mutex> guard(f2_mutex);



        // compute testll for all families/links
        if (fam2 == "binomial") {
          //            if (lnk2 == "logit")      testll = f2_binomial_logit(btemp,y,x,mu,P,alpha,wt,0);
          if (lnk2 == "logit")
          {
            testll2 = f2_binomial_logit_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);

            //Rcpp::Rcout << "rmat version: " << testll2  << "\n";

            //testll2 = f2_binomial_logit_arma(btemp,y,x,mu,P,alpha,wt,0);

            //              Rcpp::Rcout << "arma version: " << testll2  << "\n";

            //              testll2 = f2_binomial_logit(btemp,y,x,mu,P,alpha,wt,0);
            //              Rcpp::Rcout << "original version: " << testll2  << "\n";

          }
          //if (lnk2 == "logit")      testll = f2_binomial_logit_arma(btemp,y,x,mu,P,alpha,wt,0);

          //else if (lnk2 == "probit") testll = f2_binomial_probit(btemp,y,x,mu,P,alpha,wt,0);
          //                    else if (lnk2 == "probit") testll = f2_binomial_probit_arma(btemp,y,x,mu,P,alpha,wt,0);
          else if (lnk2 == "probit")
          {

            testll2 = f2_binomial_probit_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
            //                      Rcpp::Rcout << "rmat version: " << testll2  << "\n";

            //                      testll2 = f2_binomial_probit_arma(btemp,y,x,mu,P,alpha,wt,0);

            //                      Rcpp::Rcout << "arma version: " << testll2  << "\n";

            //                      testll2 = f2_binomial_probit_arma(btemp,y,x,mu,P,alpha,wt,0);
          }
          //                    else                       testll = f2_binomial_cloglog(btemp,y,x,mu,P,alpha,wt,0);
          //                    else                       testll = f2_binomial_cloglog_arma(btemp,y,x,mu,P,alpha,wt,0);
          else
          {

            testll2 = f2_binomial_cloglog_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
            //                                            Rcpp::Rcout << "rmat version: " << testll2  << "\n";

            //                      testll2 = f2_binomial_cloglog_arma(btemp,y,x,mu,P,alpha,wt,0);

            //                                            Rcpp::Rcout << "arma version: " << testll2  << "\n";

          }
        }
        else if (fam2 == "quasibinomial") {
          //            if (lnk2 == "logit")      testll = f2_binomial_logit(btemp,y,x,mu,P,alpha,wt,0);
          if (lnk2 == "logit")

          {
            testll2 = f2_binomial_logit_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
            //              testll2 = f2_binomial_logit_arma(btemp,y,x,mu,P,alpha,wt,0);
            //              testll2 = f2_binomial_logit(btemp,y,x,mu,P,alpha,wt,0);

          }
          //            else if (lnk2 == "probit") testll = f2_binomial_probit(btemp,y,x,mu,P,alpha,wt,0);
          //            else if (lnk2 == "probit") testll = f2_binomial_probit_arma(btemp,y,x,mu,P,alpha,wt,0);
          else if (lnk2 == "probit")

          {
            //              Rcout << "Enter f2"  << std::endl;

            testll2 = f2_binomial_probit_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
            //                        Rcpp::Rcout << "rmat version: " << testll2  << "\n";

            //              testll2 = f2_binomial_probit_arma(btemp,y,x,mu,P,alpha,wt,0);

            //                          Rcpp::Rcout << "arma version: " << testll2  << "\n";

            //            Rcout << "Exit f2"  << std::endl;
          }

          //            else                       testll = f2_binomial_cloglog(btemp,y,x,mu,P,alpha,wt,0);
          //            else                       testll = f2_binomial_cloglog_arma(btemp,y,x,mu,P,alpha,wt,0);
          else

          {
            testll2 = f2_binomial_cloglog_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);

            //              testll2 = f2_binomial_cloglog_arma(btemp,y,x,mu,P,alpha,wt,0);

          }
        }
        else if (fam2 == "poisson"   || fam2 == "quasipoisson") {

          //            testll = f2_poisson(btemp,y,x,mu,P,alpha,wt,0);
          //            testll = f2_poisson_arma(btemp,y,x,mu,P,alpha,wt,0);

          testll2 = f2_poisson_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);

          //            Rcpp::Rcout << "rmat version v2: " << testll2  << "\n";


          //            testll2 = f2_poisson_rmat(btemp,y,x,mu,P,alpha,wt,0);

          //            Rcpp::Rcout << "rmat version: " << testll2  << "\n";

          //            testll2 = f2_poisson_arma(btemp,y,x,mu,P,alpha,wt,0);

          //            Rcpp::Rcout << "arma version: " << testll2  << "\n";


          //            testll[0]=testll2[0];
        }
        else if (fam2 == "Gamma") {
          //            testll = f2_gamma(btemp,y,x,mu,P,alpha,wt,0);
          //            testll = f2_gamma_arma(btemp,y,x,mu,P,alpha,wt,0);
          testll2 = f2_gamma_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
          //                        Rcpp::Rcout << "rmat version v2: " << testll2  << "\n";
          //            testll2 = f2_gamma_arma(btemp,y,x,mu,P,alpha,wt,0);
          //                        Rcpp::Rcout << "arma version: " << testll2  << "\n";
        }
        else { // gaussian
          //            testll = f2_gaussian(btemp,y,x,mu,P,alpha,wt);
          //            testll = f2_gaussian_arma(btemp,y,x,mu,P,alpha,wt);

          //  Note: This Envelope based sampling method for the Gaussian
          //        is not currently used. May implement future option
          //        to use as this is of theoretica interest
          //        and can be used to validate upper bounds
          testll2 = f2_gaussian_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
          //  Rcpp::Rcout << "rmat version: " << testll2  << "\n";
          //            testll2 = f2_gaussian_arma(btemp,y,x,mu,P,alpha,wt);
          //            Rcpp::Rcout << "arma version: " << testll2  << "\n";

        }

        // calculate and print the acceptance statistic
        //          double test = LLconst[J]+ testtemp2(0,0) - std::log(U2)- testll[0];
        double test = LLconst[J]+ testtemp2(0,0) - std::log(U2)- testll2[0];

        // 5) Accept/reject logic



        if (test >= 0.0) {

          a1 = 1.0;            // accept

        } else {

          // keep existing behavior: increment trial count
          draws[i] = draws[i] + 1.0;

          // effective cap: use max_draws when provided, otherwise use legacy 1000 for diagnostic
          //            int cap = (max_draws >= 0) ? max_draws : 1000;

          // print exactly once when we hit the cap (use your existing mutex for thread-safety)
          //            if (static_cast<int>(draws[i]) == cap) {
          if (max_draws>0 && static_cast<int>(draws[i]) >= max_draws) {
            tbb::mutex::scoped_lock lock(f2_mutex);
            Rcpp::Rcout << "[WARN] index=" << i << " reached draws=" << draws[i]
                        << " (cap=" << max_draws << ") — forcing a1=1.0 to avoid infinite loop\n";

            Rcpp::Rcout << "[DEBUG] Acceptance test breakdown:\n";
            Rcpp::Rcout << "  LLconst[" << J << "] = " << LLconst[J] << "\n";
            Rcpp::Rcout << "  testtemp2(0,0) = " << testtemp2(0,0) << "\n";
            Rcpp::Rcout << "  log(U2) = " << std::log(U2) << "\n";
            Rcpp::Rcout << "  testll2[0] = " << testll2[0] << "\n";
            Rcpp::Rcout << "  test = " << test << "\n";


          }



          // when cap reached or exceeded, set the atomic flag (if provided) and force exit
          //            if (static_cast<int>(draws[i]) >= cap) {
          if (max_draws>0 && static_cast<int>(draws[i]) >= max_draws) {
            if (any_maxdraw_flag) {
              any_maxdraw_flag->store(1, std::memory_order_relaxed);
            }
            a1 = 1.0;   // force acceptance / break out of while loop
          }

        }



      }





    } // while(a1)
  }   // for(i)
  //  Rcpp::Rcout << "Exiting Worker: " << end << std::endl;

}     // operator()






namespace {

// Build pilot return list without Rcpp::Named / CharacterVector::create (rchk-sensitive paths).
Rcpp::List make_rcppparallel_pilot_result(
    SEXP out_sexp,
    SEXP draws_sexp,
    int any_maxdraw_value,
    bool include_message,
    const std::string& message,
    bool include_estimates,
    double est_total,
    double cal_elapsed
) {
  constexpr R_len_t k_n = 6;
  static const char* k_names[k_n] = {
    "out", "draws", "any_maxdraw", "message", "est_total_sec", "cal_elapsed_sec"
  };

  SEXP res = PROTECT(Rf_allocVector(VECSXP, k_n));
  SET_VECTOR_ELT(res, 0, out_sexp);
  SET_VECTOR_ELT(res, 1, draws_sexp);
  SET_VECTOR_ELT(res, 2, Rf_ScalarInteger(any_maxdraw_value));
  if (include_message) {
    SET_VECTOR_ELT(res, 3, Rf_mkString(message.c_str()));
  } else {
    SET_VECTOR_ELT(res, 3, R_NilValue);
  }
  if (include_estimates) {
    SET_VECTOR_ELT(res, 4, Rf_ScalarReal(est_total));
    SET_VECTOR_ELT(res, 5, Rf_ScalarReal(cal_elapsed));
  } else {
    SET_VECTOR_ELT(res, 4, Rf_ScalarReal(NA_REAL));
    SET_VECTOR_ELT(res, 5, Rf_ScalarReal(NA_REAL));
  }

  SEXP nms = PROTECT(Rf_allocVector(STRSXP, k_n));
  for (R_len_t i = 0; i < k_n; ++i) {
    SET_STRING_ELT(nms, i, Rf_mkChar(k_names[i]));
  }
  Rf_setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(1); // nms (now reachable from attributes of res)

  Rcpp::List out(res);
  UNPROTECT(1); // res
  return out;
}

} // namespace


// Keep this helper internal (no // [[Rcpp::export]])
// It uses your existing rNormalGLM_worker objects and the original Rcpp
// out/draws containers directly, avoiding any copying or conversions.
Rcpp::List run_rcppparallel_pilot(
    int n,
    rNormalGLM_worker& test_worker,                  // single-threaded pilot worker
    rNormalGLM_worker& worker,                       // parallel worker
    Rcpp::NumericVector& draws,                      // original Rcpp draws (shared with views)
    Rcpp::NumericMatrix& out,                        // original Rcpp out (shared with views)
    std::shared_ptr<std::atomic<int>> any_flag,      // shared atomic flag set by workers
    double E_draws,                                  // expected candidates per acceptance (scalar)
    bool verbose = true
) {

  // --- single-threaded test run with timing and diagnostic print
  int m_test = 1; // deterministic single-threaded test
  auto t0 = std::chrono::steady_clock::now();
  test_worker(0, m_test);            // invoke worker in serial mode
  auto t1 = std::chrono::steady_clock::now();
  auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  
  // inspect the flag after the test
  int any_hit_after_test = any_flag->load(std::memory_order_relaxed);
  
  if (any_hit_after_test != 0) {
    Rcpp::Rcout
    << "[WARN] One or more indices reached the max_draws cap during the test phase "
    << "with zero accepted draws.\n"
    << "This indicates that the envelope was insufficiently tight overall.\n"
    << "Complete non-acceptance is a strong indicator of posterior non-normality.\n"
    << "\n"
    << "Important note: Once you continue, the full run has no max_draws cap.\n"
    << "Because this code uses RcppParallel, the run cannot be interrupted.\n"
    << "If acceptance remains zero, the simulation may appear to 'hang' indefinitely.\n"
    << "\n"
    << "Recommended actions:\n"
    << "  - Set use_opencl = TRUE or increase the requested sample size (number of draws);\n"
    << "    both of these lead EnvelopeOpt to favor tighter envelopes, though they do not guarantee it.\n"
    << "  - Try a different Gridtype setting to force a tighter envelope.\n"
    << "  - Strengthen the prior to stabilize behavior in the tails.\n"
    << std::endl;
    
    // interactive prompt (use base::interactive(); Rf_interactive is not portable on Windows)
    Rcpp::Function r_interactive("interactive");
    Rcpp::Shield<SEXP> interactive_sexp(r_interactive());
    bool is_interactive = Rcpp::as<bool>(interactive_sexp);
    
    if (is_interactive) {
      Rcpp::Function readline("readline");
      std::string prompt = "Enter 1 to continue full run, 2 to stop and return partial results: ";
      
      while (true) {
        Rcpp::Shield<SEXP> prompt_sexp(Rf_mkString(prompt.c_str()));
        Rcpp::Shield<SEXP> ans_sexp(readline(prompt_sexp));
        std::string ans = Rcpp::as<std::string>(ans_sexp);
        // trim whitespace
        auto ltrim = [](std::string &s) {
          s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                          [](unsigned char ch){ return !std::isspace(ch); }));
        };
        auto rtrim = [](std::string &s) {
          s.erase(std::find_if(s.rbegin(), s.rend(),
                               [](unsigned char ch){ return !std::isspace(ch); }).base(), s.end());
        };
        ltrim(ans); rtrim(ans);
        
        if (ans == "1" || ans == "continue" || ans == "y" || ans == "yes") {
          Rcpp::Rcout << "[INFO] User chose to continue full run.\n";
          break; // fall through to construct/run full worker
        } else if (ans == "2" || ans == "stop" || ans == "n" || ans == "no") {
          Rcpp::Rcout << "[INFO] User chose to stop. Returning partial test results.\n";
          return make_rcppparallel_pilot_result(
            out,
            draws,
            any_hit_after_test,
            true,
            std::string("Stopped by user after test"),
            false,
            NA_REAL,
            NA_REAL
          );
        } else {
          Rcpp::Rcout << "Invalid input. Please enter 1 (continue) or 2 (stop).\n";
        }
      } // end prompt loop
    } else {
      // Non-interactive: proceed automatically
      Rcpp::Rcout << "[NOTE] Non-interactive session: proceeding automatically.\n";
    }
  }
  
  // --- runtime estimate based on single-sample pilot ---
  int candidates_used = static_cast<int>(draws[0]);     // from pilot sample
  double est_total_sec = std::numeric_limits<double>::quiet_NaN();
  
  if (candidates_used > 0 && !Rcpp::NumericVector::is_na(E_draws)) {
    double per_candidate_ms = static_cast<double>(elapsed_ms) / candidates_used;
    double est_total_ms = per_candidate_ms * E_draws * static_cast<double>(n);
    
    // auto fmt_hms = [](double seconds) {
    //   int s = static_cast<int>(std::round(seconds));
    //   int h = s / 3600; s %= 3600;
    //   int m = s / 60;   s %= 60;
    //   std::ostringstream oss;
    //   oss << h << "h " << m << "m " << s << "s";
    //   return oss.str();
    // };
    
    est_total_sec = est_total_ms / 1000.0;
    
    // if(verbose){
    // Rcpp::Rcout << "[rNormalGLM:Pilot] Estimated simulation time (" << n << " draws): "
    //             << est_total_sec << " seconds (" << fmt_hms(est_total_sec) << ").\n"
    //             << "Note: this phase uses RcppParallel and cannot be safely interrupted.\n";
    // }
  
  if (verbose) {
    Rcpp::Rcout
    << "[rNormalGLM:Pilot] Estimated simulation time ("
    << glmbayes::progress::format_int_with_commas(static_cast<long long>(n))
    << " draws): "
    << glmbayes::progress::format_hms(est_total_sec)
    << "\n"
    << "[rNormalGLM:Pilot] Note: this phase uses RcppParallel and cannot be safely interrupted.\n";
  }
  
  }
  
  // --- conservative calibration sizing based on serial bound ---
  double per_candidate_ms_serial = static_cast<double>(elapsed_ms) / std::max(1, candidates_used);
  double est_per_draw_ms_serial  = per_candidate_ms_serial * E_draws; // conservative bound per draw
  // double est_total_ms_serial     = est_per_draw_ms_serial * static_cast<double>(n); // optional
  
  int m1 = std::max(1, (int)std::ceil(0.01 * (double)n));
  int m2 = std::max(1, (int)std::floor(300000.0 / std::max(1.0, est_per_draw_ms_serial))); // 300k ms = ~5 minutes
  int m_stage = std::min(m1, m2);
  
  // if(verbose){
  // Rcpp::Rcout << "Calibrating simulation time estimate using " << m_stage
  //             << " draws at "
  //             << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")()))
  //             << "\n";
  // }
  
  if (verbose) {
    Rcpp::Rcout
    << "[rNormalGLM:Pilot] Calibrating simulation time estimate using "
    << glmbayes::progress::format_int_with_commas(static_cast<long long>(m_stage))
    << " draws at "
    << glmbayes::progress::now_hms()
    << "\n";
  }
  
  // --- calibration run for m_stage draws ---
  auto t_cal0 = std::chrono::steady_clock::now();
  RcppParallel::parallelFor(0, m_stage, worker);
  auto t_cal1 = std::chrono::steady_clock::now();
  double cal_elapsed_sec = std::chrono::duration<double>(t_cal1 - t_cal0).count();
  
  // Sum over all pilot draws (read from original draws container)
  long long total_candidates = 0;
  for (int i = 0; i < m_stage; ++i) {
    total_candidates += static_cast<int>(draws[i]);
  }
  
  // Average candidates per accepted draw (empirical)
  double avg_candidates_per_draw = static_cast<double>(total_candidates) / m_stage;
  
  // Per‑candidate cost from calibration
  double per_candidate_sec = cal_elapsed_sec / std::max(1.0, static_cast<double>(total_candidates));
  
  // Per‑draw cost: scale by expected candidates per acceptance
  double est_per_draw_sec = per_candidate_sec * E_draws;
  
  // Total estimate for n draws
  est_total_sec = est_per_draw_sec * static_cast<double>(n);
  
  // --- print diagnostics ---
  
  // if (verbose){
  // Rcpp::Rcout << "[CALIB] Calibration elapsed = " << cal_elapsed_sec
  //             << " s for " << m_stage << " accepted draws using "
  //             << total_candidates << " candidates.\n";
  // 
  // Rcpp::Rcout << "[CALIB] avg_candidates_per_draw (empirical) = "
  //             << avg_candidates_per_draw
  //             << " vs E_draws = " << E_draws << "\n";
  // 
  // Rcpp::Rcout << "[CALIB] per_candidate_sec = " << per_candidate_sec
  //             << " s, est_per_draw_sec = " << est_per_draw_sec << " s\n";
  // }
  
  if (verbose) {
    Rcpp::Rcout
    << "[rNormalGLM:Pilot] Calibration elapsed = "
    << cal_elapsed_sec << " s for "
    << glmbayes::progress::format_int_with_commas(static_cast<long long>(m_stage))
    << " accepted draws using "
    << glmbayes::progress::format_int_with_commas(static_cast<long long>(total_candidates))
    << " candidates.\n";
    
    Rcpp::Rcout
    << "[rNormalGLM:Pilot] avg_candidates_per_draw (empirical) = "
    << avg_candidates_per_draw
    << " vs E_draws = " << E_draws << "\n";
    
    Rcpp::Rcout
    << "[rNormalGLM:Pilot] per_candidate_sec = "
    << per_candidate_sec
    << " s, est_per_draw_sec = "
    << est_per_draw_sec << " s\n";
  }
  
  
  // auto fmt_hms2 = [](double seconds) {
  //   long long s = static_cast<long long>(std::round(seconds));
  //   long long h = s / 3600; s %= 3600;
  //   long long m = s / 60;   s %= 60;
  //   std::ostringstream oss;
  //   if (h > 0) oss << h << "h ";
  //   if (m > 0 || h > 0) oss << m << "m ";
  //   oss << s << "s";
  //   return oss.str();
  // };
  // 
  // if (verbose) {Rcpp::Rcout << "Refined simulation time estimate (" << n << " draws): "
  //             << est_total_sec << " seconds ("
  //             << fmt_hms2(est_total_sec) << ").\n";
  // }
  
  if (verbose) {
    Rcpp::Rcout
    << "[rNormalGLM:Pilot] Refined simulation time estimate ("
    << glmbayes::progress::format_int_with_commas(static_cast<long long>(n))
    << " draws): "
    << glmbayes::progress::format_hms(static_cast<long>(std::round(est_total_sec)))
    << "\n";
  }
  
  // --- yes/no option if estimate exceeds 5 minutes ---
  if (est_total_sec > 300.0) {
    
    // Blank line before the prompt (matches EnvelopeBuild style)
    Rcpp::Rcout << "\nEstimated simulation exceeds 5 minutes ("
                << glmbayes::progress::format_hms(static_cast<long>(std::round(est_total_sec)))
                << ").\n";
    
    std::string prompt = "Do you want to continue? [y/N]: ";
    
    Rcpp::Function r_interactive("interactive");
    Rcpp::Shield<SEXP> interactive_sexp(r_interactive());
    bool is_interactive = Rcpp::as<bool>(interactive_sexp);
    
    if (is_interactive) {
      
      Rcpp::Function readline("readline");
      
      while (true) {
        Rcpp::Shield<SEXP> prompt_sexp(Rf_mkString(prompt.c_str()));
        Rcpp::Shield<SEXP> ans_sexp(readline(prompt_sexp));
        std::string ans = Rcpp::as<std::string>(ans_sexp);
        
        // trim whitespace
        auto ltrim = [](std::string &s) {
          s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                          [](unsigned char ch){ return !std::isspace(ch); }));
        };
        auto rtrim = [](std::string &s) {
          s.erase(std::find_if(s.rbegin(), s.rend(),
                               [](unsigned char ch){ return !std::isspace(ch); }).base(), s.end());
        };
        ltrim(ans); rtrim(ans);
        
        if (ans == "y" || ans == "yes" || ans == "1" || ans == "continue") {
          
          Rcpp::Rcout << "\nUser chose to continue full run.\n";
          
          Rcpp::Rcout << "[rNormalGLM:Sampler] Running full parallel sampler: "
                      << glmbayes::progress::timestamp_cpp()
                      << "\n";
          
          break;
          
        } else if (ans == "n" || ans == "no" || ans == "2" || ans.empty()) {
          
          Rcpp::Rcout << "User declined. Stopping simulation.\n";
          Rcpp::stop("Simulation stopped by user after time estimate.");
          
        } else {
          Rcpp::Rcout << "Invalid input. Please enter y (continue) or N (stop).\n";
        }
      }
      
    } else {
      
      // Non-interactive session (e.g., CI/CRAN): auto-approve
      Rcpp::Rcout << "[NOTE] Non-interactive session: proceeding automatically.\n";
      Rcpp::Rcout << "Proceeding with full run.\n";
      
      Rcpp::Rcout << "[rNormalGLM:Sampler] Running full parallel sampler: "
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
    }
  }
  
  
  // Return summary with zero-copy references to original out/draws
  return make_rcppparallel_pilot_result(
    out,
    draws,
    any_hit_after_test,
    false,
    std::string(),
    true,
    est_total_sec,
    cal_elapsed_sec
  );
}


List rNormalGLM_std_parallel(
    int                   n,
    NumericVector         y,
    NumericMatrix         x,
    NumericMatrix         mu,
    NumericMatrix         P,
    NumericVector         alpha,
    NumericVector         wt,
    Function              f2,
    List                  Envelope,
    CharacterVector       family,
    CharacterVector       link,
    int                   progbar = 1,
    bool verbose = false
  
) {
  
  // local debug toggle (temporary)
  // const bool debug = true;
  
  // auto mb = [](double bytes) {
  //   return bytes / (1024.0 * 1024.0);
  // };
  
  // 0) allocate output buffers (always a big one)
  int p = mu.nrow();

  
  // allocate output buffers
  NumericMatrix out(n, p);
  NumericVector draws(n);
  

  Rcpp::NumericVector PLSD    = Envelope["PLSD"];
  Rcpp::NumericMatrix loglt   = Envelope["loglt"];
  Rcpp::NumericMatrix logrt   = Envelope["logrt"];
  Rcpp::NumericMatrix cbars   = Envelope["cbars"];
  Rcpp::NumericVector LLconst = Envelope["LLconst"];

  double E_draws = NA_REAL;
  if (Envelope.containsElementNamed("E_draws")) {
    E_draws = Rcpp::as<double>(Envelope["E_draws"]);
  }
  

  // if (verbose)  Rcpp::Rcout << "[rNormalGLM_std] Estimated draws per Acceptance: " << E_draws << "\n";  



  arma::vec PLSD2(PLSD.begin(), PLSD.size(), false, true);
  arma::vec LLconst2(LLconst.begin(), LLconst.size(), false, true);
  arma::mat loglt2(loglt.begin(), loglt.nrow(), loglt.ncol(), false, true);
  arma::mat logrt2(logrt.begin(), logrt.nrow(), logrt.ncol(), false, true);
  arma::mat cbars2(cbars.begin(), cbars.nrow(), cbars.ncol(), false, true);
  
  

  // Create thread-safe views from R-native containers
  RcppParallel::RVector<double> y_r(y);
  RcppParallel::RMatrix<double> x_r(x);
  RcppParallel::RMatrix<double> mu_r(mu);
  RcppParallel::RMatrix<double> P_r(P);
  RcppParallel::RVector<double> alpha_r(alpha);
  RcppParallel::RVector<double> wt_r(wt);
  RcppParallel::RMatrix<double> out_r(out);
  RcppParallel::RVector<double> draws_r(draws);
  
  

  // create shared atomic flag (init 0)
  auto any_flag = std::make_shared<std::atomic<int>>(0);
  
  
  double p_max_draws = 0.001;
  double p_accept = 1.0 / E_draws;
  
  double max_draws = std::ceil(std::log(p_max_draws) / std::log(1.0 - p_accept));
  

  
  RVector<double> PLSD_r(PLSD);
  RVector<double> LLconst_r(LLconst);
  RMatrix<double> loglt_r(loglt);
  RMatrix<double> logrt_r(logrt);
  RMatrix<double> cbars_r(cbars);
  
  


  rNormalGLM_worker test_worker(
      n, y_r, x_r, mu_r, P_r, alpha_r, wt_r,
      PLSD_r, LLconst_r, loglt_r, logrt_r, cbars_r,
      family, link, progbar, out_r, draws_r,
      any_flag, max_draws
  );
  
  

  rNormalGLM_worker worker(
      n,
      y_r, x_r, mu_r, P_r,
      alpha_r, wt_r,
      PLSD_r, LLconst_r, loglt_r, logrt_r, cbars_r,
      family, link,
      progbar,
      out_r, draws_r,
      any_flag,    // shared atomic flag
      -1
  );
  
  

  if (p >= 14) {
   
   auto pilot_res = run_rcppparallel_pilot(
     n,
     test_worker,   // your serial test worker
     worker,        // your parallel worker
     draws,         // original NumericVector
     out,           // original NumericMatrix
     any_flag,      // shared atomic flag
     E_draws,       // scalar double
     verbose
   );
    
  }
  
  

  ////////////////////////////////////////////////////////////////
  
  
  

  
  RcppParallel::parallelFor(0, n, worker);  // grain size == n → serial chunk
//      worker(0, n);  // Call serially


    //  RcppParallel::parallelFor(0, n, worker, n);  // grain size == n → serial chunk
    //int cores = std::thread::hardware_concurrency();
    //int grainSize = std::max(n / cores, 1);
    //parallelFor(0, n, worker, grainSize);
    //parallelFor(0, n, worker);
    

  // return complete out + draws
  return List::create(
    Named("out")   = out,
    Named("draws") = draws
  );
}





///////////////////////////////////////////////////////////////////////////

namespace glmbayes {

namespace sim {



Rcpp::List  rNormalGLM_std(int n,NumericVector y,NumericMatrix x,
                               NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,
                               Function f2,Rcpp::List  Envelope,Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, 
                               int progbar,
                               bool verbose                                  
)
{
  
  
  
  //RNGScope scope;
  int l1 = mu.nrow();
  //  int l2=pow(3,l1);
  
  std::string family2 = Rcpp::as<std::string>(family);
  std::string link2 = Rcpp::as<std::string>(link);  
  
  
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(n,l1);
  double a1=0;
  double a2=0;
  double U=0;
  double test=0;
  double U2=0;
  
  //out(0,0)=1;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
  NumericMatrix cbars=Envelope["cbars"];
  NumericVector LLconst=Envelope["LLconst"]; 
  
  NumericVector outtemp=out(0,_);
  arma::rowvec outtemp2(outtemp.begin(),l1,false);
  NumericVector cbartemp=cbars(0,_);
  arma::rowvec cbartemp2(cbartemp.begin(),l1,false);
  NumericMatrix testtemp(1,1);
  arma::mat testtemp2(testtemp.begin(),1,1,false);
  NumericMatrix btemp(l1,1);
  arma::mat btemp2(btemp.begin(),l1,1,false); 
  NumericVector testll(1);
  
  
  
  
  //Rcout << "PLSD"  << PLSD << std::endl;
  //Rcout << "cbars"  << cbars << std::endl;
  //Rcout << "LLconst"  << LLconst << std::endl;
  
  //  Rcpp::stop("Envelope Components Above");
  
  if(progbar==1){ Rcpp::Rcout << "Starting Simulation:" << std::endl;  };
  for(int i=0;i<n;i++){
    
    Rcpp::checkUserInterrupt();
    if(progbar==1){
      progress_bar(i, n-1);
      if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
    }
    
    
    a1=0;
    draws(i)=1;
    while(a1==0){
      
      U=R::runif(0.0, 1.0);
      a2=0;
      J(i)=0;    
      while(a2==0){
        if(U<=PLSD(J(i))) a2=1;
        if(U>PLSD(J(i))){ 
          U=U-PLSD(J(i));
          J(i)=J(i)+1;
          
        }
        //a2=1; 
      }
      
      
      
      
      for(int j=0;j<l1;j++){  
        
        out(i,j)=rnorm_ct(logrt(J(i),j),loglt(J(i),j),-cbars(J(i),j),1.0);    
        
        
      }
      
      outtemp=out(i,_);
      cbartemp=cbars(J(i),_);
      testtemp2=outtemp2 * trans(cbartemp2);
      U2=R::runif(0.0, 1.0);
      btemp2=trans(outtemp2);    
      
      // Need to modify to call correct f2 function based on family and link
      
      
      
      if(family2=="binomial"){
        if(link2=="logit"){  
          testll=f2_binomial_logit(btemp,y, x,mu,P,alpha,wt,0);
        }
        if(link2=="probit"){  
          testll=f2_binomial_probit(btemp,y, x,mu,P,alpha,wt,0);
        }
        if(link2=="cloglog"){  
          testll=f2_binomial_cloglog(btemp,y, x,mu,P,alpha,wt,0);
        }
      }
      
      if(family2=="quasibinomial"){
        if(link2=="logit"){  
          testll=f2_binomial_logit(btemp,y, x,mu,P,alpha,wt,0);
        }
        if(link2=="probit"){  
          testll=f2_binomial_probit(btemp,y, x,mu,P,alpha,wt,0);
        }
        if(link2=="cloglog"){  
          testll=f2_binomial_cloglog(btemp,y, x,mu,P,alpha,wt,0);
        }
      }
      
      
      if(family2=="poisson"){  
        testll=f2_poisson(btemp,y, x,mu,P,alpha,wt,0);
      }
      if(family2=="quasipoisson"){  
        testll=f2_poisson(btemp,y, x,mu,P,alpha,wt,0);
      }
      
      if(family2=="Gamma"){  
        testll=f2_gamma(btemp,y, x,mu,P,alpha,wt,0);
      }
      
      if(family2=="gaussian"){  
        testll=f2_gaussian(btemp,y, x,mu,P,alpha,wt);
      }
      
      ////
      
      
      test=LLconst(J(i))+testtemp(0,0)-log(U2)-testll(0);
      
      
      if(test>=0) 
      {        a1=1;      }
      if(test<0) draws(i)=draws(i)+1;
      
    }
    
    
  }
  
  //return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws,Rcpp::Named("J")=J,Rcpp::Named("PLSD")=PLSD,Rcpp::Named("famout")=family);
  return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws);
  
}



Rcpp::List rNormalGLM(int n,NumericVector y,NumericMatrix x, 
                          NumericVector mu,NumericMatrix P,NumericVector offset,NumericVector wt,
                          double dispersion,
                            Function f2,
                            Function f3,
                            NumericVector start,
                            std::string family,
                            std::string link,
                            int Gridtype,
                            int n_envopt       ,  // NEW: envelope sizing proxy,
                            bool use_parallel ,       // Enables parallel simulation
                            bool use_opencl ,        // Enables OpenCL acceleration during envelope construction
                            bool verbose 

                            ) {

  if (n_envopt < 0) {
    n_envopt = n; // default fallback
  }
  
  //                          Rcpp::List  famfunc,
  //                            Function f1,
  
  NumericVector offset2=offset;
  
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function asVec("as.vector");
  int l1=x.ncol();
  int l2=x.nrow();
  
  int l1b=mu.length();
  int l1c=P.ncol();
  int l1d=P.nrow();
  
  if(l1b!=l1) Rcpp::stop("Number of rows in mu not consistent with number of columns in matrix x");
  if(l1c!=l1) Rcpp::stop("Number of columns in matrix P not consistent with number of columns in matrix x");
  if(l1d!=l1) Rcpp::stop("Number of rows in matrix P not consistent with number of columns in matrix x");
  
  int l2b=y.length();
  int l2c=offset2.length();
  int l2d=wt.length();
  
  if(l2b!=l2) Rcpp::stop("Number of rows in y not consistent with number of rows in matrix x");
  if(l2c!=l2) Rcpp::stop("Number of rows in offset2 vector not consistent with number of rows in matrix x");
  if(l2d!=l2) Rcpp::stop("Number of rows in wt vector not consistent with number of rows in matrix x");
  
  
  double dispersion2;
  NumericVector alpha(l2);
  NumericMatrix mu2a=asMat(mu);
  NumericMatrix out(l1,n);   
  NumericVector LL(n);   
  
  
  arma::mat x2(x.begin(), l2, l1, false);
  arma::vec alpha2(alpha.begin(),l2,false);  
  arma::vec offset2b(offset2.begin(),l2,false);  
  arma::mat mu2(mu2a.begin(), mu2a.nrow(), mu2a.ncol(), false);
  
  NumericMatrix x2b(clone(x));
  arma::mat P2(P.begin(), P.nrow(), P.ncol(), false);
  
  if(family=="poisson"||family=="binomial")dispersion2=1;
  else dispersion2=dispersion;
  
  int i;  // This can be likely be shifted towards top of function
  
  NumericVector  wt2=wt/dispersion2; // Adjusts weight for dispersion
  arma::vec wt3(wt2.begin(), x.nrow());
  
  // Step 1: Shifts mean and offset to alpha --> Modified offset - Set inputs for optimization
  // Transformed Model now has prior mean 0
  
  alpha2=x2*mu2+offset2b; 
  NumericVector parin=start-mu;  // Starting value for optimization is now start - mu
  NumericVector mu1=mu-mu;       // new prior means are zero
  Rcpp::Function optfun("optim");
  Rcpp::Function tryfun("try");   // add this
  
  
  arma::vec mu1b(mu1.begin(),l2,false);
  
  // Step 2: Run posterior optimization with log-posterior function and gradient functions
  // Note: May eventually replace this with use of a call to modified version of glm.fit
  // Likely would require writing modified family functions that add the prior components
  // This is a bit complex - may make a difference for larger problems or problems
  // where BFGS method for other reasons fails to find True optimmum.

  //NumericVector qc=f2_poisson(parin,y,x,mu1,P,alpha,wt2);

  //Rcout << "Entering optimization" << std::endl;
  
    
  // List opt=optfun(_["par"]=parin,_["fn"]=f2, _["gr"]=f3,_["y"]=y,
  //                 _["x"]=x,
  //                 _["mu"]=mu1,_["P"]=P,_["alpha"]=alpha,_["wt"]=wt2,_["method"]="BFGS",_["hessian"]=true);
  
   ///Trying safe optimization
   
   // Wrap optim() in try()
  Rcpp::Shield<SEXP> optSEXP(tryfun(
     optfun(
       Rcpp::_["par"]     = parin,
       Rcpp::_["fn"]      = f2,
       Rcpp::_["gr"]      = f3,
       Rcpp::_["y"]       = y,
       Rcpp::_["x"]       = x,
       Rcpp::_["mu"]      = mu1,
       Rcpp::_["P"]       = P,
       Rcpp::_["alpha"]   = alpha,
       Rcpp::_["wt"]      = wt2,
       Rcpp::_["method"]  = "BFGS",
       Rcpp::_["hessian"] = true
     ),
    Rcpp::_["silent"] = true
  ));
   
   // Check for try-error using R API
   if (Rf_inherits(optSEXP, "try-error")) {
     Rcpp::stop("Optimization failed in rNormalGLM");
   }
   
   // Safe to convert to List
   Rcpp::List opt(optSEXP);   
  

  // Extract convergence code
  int conv = Rcpp::as<int>(opt["convergence"]);
  
  // If not converged, print the code and message
  if (conv != 0) {
    Rcpp::Rcout << "optim() returned non-zero convergence code: "
                << conv << std::endl;
    
    if (opt.containsElementNamed("message")) {
      Rcpp::Rcout << "optim() message: "
                  << Rcpp::as<std::string>(opt["message"])
                  << std::endl;
    }
  }
  
  //Rcout << "Completed optimization"  << std::endl;
  
  NumericMatrix b2a=asMat(opt[0]);  // optimized value
  NumericVector min1=opt[1]; // Not clear this is used - should be minimum
  int conver1=opt[3]; // check on convergence
  
  // Approximate hessian - Should consider replacing with Hessian based on 
  // known Hessian formula (when available)
  // This could be a source of error 
  
  NumericMatrix A1=opt[5]; 
  
  // Return Error if Optimizaton failed
  
  if(conver1>0){Rcpp::stop("Posterior Optimization failed");}
  
  // Step 3: Standardize the model 
  
//  Rcpp::Rcout << "Standardizing the model:" << std::endl;

//Rcout << "Standardizing Model" <<  std::endl;
  
  
  Rcpp::List Standard_Mod=glmb_Standardize_Model(
    y, 
    x,   // Original design matrix (to be adjusted)
    P,   // Prior Precision Matrix (to be adjusted)
    b2a, // Posterior Mode from optimization (to be adjusted)
    A1  // Precision for Log-Posterior at posterior mode (to be adjusted)
  );

  //Rcout << "Finished Standardizing Model"  << std::endl;
  
  // Get output from call to glmb_Standardize_Model (not sure if they really need to be allocated)     
  // Advantage of allocating may be due to clarity of code in below
  
  NumericVector bstar2_temp=Standard_Mod[0];
  NumericMatrix A_temp=Standard_Mod[1];   //This should be the precision for the Data
  NumericMatrix x2_temp=Standard_Mod[2];
  NumericMatrix mu2_temp=Standard_Mod[3];
  NumericMatrix P2_temp=Standard_Mod[4];   /// This should be the part of prior being shifted to the log-likelihood
  arma::mat L2Inv_temp=Standard_Mod[5];
  arma::mat L3Inv_temp=Standard_Mod[6];
  
  // Calls seem to be different because of setting for gridsort component!
  // When doing samples of just 1, sorting grid is slow when running many samples,
  // using a sorted grid is much faster
  // Most applications where n=1 are likely to be Gibbs samplers
  // There are likely to be few instances where someone runs a small 
  // number of samples greater than 1
  
  // Step 4: Build the Envelope required to simulate from the Standardized Model
  
  //Rcpp::Rcout << "Starting Envelope Creation:" << std::endl;

  Rcpp::List Envelope; // Can move this towards top of the function

  // sortgrid: false when n==1, true when n>1. Sorting orders faces by decreasing
  // probability so simulation checks highest-probability faces first; with many draws
  // only a small number of faces need to be assessed (e.g., not all 47M). The benefit
  // of pre-sorting grows with sample size. For a single draw (e.g., Gibbs step),
  // skipping the sort is typically faster overall. In practice: n==1 for Gibbs, n>>1
  // for iid sampling (default 1000).
  if(n==1){
    Envelope=EnvelopeBuild(bstar2_temp, A_temp,y, x2_temp,mu2_temp,
                             P2_temp,alpha,wt2,family,link,Gridtype, n,n_envopt,false,use_opencl,verbose);
  }
  if(n>1){
    Envelope=EnvelopeBuild(bstar2_temp, A_temp,y, x2_temp,mu2_temp,
                             P2_temp,alpha,wt2,family,link,Gridtype, n,n_envopt,true,use_opencl,verbose);
  }
  
  //  Rcpp::Rcout << "Finished Envelope Creation:" << std::endl;
  
  double E_draws = NA_REAL;
  if (Envelope.containsElementNamed("E_draws")) {
    E_draws = Rcpp::as<double>(Envelope["E_draws"]);
  }

  
  
  if (verbose)  Rcpp::Rcout << "[rNormalGLM] Estimated draws per Acceptance: " << E_draws << "\n";  
  
  
  
  // Step 5: Run the simulation 

  // Rcout << "Starting Simulation"  << std::endl;
  
  int progbar=0;
  
//  Rcpp::List sim=rNormalGLM_std(n,y,x2_temp,mu2_temp,P2_temp,alpha,wt2,
//                                    f2,Envelope,family,link,progbar);
  
  Rcpp::List sim;
  
//  if (n == 1) {
//    sim = rNormalGLM_std(n, y, x2_temp, mu2_temp, P2_temp, alpha, wt2,
//                             f2, Envelope, family, link, progbar);
//  } else {

    
//    sim = rNormalGLM_std_parallel(n, y, x2_temp, mu2_temp, P2_temp, alpha, wt2, f2, Envelope, family, link, progbar);
//      }
  
// Step 5: Run the simulation  


  // Choose serial vs. parallel sampler  
  if (!use_parallel || n == 1) {  
    if (verbose) Rcpp::Rcout << "[rNormalGLM] Running serial sampler (use_parallel=FALSE or n=1): "
                             << glmbayes::progress::timestamp_cpp()
      << "\n";
    
    sim = rNormalGLM_std( n, y, x2_temp, mu2_temp, P2_temp, alpha, wt2,  f2, Envelope, family, link, progbar,verbose);  
  }
  else {  
    if (verbose) Rcpp::Rcout << "[rNormalGLM] Running parallel sampler (use_parallel=TRUE and n>1):"
                             << glmbayes::progress::timestamp_cpp()
                             << "\n";
      
    sim = rNormalGLM_std_parallel(n, y, x2_temp, mu2_temp, P2_temp, alpha, wt2,  f2, Envelope, family, link, progbar,verbose);  
  }  
  

  // if (verbose) Rcpp::Rcout << ">>> Finished Simulation: "
  //                          << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
  //                          << "\n";    
  
  if (verbose) {
    Rcpp::Rcout << "[rNormalGLM:Sampler] Finished simulation: "
                << glmbayes::progress::timestamp_cpp()
                << "\n";
  }    
  

  // Rcout << "Finished Simulation"  << std::endl;
  
  //  Step 6: Undo standardization and do some post processing
  
  //  1) Undo-Standardization of Posterior Precision
  //  2) Undo shifting of prior mean to offset
  //  3) Calculate Log_likelihood (used in model diagnostics)
  
  // These two can likely be shifted towards top of function to make code simpler
  
  NumericMatrix sim2=sim[0];
  
  // Arma matrices pointing to matrices above 
  
  arma::mat sim2b(sim2.begin(), sim2.nrow(), sim2.ncol(), false);
  arma::mat out2(out.begin(), out.nrow(), out.ncol(), false);
  
  out2=L2Inv_temp*L3Inv_temp*trans(sim2b); // reverse transformation
  
  // Add mean back in and compute LL (for post processing)
  // Can add option to not compute LL
  // Note: LL does not seem to be used by downstream functions so can likely be edited out and removed 
  // From output - It is recomputed by summary functions
  
  
  for(i=0;i<n;i++){
    out(_,i)=out(_,i)+mu;  // Add mean vector back 
//    LL[i]=as<double>(f1(_["b"]=out(_,i),_["y"]=y,_["x"]=x,offset2,wt2)); // Calculate log_likelihood
  }
  
  //Rcout << "Leaving *.cpp function"  << std::endl;
  
    
  Rcpp::List Prior=Rcpp::List::create(Rcpp::Named("mean")=mu,Rcpp::Named("Precision")=P);  
  
  
  Rcpp::List outlist=Rcpp::List::create(
    Rcpp::Named("coefficients")=trans(out2),
    Rcpp::Named("coef.mode")=b2a+mu,
    Rcpp::Named("dispersion")=dispersion2,
    Rcpp::Named("Prior")=Prior,
    Rcpp::Named("offset")=offset,
    Rcpp::Named("prior.weights")=wt,
    Rcpp::Named("y")=y,
    Rcpp::Named("x")=x,
    Rcpp::Named("fit")=opt,
    Rcpp::Named("iters")=sim[1],
    Rcpp::Named("Envelope")=Envelope
//  ,  Rcpp::Named("loglike")=LL
  );  
  
  return(outlist);
  
}

}
}

