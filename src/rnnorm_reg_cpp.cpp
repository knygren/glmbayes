// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppParallel)]]
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "RcppParallel.h"

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include "Set_Grid.h"
#include <math.h>
#include <tbb/mutex.h>
#include <thread>
#include "rng_utils.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;


// Outside Function declarations (move or delete)

void  f4_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, NumericVector NegLL, NumericMatrix cbars, int progbar=0);


void progress_bar2(double x, double N)
{
  // how wide you want the progress meter to be
  int totaldotz=40;
  double fraction = x / N;
  // part of the progressmeter that's already "full"
  int dotz = round(fraction * totaldotz);
  
  Rcpp::Rcout.precision(3);
  Rcout << "\r                                                                 " << std::flush ;
  Rcout << "\r" << std::flush ;
  Rcout << std::fixed << fraction*100 << std::flush ;
  Rcout << "% [" << std::flush ;
  int ii=0;
  for ( ; ii < dotz;ii++) {
    Rcout << "=" << std::flush ;
  }
  // remaining part (spaces)
  for ( ; ii < totaldotz;ii++) {
    Rcout << " " << std::flush ;
  }
  // and back to line begin 
  
  Rcout << "]" << std::flush ;
  
  // and back to line begin 
  
  Rcout << "\r" << std::flush ;
  
}

// [[Rcpp::export(".glmb_Standardize_Model_cpp")]]


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


//-----------------------------------------------------------------------------
// Forward-declare your serial sampler (no change)
//-----------------------------------------------------------------------------
Rcpp::List rnnorm_reg_std_cpp(
    int                   n,
    Rcpp::NumericVector   y,
    Rcpp::NumericMatrix   x,
    Rcpp::NumericMatrix   mu,
    Rcpp::NumericMatrix   P,
    Rcpp::NumericVector   alpha,
    Rcpp::NumericVector   wt,
    Rcpp::Function        f2,
    Rcpp::List            Envelope,
    Rcpp::CharacterVector family,
    Rcpp::CharacterVector link,
    int                   progbar 
);



// mutex to protect Rcpp calls
static std::mutex f2_mutex;


//static tbb::mutex rng_mutex;

//inline double safe_runif() {
//  tbb::mutex::scoped_lock lock(rng_mutex);
//  return R::runif(0.0, 1.0);
//}



//-----------------------------------------------------------------------------
// rnnorm_reg_worker: full sampler worker with debug print of `test`
//-----------------------------------------------------------------------------
struct rnnorm_reg_worker : public Worker {
  // inputs
  int                   n;
  // NumericVector         y;
  // NumericMatrix         x;
  // NumericMatrix         mu;
  // NumericMatrix         P;
  // NumericVector         alpha;
  // NumericVector         wt;
  
  RcppParallel::RVector<double> y_r;
  RcppParallel::RMatrix<double> x_r;
  RcppParallel::RMatrix<double> mu_r;
  RcppParallel::RMatrix<double> P_r;
  RcppParallel::RVector<double> alpha_r;
  RcppParallel::RVector<double> wt_r;
  
//  Function              f2;
  //List                  Envelope;
  // NEW: envelope components in arma format
  arma::vec PLSD, LLconst;
  arma::mat loglt, logrt, cbars;
  
    CharacterVector       family;
  CharacterVector       link;
  int                   progbar;
  
  
  
  // outputs
  RMatrix<double>       out;
  RVector<double>       draws;
  int                   ncol;
  
  // constructor
  rnnorm_reg_worker(
    int                   n_,
    // const NumericVector&  y_,
    // const NumericMatrix&  x_,
    // const NumericMatrix&  mu_,
    // const NumericMatrix&  P_,
    // const NumericVector&  alpha_,
    // const NumericVector&  wt_,

    const RcppParallel::RVector<double>& y_r_,
    const RcppParallel::RMatrix<double>& x_r_,
    const RcppParallel::RMatrix<double>& mu_r_,
    const RcppParallel::RMatrix<double>& P_r_,
    const RcppParallel::RVector<double>& alpha_r_,
    const RcppParallel::RVector<double>& wt_r_,
    
//    const Rcpp::List& Envelope_,
    const arma::vec& PLSD_,
    const arma::vec& LLconst_,
    const arma::mat& loglt_,
    const arma::mat& logrt_,
    const arma::mat& cbars_,
    
//    const Function&       f2_,
    const CharacterVector& family_,
    const CharacterVector& link_,
    int                   progbar_,
    RcppParallel::RMatrix<double>& out_,
    RcppParallel::RVector<double>& draws_
  
  // RcppParallel::RMatrix<double>& out_,
  // RcppParallel::RVector<double>& draws_
  )
//     : n(n_), y(y_), x(x_), mu(mu_), P(P_),
//       alpha(alpha_), wt(wt_),
// //      Envelope(Envelope_),
//       PLSD(PLSD_),  LLconst(LLconst_),
//       loglt(loglt_), logrt(logrt_), cbars(cbars_),
//       //f2(f2_),
//        family(family_), link(link_),
//       progbar(progbar_),
//       out(out_), draws(draws_),
//       ncol(out_.ncol())
    : n(n_),
      y_r(y_r_), x_r(x_r_), mu_r(mu_r_), P_r(P_r_),
      alpha_r(alpha_r_), wt_r(wt_r_),
      PLSD(PLSD_), LLconst(LLconst_),
      loglt(loglt_), logrt(logrt_), cbars(cbars_),
      family(family_), link(link_), progbar(progbar_),
      out(out_), draws(draws_), ncol(out_.ncol())

    
  {}
  


  // operator() implements the parallel loop
  void operator()(std::size_t begin, std::size_t end) {
//    Rcpp::RNGScope scope;  // enable RNG in threads


    // Wrap R-native inputs into thread-safe views
    // RcppParallel::RVector<double> y_r(y);           // observed counts
    // RcppParallel::RMatrix<double> x_r(x);           // design matrix
    // RcppParallel::RMatrix<double> mu_r(mu);         // mode vector
    // RcppParallel::RMatrix<double> P_r(P);           // precision matrix
    // RcppParallel::RVector<double> alpha_r(alpha);   // predictor offset
    // RcppParallel::RVector<double> wt_r(wt);         // observation weights
    

        
    // Convert NumericMatrix and NumericVector inputs to Armadillo
//    arma::vec y2(y.begin(), y.size(), false);
//    arma::vec alpha2(alpha.begin(), alpha.size(), false);
//    arma::vec wt2(wt.begin(), wt.size(), false);
    
//    arma::mat x2(x.begin(), x.nrow(), x.ncol(), false);
//    arma::mat mu2(mu.begin(), mu.nrow(), mu.ncol(), false);
//    arma::mat P2(P.begin(), P.nrow(), P.ncol(), false);        

        
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

    
    
    //NumericVector cbartemp=cbars(0,_);
    //arma::rowvec cbartemp2(cbartemp.begin(),l1,false);
    

            
//    NumericMatrix       btemp(l1,1);
//    arma::mat           btemp2(btemp.begin(),     l1,1,false);
    
//    Rcpp::NumericMatrix btemp(l1,1);                     // stays for interface
//    RcppParallel::RMatrix<double> btemp_r(btemp);        // thread-safe wrapper
//    arma::mat btemp2(btemp_r.begin(), l1, 1, false);     // links directly to RMatrix view
    
    std::vector<double> btemp_buf(l1);
    arma::mat btemp2(btemp_buf.data(), l1, 1, false);
    RcppParallel::RMatrix<double> btemp_r(btemp_buf.data(), l1, 1); // optional: only if still needed
    
    
    arma::mat testtemp2(1, 1);  // Allocated directly on the heap
//    NumericVector       testll(1);
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
        double U = safe_runif();
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
          out(i, j) = ctrnorm_cpp(logrt(J, j),loglt(J, j),-cbars(J, j), 1.0 );
        //  out(i, j) = ctrnorm_cpp(logrt2(J, j),loglt2(J, j),-cbars2(J, j), 1.0 );
        }
        
        // 3) prepare for test
        for (int j = 0; j < l1; ++j) {
          outtemp_buf[j]  = out(i, j);
          cbartemp_buf[j] = cbars(J, j);
          //cbartemp_buf[j] = cbars2(J, j);
          
          
        }
        testtemp2 = outtemp2 * trans(cbartemp2);
  //      double U2 = R::runif(0.0, 1.0);

        double U2 = safe_runif();
        
          
        btemp2   = trans(outtemp2);
        
        // declare test here so it’s in scope below
        //double test;
        
        
        
        // 4) compute log‐lik and print test under lock
        {
          std::lock_guard<std::mutex> guard(f2_mutex);

          
                    
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
              
              draws[i]=draws[i]+1.0;     // reject and try again
              
                    }
            

        }
        
       
                  
         

             } // while(a1)
    }   // for(i)
  //  Rcpp::Rcout << "Exiting Worker: " << end << std::endl;
    
  }     // operator()
};







// [[Rcpp::export(".rnnorm_reg_std_cpp")]]

Rcpp::List  rnnorm_reg_std_cpp(int n,NumericVector y,NumericMatrix x,
NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,
Function f2,Rcpp::List  Envelope,Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, int progbar=1)
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
        progress_bar2(i, n-1);
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
          
          out(i,j)=ctrnorm_cpp(logrt(J(i),j),loglt(J(i),j),-cbars(J(i),j),1.0);    

          
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
      
      test=LLconst(J(i))+testtemp(0,0)-log(U2)-testll(0);

      if(test>=0) a1=1;
	    if(test<0) draws(i)=draws(i)+1;
      
      }
    
    
    }

//return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws,Rcpp::Named("J")=J,Rcpp::Named("PLSD")=PLSD,Rcpp::Named("famout")=family);
return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws);

}





//-----------------------------------------------------------------------------
// test_all_args: wrapper matching rnnorm_reg_std_cpp signature
//-----------------------------------------------------------------------------
// [[Rcpp::export("test_all_args")]]
List rnnorm_reg_std_cpp_parallel(
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
    int                   progbar = 1
) {
  // allocate output buffers
  int p = mu.nrow();
  NumericMatrix out(n, p);
  NumericVector draws(n);
  
  // Extract from Envelope List

  Rcpp::NumericVector PLSD    = Envelope["PLSD"];
  Rcpp::NumericMatrix loglt   = Envelope["loglt"];
  Rcpp::NumericMatrix logrt   = Envelope["logrt"];
  Rcpp::NumericMatrix cbars   = Envelope["cbars"];
  Rcpp::NumericVector LLconst = Envelope["LLconst"];
  
  // Convert to Armadillo (deep-safe versions)
  arma::vec PLSD2    = Rcpp::as<arma::vec>(PLSD);
  arma::vec LLconst2 = Rcpp::as<arma::vec>(LLconst);
  arma::mat loglt2   = Rcpp::as<arma::mat>(loglt);
  arma::mat logrt2   = Rcpp::as<arma::mat>(logrt);
  arma::mat cbars2   = Rcpp::as<arma::mat>(cbars);
  
  
//  Rcpp::Rcout << " 1.0 Launching Worker \n" <<;
//  Rcpp::Rcout << "1.0 Launching Worker " <<  std::endl;
  
  // Create thread-safe views from R-native containers
  RcppParallel::RVector<double> y_r(y);
  RcppParallel::RMatrix<double> x_r(x);
  RcppParallel::RMatrix<double> mu_r(mu);
  RcppParallel::RMatrix<double> P_r(P);
  RcppParallel::RVector<double> alpha_r(alpha);
  RcppParallel::RVector<double> wt_r(wt);
  RcppParallel::RMatrix<double> out_r(out);
  RcppParallel::RVector<double> draws_r(draws);
  
  // launch parallel worker
  
  // rnnorm_reg_worker worker(
  //     n, y, x, mu, P,
  //     alpha, wt,
  //   //  Envelope,
  //     PLSD2,LLconst2, loglt2, logrt2, cbars2, 
  //     family, link,
  //     progbar,
  //     out, draws
  // );
  
  rnnorm_reg_worker worker(
      n,
      y_r, x_r, mu_r, P_r,
      alpha_r, wt_r,
      PLSD2, LLconst2, loglt2, logrt2, cbars2, 
      family, link,
      progbar,
      out_r, draws_r
  );
  
     /// Calling the workers (Paralle or - for testing - serially)
    
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


// [[Rcpp::export(".rnnorm_reg_cpp")]]
Rcpp::List rnnorm_reg_cpp(int n,NumericVector y,NumericMatrix x, 
                          NumericVector mu,NumericMatrix P,NumericVector offset,NumericVector wt,
                          double dispersion,
                            Function f2,Function f3,NumericVector start,
                            std::string family="binomial",
                            std::string link="logit",
                            int Gridtype=2      
) {

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
  
  arma::vec mu1b(mu1.begin(),l2,false);
  
  // Step 2: Run posterior optimization with log-posterior function and gradient functions
  // Note: May eventually replace this with use of a call to modified version of glm.fit
  // Likely would require writing modified family functions that add the prior components
  // This is a bit complex - may make a difference for larger problems or problems
  // where BFGS method for other reasons fails to find True optimmum.

  //NumericVector qc=f2_poisson(parin,y,x,mu1,P,alpha,wt2);

  //Rcout << "Entering optimization" << std::endl;
  
    
  List opt=optfun(_["par"]=parin,_["fn"]=f2, _["gr"]=f3,_["y"]=y,
                  _["x"]=x,
                  _["mu"]=mu1,_["P"]=P,_["alpha"]=alpha,_["wt"]=wt2,_["method"]="BFGS",_["hessian"]=true);
  

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
  NumericMatrix A_temp=Standard_Mod[1];
  NumericMatrix x2_temp=Standard_Mod[2];
  NumericMatrix mu2_temp=Standard_Mod[3];
  NumericMatrix P2_temp=Standard_Mod[4];
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
  
  
  if(n==1){
    Envelope=EnvelopeBuild_c(bstar2_temp, A_temp,y, x2_temp,mu2_temp,
                             P2_temp,alpha,wt2,family,link,Gridtype, n,false);


    
      }
  
  if(n>1){
    Envelope=EnvelopeBuild_c(bstar2_temp, A_temp,y, x2_temp,mu2_temp,
                             P2_temp,alpha,wt2,family,link,Gridtype, n,true);
  }
  
  //  Rcpp::Rcout << "Finished Envelope Creation:" << std::endl;
  
  // Step 5: Run the simulation 

  // Rcout << "Starting Simulation"  << std::endl;
  
  int progbar=0;
  
//  Rcpp::List sim=rnnorm_reg_std_cpp(n,y,x2_temp,mu2_temp,P2_temp,alpha,wt2,
//                                    f2,Envelope,family,link,progbar);
  
  Rcpp::List sim;
  
  if (n == 1) {
    sim = rnnorm_reg_std_cpp(n, y, x2_temp, mu2_temp, P2_temp, alpha, wt2,
                             f2, Envelope, family, link, progbar);
  } else {
//    sim = rnnorm_reg_std_parallel(n, y, x2_temp, mu2_temp, P2_temp, alpha, wt2,
//                                  f2, Envelope, family, link, progbar);

   // sim = test_parallel(n, y, x2_temp, mu2_temp, P2_temp, alpha, wt2,
  //                                f2, Envelope, family, link, progbar);
    
    
    sim = rnnorm_reg_std_cpp_parallel(n, y, x2_temp, mu2_temp, P2_temp, alpha, wt2, f2, Envelope, family, link, progbar);
  //  test_all_args
    //sim=test_parallel( n, l1);
      }
  
  

  // Rcout << "Finished Simulation"  << std::endl;
  
  //  Step 6: Undo standaridzation and do some post processing
  
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
  // Note: LL does not seem to be used by downstream functins so can likely be edited out and removed 
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



















