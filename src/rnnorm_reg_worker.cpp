#include "rnnorm_reg_worker.h"

// Only add more includes if strictly needed
#include <cmath>         // for std::log or std::exp if used
#include "famfuncs.h"
#include "Set_Grid.h"


// mutex to protect Rcpp calls
static std::mutex f2_mutex;


  // operator() implements the parallel loop
  void rnnorm_reg_worker::operator()(std::size_t begin, std::size_t end) {
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


