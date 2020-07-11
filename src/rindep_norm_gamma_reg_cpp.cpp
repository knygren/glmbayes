// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include "Set_Grid.h"
#include <math.h>


using namespace Rcpp;



// [[Rcpp::export(".rindep_norm_gamma_reg_std_cpp")]]

Rcpp::List  rindep_norm_gamma_reg_std_cpp(int n,NumericVector y,NumericMatrix x,
                                          NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,
                                          Function f2,Rcpp::List  Envelope,Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, int progbar=1)
{
  RNGScope scope;
  int l1 = mu.nrow();
  //  int l2=pow(3,l1);
  
  std::string family2 = Rcpp::as<std::string>(family);
  std::string link2 = Rcpp::as<std::string>(link);  
  
  
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(n,l1);
  //double a1=0;
  double a2=0;
  double U=0;
  double test=0;
//  double test_int=0;
//  double test_data=0;
  double U2=0;
  
  //out(0,0)=1;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
  NumericMatrix cbars=Envelope["cbars"];
  

//  NumericMatrix cbars_int=Envelope["cbars_int"];
  NumericVector LLconst=Envelope["LLconst"]; 
//  NumericVector LLconst_int=Envelope["LLconst_int"]; 

  NumericVector outtemp=out(0,_);
  arma::rowvec outtemp2(outtemp.begin(),l1,false);
  NumericVector cbartemp=cbars(0,_);
  arma::rowvec cbartemp2(cbartemp.begin(),l1,false);

  // Add Corresponding prior terms here
  
//  NumericVector cbartemp_int=cbars_int(0,_);
//  arma::rowvec cbartemp_int2(cbartemp_int.begin(),l1,false);
  

    NumericMatrix testtemp(1,1);
  arma::mat testtemp2(testtemp.begin(),1,1,false);

  NumericMatrix testtemp_int(1,1);
  arma::mat testtemp_int2(testtemp_int.begin(),1,1,false);
  
  
    NumericMatrix btemp(l1,1);
  arma::mat btemp2(btemp.begin(),l1,1,false); 
  NumericVector testll(1);
  NumericVector testll_data(1);


    
  if(progbar==1){ Rcpp::Rcout << "Starting Simulation:" << std::endl;  };

    Rcpp::checkUserInterrupt();
    if(progbar==1){
//      progress_bar2(i, n-1);
//      if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
    }
    

    //a1=0;

      U=R::runif(0.0, 1.0);
      a2=0;
      J(0)=0;    
      while(a2==0){
        if(U<=PLSD(J(0))) a2=1;
        if(U>PLSD(J(0))){ 
          U=U-PLSD(J(0));
          J(0)=J(0)+1;
          
        }
        //a2=1; 
      }
      

      for(int j=0;j<l1;j++){  
        
        // Switch to using thetabars here
        
        out(0,j)=ctrnorm_cpp(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);    
        
        
      }


            
      // cbars_int2 holds the intercept part
      
      outtemp=out(0,_);   // Hopefully this doew not break link to outtemp2
      cbartemp=cbars(J(0),_);
      testtemp2=outtemp2 * trans(cbartemp2);
//      testtemp_int2=outtemp2 * trans(cbartemp_int2);
      U2=R::runif(0.0, 1.0);
      btemp2=trans(outtemp2);    
      
      // Need to modify to call correct f2 function based on family and link
      
      
      if(family2=="Gamma"){  
        testll=f2_gamma(btemp,y, x,mu,P,alpha,wt,0);
      }
      
      if(family2=="gaussian"){  
        testll=f2_gaussian(btemp,y, x,mu,P,alpha,wt);
        testll_data=f1_gaussian(btemp,y, x,alpha,wt);
      }
      

      test=LLconst(J(0))+testtemp(0,0)-testll(0);


//      test_int=LLconst_int(J(0))+testtemp_int(0,0)-(testll(0)-testll_data(0));
//      test_data=test-test_int;


  return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws,Rcpp::Named("test")=test,Rcpp::Named("J")=J,
                                        Rcpp::Named("log_U2")=log(U2)
//                                      ,  Rcpp::Named("test_int")=test_int,
//                                        Rcpp::Named("test_data")=test_data
                              );      
  
}

