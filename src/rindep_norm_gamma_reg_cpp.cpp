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
  
  if(progbar==1){ Rcpp::Rcout << "Starting Simulation:" << std::endl;  };
  for(int i=0;i<n;i++){
    
    Rcpp::checkUserInterrupt();
    if(progbar==1){
//      progress_bar2(i, n-1);
//      if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
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

      
      return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws,Rcpp::Named("test")=test);      
      
      if(test>=0) a1=1;
      if(test<0) draws(i)=draws(i)+1;
      
    }
    
    
  }
  
  //return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws,Rcpp::Named("J")=J,Rcpp::Named("PLSD")=PLSD,Rcpp::Named("famout")=family);
  return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws);
  
}

