// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <math.h>


using namespace Rcpp;

double ctrnorm_cpp(double lgrt,double lglt,double mu,double sigma){
  RNGScope scope;
  
  double U=0;
  double out=0;
  double lgU2=0;
  
  if(lgrt>=lglt){
    U=R::runif(0.0, 1.0);
    
    double  u1=1-exp(lgrt);
    double lgu1=log(u1);
    lgU2=log(U)+lglt+log(1-exp(lgu1-lglt));
    double lgU3=lgU2+log(1+exp(lgu1-lgU2));
    out=R::qnorm(lgU3,mu,sigma,TRUE,TRUE);
    
    
  }  
  
  if(lgrt<lglt){
    U=R::runif(0.0, 1.0);
    double e1mu2=1-exp(lglt);
    double lg1mu2=log(e1mu2);
    double lgU2=log(U)+lgrt+log(1-exp(lg1mu2-lgrt));
    double lgU3=lgU2+log(1+exp(lg1mu2-lgU2));
    out=R::qnorm(lgU3,mu,sigma,FALSE,TRUE);
    
  }
  
  
  return out;
  
}
