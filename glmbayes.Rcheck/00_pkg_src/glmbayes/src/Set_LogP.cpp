// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"

using namespace Rcpp;



void setlogP_C2(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst){
  
  int n = logP.nrow(), k = logP.ncol();
  int l1 =cbars.ncol();
  
  arma::mat logP2(logP.begin(), n, k, false); 
  NumericVector cbartemp=cbars(0,_);  
  NumericVector G3temp=G3(0,_);  
  
  arma::colvec cbarrow(cbartemp.begin(),l1,false);
  arma::colvec G3row(G3temp.begin(),l1,false);
  
  
  for(int i=0;i<n;i++){
    cbartemp=cbars(i,_);  
    G3temp=G3(i,_);  
    
    logP(i,1)=logP(i,0)-NegLL(i)+0.5*arma::as_scalar(cbarrow.t() * cbarrow)+arma::as_scalar(G3row.t() * cbarrow);
    
    LLconst(i,0)=NegLL(i)-arma::as_scalar(G3row.t() * cbarrow);
  }
  
  
}


// [[Rcpp::export(".setlogP_cpp")]]


Rcpp::List   setlogP(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3) {
  
  int n = logP.nrow(), k = logP.ncol();
  int l1 =cbars.ncol();
  //    int l2=cbars.nrow();
  
  arma::mat logP2(logP.begin(), n, k, false); 
  NumericVector cbartemp=cbars(0,_);  
  NumericVector G3temp=G3(0,_);  
  Rcpp::NumericMatrix LLconst(n,1);
  
  arma::colvec cbarrow(cbartemp.begin(),l1,false);
  arma::colvec G3row(G3temp.begin(),l1,false);
  
  //    double v = arma::as_scalar(cbarrow.t() * cbarrow);
  //    LLconst[j]<--t(as.matrix(cbars[j,1:l1]))%*%t(as.matrix(G3[j,1:l1]))+NegLL[j]    
  
  for(int i=0;i<n;i++){
    cbartemp=cbars(i,_);  
    G3temp=G3(i,_);  
    logP(i,1)=logP(i,0)-NegLL(i)+0.5*arma::as_scalar(cbarrow.t() * cbarrow)+arma::as_scalar(G3row.t() * cbarrow);
    LLconst(i,0)=NegLL(i)-arma::as_scalar(G3row.t() * cbarrow);
  }
  
  
  //    return logP;
  return Rcpp::List::create(Rcpp::Named("logP")=logP,Rcpp::Named("LLconst")=LLconst);
  
}




//////////////////////////////////////////////////////////////////////////////

Rcpp::List   setlogP_C(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst) {
  
  int n = logP.nrow(), k = logP.ncol();
  int l1 =cbars.ncol();
  
  arma::mat logP2(logP.begin(), n, k, false); 
  NumericVector cbartemp=cbars(0,_);  
  NumericVector G3temp=G3(0,_);  
  
  arma::colvec cbarrow(cbartemp.begin(),l1,false);
  arma::colvec G3row(G3temp.begin(),l1,false);
  
  
  for(int i=0;i<n;i++){
    cbartemp=cbars(i,_);  
    G3temp=G3(i,_);  

    // Remark 6 in Nygren and Nygren (2006)
    // logP is log_density for component
    // -NegLL (is g())
    // last term is log of denominator 
    // 3rd term is MGF from Claim1
    
    logP(i,1)=logP(i,0)-NegLL(i)+0.5*arma::as_scalar(cbarrow.t() * cbarrow)+arma::as_scalar(G3row.t() * cbarrow);
    
    LLconst(i,0)=NegLL(i)-arma::as_scalar(G3row.t() * cbarrow);
  }
  
  
  return Rcpp::List::create(Rcpp::Named("logP")=logP,Rcpp::Named("LLconst")=LLconst);
  
}



