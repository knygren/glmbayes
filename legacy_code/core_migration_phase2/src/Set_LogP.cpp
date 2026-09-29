// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"

using namespace Rcpp;


// 
// void setlogP_C2(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst){
//   
//   int n = logP.nrow(), k = logP.ncol();
//   int l1 =cbars.ncol();
//   
//   arma::mat logP2(logP.begin(), n, k, false); 
//   NumericVector cbartemp=cbars(0,_);  
//   NumericVector G3temp=G3(0,_);  
//   
//   arma::colvec cbarrow(cbartemp.begin(),l1,false);
//   arma::colvec G3row(G3temp.begin(),l1,false);
//   
//   
//   for(int i=0;i<n;i++){
//     cbartemp=cbars(i,_);  
//     G3temp=G3(i,_);  
//     
//     logP(i,1)=logP(i,0)-NegLL(i)+0.5*arma::as_scalar(cbarrow.t() * cbarrow)+arma::as_scalar(G3row.t() * cbarrow);
//     
//     LLconst(i,0)=NegLL(i)-arma::as_scalar(G3row.t() * cbarrow);
//   }
//   
//   
// }


void setlogP_C2(NumericMatrix logP,
                NumericVector NegLL,
                NumericMatrix cbars,
                NumericMatrix G3,
                NumericMatrix LLconst)
{
  int n  = logP.nrow();
  int k  = logP.ncol();
  int l1 = cbars.ncol();
  
  arma::mat logP2(logP.begin(), n, k, false);
  
  NumericVector cbartemp = cbars(0, _);
  NumericVector G3temp   = G3(0, _);
  
  arma::colvec cbarrow(cbartemp.begin(), l1, false);
  arma::colvec G3row(G3temp.begin(),   l1, false);
  
  // --- Track min/max of logP1 ---
  double minLogP1 = R_PosInf;
  double maxLogP1 = R_NegInf;
  int minIndex = -1;
  int maxIndex = -1;
  
  for (int i = 0; i < n; i++) {
    
    cbartemp = cbars(i, _);
    G3temp   = G3(i, _);
    
    // Armadillo views update automatically
    double quad = 0.5 * arma::as_scalar(cbarrow.t() * cbarrow);
    double lin  = arma::as_scalar(G3row.t() * cbarrow);
    
    logP(i,1) = logP(i,0) - NegLL(i) + quad + lin;
    
    LLconst(i,0) = NegLL(i) - lin;
    
    // Track min/max
    double v = logP(i,1);
    if (v < minLogP1) { minLogP1 = v; minIndex = i; }
    if (v > maxLogP1) { maxLogP1 = v; maxIndex = i; }
  }
  
  // --- Underflow threshold check ---
  double spread = minLogP1 - maxLogP1;   // this will be negative
  
  if (spread < -745) {
    Rcout << "[ERROR] logP1 spread too large: exp(min-max)=0\n"
          << "  minIndex=" << minIndex << "  minLogP1=" << minLogP1 << "\n"
          << "  maxIndex=" << maxIndex << "  maxLogP1=" << maxLogP1 << "\n"
          << "  spread=" << spread << "\n";
    
          {
            NumericVector cb = cbars(minIndex, _);
            NumericVector g3 = G3(minIndex, _);
            arma::colvec c(cb.begin(), l1, false);
            arma::colvec g(g3.begin(), l1, false);
            
            double quad = 0.5 * arma::as_scalar(c.t() * c);
            double lin  = arma::as_scalar(g.t() * c);
            
            Rcout << "[DETAIL] Components at minIndex=" << minIndex << "\n"
                  << "  logP0=" << logP(minIndex,0) << "\n"
                  << "  NegLL=" << NegLL[minIndex] << "\n"
                  << "  quad="  << quad << "\n"
                  << "  lin="   << lin  << "\n"
                  << "  logP1=" << logP(minIndex,1) << "\n";
            
            // NEW: print cbars and its norm
            Rcout << "  cbar[minIndex] = " << cb << "\n";
            Rcout << "  ||cbar[minIndex]|| = " << arma::norm(c) << "\n";
          }
          
          
          {
            NumericVector cb = cbars(maxIndex, _);
            NumericVector g3 = G3(maxIndex, _);
            arma::colvec c(cb.begin(), l1, false);
            arma::colvec g(g3.begin(), l1, false);
            
            double quad = 0.5 * arma::as_scalar(c.t() * c);
            double lin  = arma::as_scalar(g.t() * c);
            
            Rcout << "[DETAIL] Components at maxIndex=" << maxIndex << "\n"
                  << "  logP0=" << logP(maxIndex,0) << "\n"
                  << "  NegLL=" << NegLL[maxIndex] << "\n"
                  << "  quad="  << quad << "\n"
                  << "  lin="   << lin  << "\n"
                  << "  logP1=" << logP(maxIndex,1) << "\n";
            
            // NEW: print cbars and its norm
            Rcout << "  cbar[maxIndex] = " << cb << "\n";
            Rcout << "  ||cbar[maxIndex]|| = " << arma::norm(c) << "\n";
          }
          
  //        stop("[EnvelopeBuild] logP1 underflow risk detected.");
  }
}



namespace glmbayes{

namespace env {

Rcpp::List   EnvelopeSet_LogP(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3) {
  
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

}
}