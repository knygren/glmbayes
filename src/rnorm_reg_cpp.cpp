// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export(".rnorm_reg_cpp")]]

Rcpp::List rnorm_reg_cpp(int n,NumericVector y,NumericMatrix x, 
                         NumericVector mu,NumericMatrix P,
                         NumericVector offset2,NumericVector wt,
                         double dispersion,
                         Rcpp::List famfunc, 
                         Function f1,Function f2,Function f3,
                           NumericVector start,
                           std::string family="binomial",
                           std::string link="logit",
                           int Gridtype=2      
) {
  
  // add checks to make sure that dimensions are consistent
  // (i) number of rows in y, offset2, and wt2 should equal rows of x
  // (ii) number of rows in mu should equal number of columns in x
  // (iii) number of rows and columns in P should equal number of columns in x
  
  
  // Need to check combination of weighting and offset working properly
  
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function asVec("as.vector");
  int l1=x.ncol();
  int l2=x.nrow();
  
  int l1b=mu.length();
  int l1c=P.ncol();
  int l1d=P.nrow();
  
  if(l1b!=l1) Rcpp::stop("Number of rows in mu not consistent with number of columns in matrix x");;
  if(l1c!=l1) Rcpp::stop("Number of columns in matrix P not consistent with number of columns in matrix x");;
  if(l1d!=l1) Rcpp::stop("Number of rows in matrix P not consistent with number of columns in matrix x");;
  
  int l2b=y.length();
  int l2c=offset2.length();
  int l2d=wt.length();
  
  if(l2b!=l2) Rcpp::stop("Number of rows in y not consistent with number of rows in matrix x");;
  if(l2c!=l2) Rcpp::stop("Number of rows in offset2 vector not consistent with number of rows in matrix x");;
  if(l2d!=l2) Rcpp::stop("Number of rows in wt vector not consistent with number of rows in matrix x");;
  
  double dispersion2=dispersion;
  //    NumericVector alpha(l2);
  NumericMatrix mu2a=asMat(mu);
  
  arma::mat x2(x.begin(), l2, l1, false);
  //    arma::vec alpha2(alpha.begin(),l2,false);  
  arma::vec offset2b(offset2.begin(),l2,false);  
  arma::mat mu2(mu2a.begin(), mu2a.nrow(), mu2a.ncol(), false);
  
  NumericMatrix x2b(clone(x));
  arma::mat x2bb(x2b.begin(), l2, l1, false);
  arma::mat P2(P.begin(), P.nrow(), P.ncol(), false);
  
  // Adjusts for dispersion here 
  
  NumericVector  wt2=wt/dispersion2;
  arma::vec wt3(wt2.begin(), x.nrow());
  
  // Is this needed (Check Offset code works properly-Both new and old code)
  
  
  // Move inside If Statement once confirmed works
  
  int i;
  
  // Should this subtract alpha2 [subtraction of offset likely makes this ok]?
  
  NumericVector  y1=y-offset2;
  arma::vec y2b(y1.begin(),l2,false);
  NumericMatrix W1(l2+l1,l1);
  arma::mat W(W1.begin(), W1.nrow(), W1.ncol(), false);
  NumericVector z1(l2+l1);
  arma::vec z(z1.begin(),l2+l1,false);
  
  for(i=0;i<l2;i++){
    x2b(i,_) =x2b(i,_)*sqrt(wt2[i]);
    y1(i)=y1(i)*sqrt(wt2[i]);
  }        
  
  
  arma::mat RA=arma::chol(P2);
  
  // Stacks rows of design matrix and RA (Cholesky of Prior)  
  
  W.rows(0,l2-1)=x2bb;
  W.rows(l2,l2+l1-1)=RA;
  
  // Stacks data and Cholesky multiplied by Cholesky (think Square root of P2)
  
  z.rows(0,l2-1)=y2b;
  z.rows(l2,l1+l2-1)=RA*mu2;

  // The function (.lm.fit should be used here is it is available)
  
   // lm.fit is R Wrappper around *.cpp function named C_Cdqrls
   // Not sure how to call C_cdqrls directly.
   // For now, set up to call lm.fit [Figure out later]
   // Likely safer to call lm.fit as it has wrappers to catch problems (although slower)
   
   Rcpp::Function lm_fit_fun("lm.fit");
   
   // Call lm_fit_fun to get the posterior mode
   // Should be numerically more accurate for large matrices that old method

   List fit=lm_fit_fun(_["x"]=W,_["y"]=z);
   
   NumericMatrix b2a=asMat(fit[0]);  // optimized value (coefficients)
   // Rcpp::Rcout << "Optimized value from lm_fit :"  << b2a << std::endl;
   
  // This should be WTW= XTX +P
  // Seems redundant to store here and then recalculate below
  
  arma::mat WTW=trans(W)*W;
  
  // This should be IR = (XTX + P)^{-1} - The posterior variance-Covariance
  
  arma::mat IR=arma::inv(trimatu(chol(trans(W)*W)));
  
  // This should be posterior mean [Keep old code in case we want to return to it]
  
  // arma::mat b2=(IR*trans(IR))*(trans(W)*z);
  // b2.print("optimized value from old code")  ; 
  
  // arma::mat b2=(IR*trans(IR))*(trans(W)*z);
  arma::mat b2(b2a.begin(), b2a.nrow(), b2a.ncol(), false);
  //b2.print("optimized value from new code")  ; 
  
  NumericMatrix out(n,l1);
  arma::mat out2(out.begin(), out.nrow(), out.ncol(), false);
  NumericVector draws(n,1.0);
  NumericVector LL(n);
  
  NumericMatrix U1(l1,n);
  arma::mat U2(U1.begin(), U1.nrow(), U1.ncol(), false);     
  
  for(i=0;i<n;i++){
    U1( _, i)=rnorm(l1);
    // Normal Draws are scaled by Posterior Variance
    out2.row(i)=trans(b2+IR*U2.col(i));

  }
  
  
  
  // Note: LL does not seem to be used by downstream functins so can likely be edited out and removed 
  // From output - It is recomputed by summary functions
  // Note: in rnnorm_reg_cpp, out is the transpose of what is here
  // Do temporary (slow) solution to try to ensure the below works
  // will come back and correct later
  
  
//    NumericMatrix out3(l1,n);   
  
//    int j;
//    for(i=0;i<n;i++){
//      for(j=0;j<l1;j++){
//      out3(j,i)=out(i,j);  // copy to new matrix (do more elegant solution later) 
//      }
//     LL[i]=as<double>(f1(_["b"]=out3(_,i),_["y"]=y,_["x"]=x,offset2,wt2)); // Calculate log_likelihood
//   }
  
  //  Rcpp::Rcout << "LL is :"  << LL << std::endl;
  
  
  
  Rcpp::List Prior=Rcpp::List::create(Rcpp::Named("mean")=mu,Rcpp::Named("Precision")=P);  
  
  Rcpp::List outlist=Rcpp::List::create(
    Rcpp::Named("coefficients")=out,
    Rcpp::Named("coef.mode")=b2,
    Rcpp::Named("dispersion")=dispersion2,
    Rcpp::Named("Prior")=Prior,
    Rcpp::Named("prior.weights")=wt,
    Rcpp::Named("y")=y,
    Rcpp::Named("x")=x,
    Rcpp::Named("famfunc")=famfunc,
    Rcpp::Named("iters")=draws,
    Rcpp::Named("Envelope")=NULL,
    Rcpp::Named("loglike")=LL
  );  
  
  return(outlist);
  
}




