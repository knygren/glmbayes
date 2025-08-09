// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#pragma once
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <cmath>
#include "RcppArmadillo.h"
#include <limits>
#include <RcppParallel.h>
#define MATHLIB_STANDALONE
#include "nmath_local.h"
#include "dpq_local.h"



using namespace Rcpp;
using namespace RcppParallel;

#include "OpenCL_helper.h"
using namespace OpenCLHelper;


#include "OpenCL_proxy_kernels.h"
using namespace OpenCLProxy;



void progress_bar2(double x, double N);


//NumericVector dbinom_glmb( NumericVector x, NumericVector N, NumericVector means, int lg){
//    int n = x.size() ;
//    NumericVector res(n) ;
//    for( int i=0; i<n; i++) res[i] = R::dbinom( round(x[i]*N[i]),round(N[i]), means[i], lg ) ;
//    return res ;
//}

//constexpr double MY_PI = 3.14159265358979323846;
//constexpr double MY_2PI = 2.0 * MY_PI;


//extern "C" {
//  double dbinom_raw(double x, double n, double p, double q, int give_log);
//}

/////////////////////////////////////////////////////////////


inline double dbinom_raw_local(double x, double n, double p, double q, int give_log) {
  // Pass-through version for now
  return dbinom_raw(x, n, p, q, give_log);
}


///////////////////////////////////////////////////////////

void neg_dbinom_glmb_rmat(const RVector<double>& x,         // success proportion
                          const RVector<double>& N,         // trial count
                          const std::vector<double>& p_vec, // success probabilities
                          std::vector<double>& res,         // output buffer
                          const int lg)                     // log=TRUE?
{
  std::size_t n = x.size();
  if (res.size() != n)
    res.resize(n);  // ensure buffer size
  
  for (std::size_t i = 0; i < n; ++i) {
    double trials  = std::round(N[i]);
    double success = std::round(x[i] * N[i]);
    double p       = p_vec[i];
    double q       = 1.0 - p;
    
    // Thread-safe backend log-density
    res[i] = -dbinom_raw_local(success, trials, p, q, lg);
  }
}



// void neg_dbinom_glmb_rmat_old(const RVector<double>& x,         // success proportion
//                           const RVector<double>& N,         // trial count
//                           const RVector<double>& p_vec,     // standard-scale probabilities
//                           RVector<double>& res,             // output buffer
//                           const int lg) {
//   std::size_t n = x.length();
//   
//   for (std::size_t i = 0; i < n; ++i) {
//     // Input transformation: round values like dbinom_glmb()
//     double trials  = std::round(N[i]);
//     double success = std::round(x[i] * N[i]);
//     
//     double p = p_vec[i];         // success probability
//     double q = 1.0 - p;          // failure probability
//   
// //  Rcpp::Rcout << "success: " << success <<  "\n";
// //  Rcpp::Rcout << "trials: " << trials <<  "\n";
// //  Rcpp::Rcout << "p: " << p <<  "\n";
// //  Rcpp::Rcout << "q: " << q <<  "\n";
//   
//     
//     // Evaluate negative log-likelihood using R-accurate approximation
// //    res[i] = -dbinom_raw_log(success, trials, p, q);
//     
// //         Rcpp::Rcout << "c++: " << res[i] <<  "\n";
//     
//     res[i] = - dbinom_raw(success, trials, p, q, 1);  // log=TRUE
// 
//     //Rcpp::Rcout << "dbinom_raw: " << res[i] <<  "\n";
//     
//     //res[i] = - dbinom_raw_local(success, trials, p, q, 1);  // log=TRUE
//     
// //    Rcpp::Rcout << "dbinom_raw_local: " << res[i] <<  "\n";
//     
//   //      res[i] = - R::dbinom(success, trials, p, true);
//   //  Rcpp::Rcout << "R: " << res[i] <<  "\n";
//     
//       }
// }
// 


////////////////////////////////////////////////////////////






NumericVector dbinom_glmb(NumericVector x, NumericVector N, NumericVector means, int lg) {
  int n = x.size();
  NumericVector res(n);
  
  for (int i = 0; i < n; i++) {
    // Round to nearest integer for trial count and success count
    int trials  = static_cast<int>(std::round(N[i]));
    int success = static_cast<int>(std::round(x[i] * N[i]));
    
    // Clamp probabilities to avoid extreme values
    double p = std::min(1.0, std::max(0.0, means[i]));
    
    // Evaluate binomial log-likelihood
    res[i] = R::dbinom(success, trials, p, lg);

    
      }
  
  return res;
}



NumericVector dbinom_glmb2( NumericVector x, NumericVector N, NumericVector means, int lg){
  
  int i;
  int n = x.size() ;
  NumericVector res =no_init(n) ;
  NumericVector y=no_init(n);
  
  double nmax=max(N);
  double nmin=min(N);
  

  if(nmax==nmin){
    
  if(nmax==1)    std::transform( x.begin(), x.end(), means.begin(), res.begin(), [=](double y1, double means1){ return R::dbinom(y1,nmax, means1,lg); }); 

  if(nmax>1){
  for(i=1;i<n;i++){y[i]=round(x[i]*N[i]);  }
  std::transform( y.begin(), y.end(), means.begin(), res.begin(), [=](double y1, double means1){ return R::dbinom(y1,N[0], means1,lg); }); 
  }
  
  if(nmax<1){
    for(i=1;i<n;i++){y[i]=round(x[i]*N[i]);  }
    std::transform( y.begin(), y.end(), means.begin(), res.begin(), [=](double y1, double means1){ return R::dbinom(y1,N[0], means1,lg); }); 
  }
  
  }  
  
  
  if(nmax>nmin){
  for(  i=0; i<n; i++) {  res[i] = R::dbinom( round(x[i]*N[i]),round(N[i]), means[i], lg  );  } 
  }
  
    return (res) ;
}


NumericMatrix dbinom_glmb3( NumericMatrix x, NumericMatrix N, NumericMatrix means, int lg){
  
  int n1 = x.rows() ;
  int n2 = x.cols() ;
  NumericMatrix res =no_init(n1,n2) ;

  double nmax=max(N);
  double nmin=min(N);

  //y=round(x);
  

  // For now, deal only with case where nmax=nmax=1
  // return negative log-likelihood
  
  if(nmax==nmin){
    

        std::transform( x.begin(), x.end(), means.begin(), res.begin(), [=](double y1, double means1){ return -R::dbinom(y1,N[0], means1,lg); });

  }  
  
  

  return (res) ;
}





NumericVector cpprbinom2(int n, double size, NumericVector prob) { 
  NumericVector v = no_init(n);
  
  std::transform( prob.begin(), prob.end(), v.begin(), [=](double p){ return R::rbinom(size, p); }); 

    return(v);}



///////////////////////// Logit Functions ///////////////////////////////////////

NumericVector  f1_binomial_logit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);


    for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
  
    
    xb2=exp(-alpha2- x2 * b2);

    for(int j=0;j<l1;j++){
    xb(j)=1/(1+xb(j));
    }

    yy=-dbinom_glmb(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);

    }
    
    return res;      
}




NumericVector  f2_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
     //   Note: Does not seem to be used-- Editing out         
  //  NumericVector invwt=1/sqrt(wt);

    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);

    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 

    double res1=0;


    for(int i=0;i<m1;i++){
      
      
      Rcpp::checkUserInterrupt();
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 

    bmu2=b2-mu2;
        
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
    xb2=exp(-alpha2- x2 * b2);

    for(int j=0;j<l1;j++){
    xb(j)=1/(1+xb(j));
    }


    yy=-dbinom_glmb(y,wt,xb,true);


    res(i) =std::accumulate(yy.begin(), yy.end(), res1);

    }
    
    return res;      
}


//#include <RcppArmadillo.h>
//using namespace Rcpp;
//using namespace arma;

// arma::vec f2_binomial_logit_arma(NumericMatrix b, NumericVector y,
//                                 NumericMatrix x, NumericMatrix mu,
//                                 NumericMatrix P, NumericVector alpha,
//                                 NumericVector wt, int progbar = 0) {
//   int l1 = x.nrow(), l2 = x.ncol();
//   int m1 = b.ncol();
//   
//   arma::mat b2full(b.begin(), l2, m1, false);  // Shallow view over b
//   
//   Rcpp::NumericMatrix b2temp(l2, 1);
//   arma::mat x2(x.begin(), l1, l2, false);
//   arma::mat alpha2(alpha.begin(), l1, 1, false);
//   
//   Rcpp::NumericVector xb(l1);
//   arma::colvec xb2(xb.begin(), l1, false);
//   
//   NumericVector yy(l1);
// //  NumericVector res(m1);
//   NumericMatrix bmu(l2, 1);
//   arma::mat mu2(mu.begin(), l2, 1, false);
//   arma::mat bmu2(bmu.begin(), l2, 1, false);
// 
//   arma::vec res2(m1);  // Owning allocation
//   
//   //std::vector<double> xb_temp(l1), yy_temp(l1);
//   //arma::colvec xb_temp2(xb_temp.data(), l1, false);  // shallow Armadillo view
//   
//   
//   
//   for (int i = 0; i < m1; i++) {
//     Rcpp::checkUserInterrupt();
//     if (progbar == 1) {
//       progress_bar2(i, m1 - 1);
//       if (i == m1 - 1) { Rcpp::Rcout << "" << std::endl; }
//     }
//     
//     b2temp = b(Range(0, l2 - 1), Range(i, i));
//     //arma::mat b2(b2temp.begin(), l2, 1, false);
//     arma::mat b2(b2full.colptr(i), l2, 1, false);  // View of column i, size l2 × 1
//     
//     arma::mat P2(P.begin(), l2, l2, false);
//     
//     bmu2 = b2 - mu2;
//     double res1 = 0.5 * arma::as_scalar(bmu2.t() * P2 * bmu2);
//     
//     xb2 = exp(-alpha2 - x2 * b2);  // eta = -alpha - Xb
//  
//     for (int j = 0; j < l1; j++) {
//       xb(j) = 1.0 / (1.0 + xb(j));  // logistic link
//     }
// 
//     
//     // Wrap existing R memory buffers
//     RcppParallel::RVector<double> y_view(y);
//     RcppParallel::RVector<double> wt_view(wt);
//     RcppParallel::RVector<double> xb_view(xb);
//     RcppParallel::RVector<double> yy_view(yy); //must be preallocated to match y.size()
//     
//     // In-place evaluation using your log-scale accurate backend
//     neg_dbinom_glmb_rmat_old(y_view, wt_view, xb_view, yy_view,1.0);
//         
// //    yy=-dbinom_glmb(y,wt,xb,true);
//     
//     
//     //res(i) =std::accumulate(yy.begin(), yy.end(), res1);
//     res2(i) = std::accumulate(yy.begin(), yy.end(), res1);
// 
//       }
//   
//   return res2;
// }


arma::vec f2_binomial_logit_rmat(
    // NumericMatrix b, NumericVector y,
    //                              NumericMatrix x, NumericMatrix mu,
    //                              NumericMatrix P, NumericVector alpha,
    //                              NumericVector wt, int progbar = 0
                                   const RMatrix<double>& b,
                                   const RVector<double>& y,
                                   const RMatrix<double>& x,
                                   const RMatrix<double>& mu,
                                   const RMatrix<double>& P,
                                   const RVector<double>& alpha,
                                   const RVector<double>& wt,
                                   const int progbar=0   
                                   ) {
//  int l1 = x.nrow(), l2 = x.ncol();
//  int m1 = b.ncol();
  
  /////////////////////////////////////////////////////////////////////////////
  
  
  
  std::size_t l1 = x.nrow();
  std::size_t l2 = x.ncol();
  std::size_t m1 = b.ncol();
  
  // Armadillo views over RMatrix memory (using pointer cast for compatibility)
  arma::mat b2full(const_cast<double*>(&*b.begin()), l2, m1, false);
  arma::mat x2(const_cast<double*>(&*x.begin()), l1, l2, false);
  arma::mat mu2(const_cast<double*>(&*mu.begin()), l2, 1, false);
  arma::mat P2(const_cast<double*>(&*P.begin()), l2, l2, false);
  arma::mat alpha2(const_cast<double*>(&*alpha.begin()), l1, 1, false);
  
  arma::vec res(m1, arma::fill::none);
  arma::mat bmu(l2, 1, arma::fill::none);
  
  
  std::vector<double> xb_temp(l1), yy_temp(l1);
  arma::colvec xb_temp2(xb_temp.data(), l1, false);  // shallow Armadillo view
  
  ////////////////////////////////////////////////////////////

  
  for (std::size_t i = 0; i < m1; ++i) {

    arma::mat b_i(b2full.colptr(i), l2, 1, false);
    
    bmu = b_i - mu2;
    
    double mahal = 0.5 * arma::as_scalar(bmu.t() * P2 * bmu);
    
  
    xb_temp2 = arma::exp(-alpha2 - x2 * b_i);
    
    for (std::size_t  j = 0; j < l1; j++) {

      xb_temp[j] = 1.0 / (1.0 + xb_temp[j]);  // logistic link
      
    }
    

    // In-place evaluation using your log-scale accurate backend
  neg_dbinom_glmb_rmat(y, wt, xb_temp, yy_temp,1.0);


    res(i) =std::accumulate(yy_temp.begin(), yy_temp.end(), mahal);

  }
  
  return res;
}




arma::mat  f3_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat y2(y.begin(), l1, 1, false);
    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
       
   
       
    NumericMatrix Ptemp(l1,l1);  
      
    for(int i=0;i<l1;i++){
     Ptemp(i,i)=wt(i); 
    }  
    
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);
    NumericMatrix out(l2,m1);


    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 
    arma::mat Ptemp2(Ptemp.begin(), l1, l1, false);
    arma::mat out2(out.begin(), l2, m1, false);
    
    NumericMatrix::Column outtemp=out(_,0);
    //NumericMatrix::Row outtempb=out(0,_);

    arma::mat outtemp2(outtemp.begin(),1,l2,false);
    //arma::mat outtempb2(outtempb.begin(),1,l2,false);
    
    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
    if(progbar==1){ 
    progress_bar2(i, m1-1);
    if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
    };  
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);

    bmu2=b2-mu2;

//  		p<-1/(1+t(exp(-alpha-x%*%b)))
//		t(x)%*%((t(p)-y)*wt)+P%*%(b-mu)


    xb2=exp(-alpha2- x2 * b2);
    
    

    for(int j=0;j<l1;j++){
    xb(j)=1/(1+xb(j));  
    xb(j)=(xb(j)-y(j))*wt(j);
    }



    outtemp2= P2 * bmu2+x2.t() * xb2;
    }
    
   // return  b;
  
    return trans(out2);      
}




///////////////////////// Probit Functions ///////////////////////////////////////

NumericVector  f1_binomial_probit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);


    for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
 
    xb2=alpha2+ x2 * b2;   
    xb=pnorm(xb,0.0,1.0);
	

    yy=-dbinom_glmb(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);

    }
    
    return res;      
}




NumericVector  f2_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
     //   Note: Does not seem to be used-- Editing out        
//    NumericVector invwt=1/sqrt(wt);

    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);

    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 

    double res1=0;


    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 

    bmu2=b2-mu2;
        
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
  
    xb2=alpha2+ x2 * b2;   
    xb=pnorm(xb,0.0,1.0);

    yy=-dbinom_glmb(y,wt,xb,true);


    res(i) =std::accumulate(yy.begin(), yy.end(), res1);

    }
    
    return res;      
}


// 
// arma::vec  f2_binomial_probit_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
// {
//   
//   // Get dimensions of x - Note: should match dimensions of
//   //  y, b, alpha, and wt (may add error checking)
//   
//   // May want to add method for dealing with alpha and wt when 
//   // constants instead of vectors
//   
//   int l1 = x.nrow(), l2 = x.ncol();
//   int m1 = b.ncol();
//   
//   //    int lalpha=alpha.nrow();
//   //    int lwt=wt.nrow();
//   
//   arma::mat b2full(b.begin(), l2, m1, false);  // Shallow view over b
//   Rcpp::NumericMatrix b2temp(l2,1);
//   
//   arma::mat x2(x.begin(), l1, l2, false); 
//   arma::mat alpha2(alpha.begin(), l1, 1, false); 
//   
//   Rcpp::NumericVector xb(l1);
//   arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
//   
//   //   Note: Does not seem to be used-- Editing out        
//   //    NumericVector invwt=1/sqrt(wt);
//   
//   // Moving Loop inside the function is key for speed
//   
//   NumericVector yy(l1);
// //  NumericVector res(m1);
//   NumericMatrix bmu(l2,1);
//   
//   arma::mat mu2(mu.begin(), l2, 1, false); 
//   arma::mat bmu2(bmu.begin(), l2, 1, false); 
//   arma::vec res2(m1);  // Owning allocation
//   
//   double res1=0;
//   
// //  Rcpp::Rcout << "m1=" << m1 <<  "\n";
//   
//     
//   for(int i=0;i<m1;i++){
//     Rcpp::checkUserInterrupt();
//     if(progbar==1){ 
//       progress_bar2(i, m1-1);
//       if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
//     };  
//     
//     
//     
//     b2temp=b(Range(0,l2-1),Range(i,i));
//     //arma::mat b2(b2temp.begin(), l2, 1, false); 
//     arma::mat b2(b2full.colptr(i), l2, 1, false);  // View of column i, size l2 × 1
//     
//     arma::mat P2(P.begin(), l2, l2, false); 
//     
//     bmu2=b2-mu2;
//     
//     res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
//     
//     
//     xb2=alpha2+ x2 * b2;   
//     
// //    Rcpp::Rcout << "xb_arma_in: " << "\n";
//     
// //    Rcpp::Rcout << xb << "\n";
//     
//     
//     xb=pnorm(xb,0.0,1.0);
// 
// //    Rcpp::Rcout << "xb_arma_out: " << "\n";
//     
// //    Rcpp::Rcout << xb << "\n";
//     
//     
//     // Wrap existing R memory buffers
//     RcppParallel::RVector<double> y_view(y);
//     RcppParallel::RVector<double> wt_view(wt);
//     RcppParallel::RVector<double> xb_view(xb);
//     RcppParallel::RVector<double> yy_view(yy); //must be preallocated to match y.size()
//     
//     // In-place evaluation using your log-scale accurate backend
//     
// //    Rcpp::Rcout << " Entering neg_dbin_glmb_rmat  " << "\n";
// 
// 
// //  Rcout << "y: " << y << std::endl;
// //  Rcout << "wt: " << wt << std::endl;
// //  Rcout << "xb: " << xb << std::endl;
// 
// 
//   
//     neg_dbinom_glmb_rmat_old(y_view, wt_view, xb_view, yy_view,1.0);
// 
//     
//     
// //    Rcout << "yy: " << yy << std::endl;
//     
//         
// //    Rcpp::Rcout << " Exiting neg_dbin_glmb_rmat  " << "\n";
//     
//     //yy=-dbinom_glmb(y,wt,xb,true);
//     
//     
//     //res(i) =std::accumulate(yy.begin(), yy.end(), res1);
//     res2(i) = std::accumulate(yy.begin(), yy.end(), res1);
//     
// //    Rcpp::Rcout << "i=" << i << " draws=" << draws[i] << "\n";
//   
// //    Rcpp::Rcout << "res2=" << res2(i) <<  "\n";
//   
//   }
// 
//   //    Rcpp::Rcout << " Returning res2  " << "\n";
//   //  Rcpp::Rcout << "m1=" << m1 <<  "\n";
//   
//   return res2;      
// }


arma::vec f2_binomial_probit_rmat(
    // NumericMatrix b, NumericVector y,
    //                              NumericMatrix x, NumericMatrix mu,
    //                              NumericMatrix P, NumericVector alpha,
    //                              NumericVector wt, int progbar = 0
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    const int progbar=0   
) {
  //  int l1 = x.nrow(), l2 = x.ncol();
  //  int m1 = b.ncol();
  
  /////////////////////////////////////////////////////////////////////////////
  
  
  
  std::size_t l1 = x.nrow();
  std::size_t l2 = x.ncol();
  std::size_t m1 = b.ncol();
  
  // Armadillo views over RMatrix memory (using pointer cast for compatibility)
  arma::mat b2full(const_cast<double*>(&*b.begin()), l2, m1, false);
  arma::mat x2(const_cast<double*>(&*x.begin()), l1, l2, false);
  arma::mat mu2(const_cast<double*>(&*mu.begin()), l2, 1, false);
  arma::mat P2(const_cast<double*>(&*P.begin()), l2, l2, false);
  arma::mat alpha2(const_cast<double*>(&*alpha.begin()), l1, 1, false);
  
  arma::vec res(m1, arma::fill::none);
  arma::mat bmu(l2, 1, arma::fill::none);
  
  
  std::vector<double> xb_temp(l1), yy_temp(l1);
  arma::colvec xb_temp2(xb_temp.data(), l1, false);  // shallow Armadillo view
  
  ////////////////////////////////////////////////////////////
  
  
  for (std::size_t i = 0; i < m1; ++i) {
    
    arma::mat b_i(b2full.colptr(i), l2, 1, false);
    
    bmu = b_i - mu2;
    
    double mahal = 0.5 * arma::as_scalar(bmu.t() * P2 * bmu);
    
    
    xb_temp2 = alpha2+  x2 * b_i;
    
    for (std::size_t j = 0; j < l1; j++) {
      xb_temp[j] = pnorm5_local(xb_temp[j], 0.0, 1.0, 1, 0);
    }
    

    // In-place evaluation using your log-scale accurate backend
    neg_dbinom_glmb_rmat(y, wt, xb_temp, yy_temp,1.0);
    
    
    res(i) =std::accumulate(yy_temp.begin(), yy_temp.end(), mahal);
    
  }
  
  return res;
}




arma::mat  f3_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat y2(y.begin(), l1, 1, false);
    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
       
   
       
    NumericMatrix Ptemp(l1,l1);  
      
    for(int i=0;i<l1;i++){
     Ptemp(i,i)=wt(i); 
    }  
    
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);
    NumericMatrix out(l2,m1);
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector d1(l1);


    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 
    arma::mat Ptemp2(Ptemp.begin(), l1, l1, false);
    arma::mat out2(out.begin(), l2, m1, false);
    
    NumericMatrix::Column outtemp=out(_,0);
    //NumericMatrix::Row outtempb=out(0,_);

    arma::mat outtemp2(outtemp.begin(),1,l2,false);
    //arma::mat outtempb2(outtempb.begin(),1,l2,false);
    
    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);

    bmu2=b2-mu2;
    
    xb2=alpha2+ x2 * b2;
    p1=pnorm(xb,0.0,1.0);
    p2=pnorm(-xb,0.0,1.0);
    d1=dnorm(xb,0.0,1.0);
    
//    -t(x)%*%as.matrix(((y*dnorm(alpha+x%*%b)/p1)-(1-y)*dnorm(alpha+x%*%b)/p2)*wt)+P%*%(b-mu)


    for(int j=0;j<l1;j++){
    xb(j)=(y(j)*d1(j)/p1(j)-(1-y(j))*d1(j)/p2(j))*wt(j);    
    }


    outtemp2= P2 * bmu2-x2.t() * xb2;
    }
    
   // return  b;
  
    return trans(out2);      
}



///////////////////////// cLOGLOG FUNCTION ///////////////////////////////////////

NumericVector  f1_binomial_cloglog(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);


    for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
 
    xb2=alpha2+ x2 * b2;   
    xb=exp(-exp(xb));
  
    for(int j=0;j<l1;i++){
    xb(j)=1-xb(j);
    }
  

    yy=-dbinom_glmb(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);

    }
    
    return res;      
}




NumericVector  f2_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
//   Note: Does not seem to be used-- Editing out      
//    NumericVector invwt=1/sqrt(wt);

    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);

    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 

    double res1=0;


    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 

    bmu2=b2-mu2;
        
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
    xb2=alpha2+ x2 * b2;   
    
    for(int j=0;j<l1;j++){
    xb(j)=1-exp(-exp(xb(j)));
    }

    yy=-dbinom_glmb(y,wt,xb,true);


    res(i) =std::accumulate(yy.begin(), yy.end(), res1);

    }
    
    return res;      
}



// 
// arma::vec  f2_binomial_cloglog_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
// {
//   
//   // Get dimensions of x - Note: should match dimensions of
//   //  y, b, alpha, and wt (may add error checking)
//   
//   // May want to add method for dealing with alpha and wt when 
//   // constants instead of vectors
//   
//   int l1 = x.nrow(), l2 = x.ncol();
//   int m1 = b.ncol();
//   
//   //    int lalpha=alpha.nrow();
//   //    int lwt=wt.nrow();
//   
//   arma::mat b2full(b.begin(), l2, m1, false);  // Shallow view over b
//   Rcpp::NumericMatrix b2temp(l2,1);
//   
//   arma::mat x2(x.begin(), l1, l2, false); 
//   arma::mat alpha2(alpha.begin(), l1, 1, false); 
//   
//   Rcpp::NumericVector xb(l1);
//   arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
//   
//   //   Note: Does not seem to be used-- Editing out      
//   //    NumericVector invwt=1/sqrt(wt);
//   
//   // Moving Loop inside the function is key for speed
//   
//   NumericVector yy(l1);
// //  NumericVector res(m1);
//   NumericMatrix bmu(l2,1);
//   
//   arma::mat mu2(mu.begin(), l2, 1, false); 
//   arma::mat bmu2(bmu.begin(), l2, 1, false); 
//   arma::vec res2(m1);  // Owning allocation
//   
//   double res1=0;
//   
//   
//   for(int i=0;i<m1;i++){
//     Rcpp::checkUserInterrupt();
//     
//     if(progbar==1){ 
//       progress_bar2(i, m1-1);
//       if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
//     };  
//     
//     
//     b2temp=b(Range(0,l2-1),Range(i,i));
//     //arma::mat b2(b2temp.begin(), l2, 1, false); 
//     arma::mat b2(b2full.colptr(i), l2, 1, false);  // View of column i, size l2 × 1
//     
//     arma::mat P2(P.begin(), l2, l2, false); 
//     
//     bmu2=b2-mu2;
//     
//     res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
//     
//     xb2=alpha2+ x2 * b2;   
//     
//     for(int j=0;j<l1;j++){
//       xb(j)=1-exp(-exp(xb(j)));
//     }
//     
//     
//     
//     // Wrap existing R memory buffers
//     RcppParallel::RVector<double> y_view(y);
//     RcppParallel::RVector<double> wt_view(wt);
//     RcppParallel::RVector<double> xb_view(xb);
//     RcppParallel::RVector<double> yy_view(yy); //must be preallocated to match y.size()
//     
//     // In-place evaluation using your log-scale accurate backend
//     neg_dbinom_glmb_rmat_old(y_view, wt_view, xb_view, yy_view,1.0);
//     
//     
//     //yy=-dbinom_glmb(y,wt,xb,true);
//     
//     
//     //res(i) =std::accumulate(yy.begin(), yy.end(), res1);
//     res2(i) = std::accumulate(yy.begin(), yy.end(), res1);
//     
//   }
//   
//   return res2;      
// }


arma::vec f2_binomial_cloglog_rmat(
    // NumericMatrix b, NumericVector y,
    //                              NumericMatrix x, NumericMatrix mu,
    //                              NumericMatrix P, NumericVector alpha,
    //                              NumericVector wt, int progbar = 0
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    const int progbar=0   
) {
  //  int l1 = x.nrow(), l2 = x.ncol();
  //  int m1 = b.ncol();
  
  /////////////////////////////////////////////////////////////////////////////
  
  
  
  std::size_t l1 = x.nrow();
  std::size_t l2 = x.ncol();
  std::size_t m1 = b.ncol();
  
  // Armadillo views over RMatrix memory (using pointer cast for compatibility)
  arma::mat b2full(const_cast<double*>(&*b.begin()), l2, m1, false);
  arma::mat x2(const_cast<double*>(&*x.begin()), l1, l2, false);
  arma::mat mu2(const_cast<double*>(&*mu.begin()), l2, 1, false);
  arma::mat P2(const_cast<double*>(&*P.begin()), l2, l2, false);
  arma::mat alpha2(const_cast<double*>(&*alpha.begin()), l1, 1, false);
  
  arma::vec res(m1, arma::fill::none);
  arma::mat bmu(l2, 1, arma::fill::none);
  
  
  std::vector<double> xb_temp(l1), yy_temp(l1);
  arma::colvec xb_temp2(xb_temp.data(), l1, false);  // shallow Armadillo view
  
  ////////////////////////////////////////////////////////////
  
  
  for (std::size_t i = 0; i < m1; ++i) {
    
    arma::mat b_i(b2full.colptr(i), l2, 1, false);
    
    bmu = b_i - mu2;
    
    double mahal = 0.5 * arma::as_scalar(bmu.t() * P2 * bmu);
    
    
    xb_temp2 = alpha2+  x2 * b_i;
    
    for (std::size_t j = 0; j < l1; j++) {
      xb_temp[j] =1-  exp(-exp(xb_temp[j]));
    }
    

//    for(int j=0;j<l1;j++){
//      xb(j)=1-exp(-exp(xb(j)));
//    }
    
        
    // In-place evaluation using your log-scale accurate backend
    neg_dbinom_glmb_rmat(y, wt, xb_temp, yy_temp,1.0);
    
    
    res(i) =std::accumulate(yy_temp.begin(), yy_temp.end(), mahal);
    
  }
  
  return res;
}




arma::mat  f3_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat y2(y.begin(), l1, 1, false);
    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
       
   
       
    NumericMatrix Ptemp(l1,l1);  
      
    for(int i=0;i<l1;i++){
     Ptemp(i,i)=wt(i); 
    }  
    
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);
    NumericMatrix out(l2,m1);
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector atemp(l1);


    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 
    arma::mat Ptemp2(Ptemp.begin(), l1, l1, false);
    arma::mat out2(out.begin(), l2, m1, false);
    
    NumericMatrix::Column outtemp=out(_,0);
    //NumericMatrix::Row outtempb=out(0,_);

    arma::mat outtemp2(outtemp.begin(),1,l2,false);
    //arma::mat outtempb2(outtempb.begin(),1,l2,false);
    
    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);

    bmu2=b2-mu2;
    
    xb2=alpha2+ x2 * b2;

    for(int j=0;j<l1;j++){
    p1(j)=1-exp(-exp(xb(j)));
    p2(j)=exp(-exp(xb(j)));
    atemp(j)=exp(xb(j)-exp(xb(j)));
    xb(j)=((y(j)*atemp(j)/p1(j))-((1-y(j))*atemp(j)/p2(j)))*wt(j);    
    }


    outtemp2= P2 * bmu2-x2.t() * xb2;
    }
    
  
    return trans(out2);      
}


////////////////////////////// Functions preparing for OpenCL implementation



// [[Rcpp::export]]
List f2_binomial_logit_prep(
    NumericMatrix b,
    NumericVector y,
    NumericMatrix x,
    NumericMatrix mu,
    NumericMatrix P,
    NumericVector alpha,
    NumericVector wt,
    int progbar
) {
  int l1  = x.nrow();
  int l2  = x.ncol();
  int m1  = b.ncol();
  
  // Armadillo “views” on the R objects
  arma::mat    X   (x.begin(),     l1, l2, false);
  arma::colvec Y   (y.begin(),     l1,   false);
  arma::colvec A   (alpha.begin(), l1,   false);
  arma::colvec Wt  (wt.begin(),    l1,   false);
  arma::colvec Mu  (mu.begin(),    l2,   false);
  arma::mat    P2  (P.begin(),      l2, l2, false);
  
  // Outputs
  NumericMatrix xb_mat(l1, m1);
  NumericVector  qf    (m1);
  
  // Temps
  arma::colvec b2   (l2);
  arma::colvec bmu2 (l2);
  arma::colvec xb2  (l1);
  
  for (int i = 0; i < m1; ++i) {
    Rcpp::checkUserInterrupt();
    if (progbar > 0) {
      progress_bar2(i, m1 - 1);
      if (i == m1 - 1) Rcpp::Rcout << std::endl;
    }
    
    // slice b[, i]
    std::copy(
      b.begin() + i * l2,
      b.begin() + (i + 1) * l2,
      b2.begin()
    );
    
    // 1) prior quadratic form
    bmu2  = b2 - Mu;
    qf[i] = 0.5 * arma::as_scalar(bmu2.t() * P2 * bmu2);
    
    // 2) logistic prep: π = 1/(1 + exp(−(α + X·b2)))
    xb2 = arma::exp(- A - (X * b2));
    for (int j = 0; j < l1; ++j) {
      xb_mat(j, i) = 1.0 / (1.0 + xb2(j));
    }
  }
  
  return List::create(
    Rcpp::Named("xb") = xb_mat,
    Rcpp::Named("qf") = qf
  );
}


List f2_binomial_logit_prep_v2(
    NumericMatrix b,        // (l2 × m1) grid of β
    NumericVector y,        // (l1) responses (unused in prep)
    NumericMatrix x,        // (l1 × l2) design matrix
    NumericMatrix mu,       // (l2 × 1) prior means
    NumericMatrix P,        // (l2 × l2) precision matrix
    NumericVector alpha,    // (l1) offsets
    NumericVector wt,       // (l1) weights (unused)
    int progbar = 0
) {
  //--------------------------------------------------------------------------
  // 1) INPUT MARSHALLING: Rcpp types → flat arrays (host buffers)
  //--------------------------------------------------------------------------
  int l1 = x.nrow();     // number of observations
  int l2 = x.ncol();     // number of coefficients
  int m1 = b.ncol();     // number of β grid points
  
  // flatten into column-major contiguous vectors
  std::vector<double> X_flat     = flattenMatrix(x);      // length = l1*l2
  std::vector<double> B_flat     = flattenMatrix(b);      // length = l2*m1
  std::vector<double> mu_flat    = flattenMatrix(mu);     // length = l2
  std::vector<double> P_flat     = flattenMatrix(P);      // length = l2*l2
  std::vector<double> alpha_flat = copyVector(alpha);     // length = l1
  std::vector<double> wt_flat    = copyVector(wt);        // length = l1
  
  // Kernel inputs are now:
  //   X_flat, B_flat, mu_flat, P_flat, alpha_flat, wt_flat
  //   plus scalars l1, l2, m1
  
  //--------------------------------------------------------------------------
  // 2) KERNEL STUB: CPU fallback, identical math to eventual OpenCL kernel
  //--------------------------------------------------------------------------
  // Allocate flat outputs
  std::vector<double> qf_flat(m1, 0.0);           // length = m1
  std::vector<double> xb_flat((size_t)l1*m1, 0.0); // length = l1*m1
  
  // Temporaries
  std::vector<double> bmu(l2);
  std::vector<double> tmp(l2);
  
  for (int j = 0; j < m1; ++j) {
    Rcpp::checkUserInterrupt();
    if (progbar > 0) {
      progress_bar2(j, m1 - 1);
      if (j == m1 - 1) Rcpp::Rcout << std::endl;
    }
    
    // pointers into flattened buffers
    const double* bj   = B_flat.data()   + (size_t)j * l2;
    const double* muj  = mu_flat.data();
    
    // -- a) Quadratic form: 0.5 * (b_j – mu)' P (b_j – mu)
    for (int k = 0; k < l2; ++k) {
      bmu[k] = bj[k] - muj[k];
    }
    // tmp = P * bmu
    for (int r = 0; r < l2; ++r) {
      double acc = 0.0;
      size_t prow = (size_t)r * l2;
      for (int c = 0; c < l2; ++c) {
        acc += P_flat[prow + c] * bmu[c];
      }
      tmp[r] = acc;
    }
    // dot(bmu, tmp)
    double qval = 0.0;
    for (int k = 0; k < l2; ++k) {
      qval += bmu[k] * tmp[k];
    }
    qf_flat[j] = 0.5 * qval;
    
    // -- b) Logistic prep: π_{i,j} = 1 / (1 + exp(α[i] + X[i,·]·b_j))
    for (int i = 0; i < l1; ++i) {
      // compute linear predictor: α_i + sum_k X[i,k]*b_j[k]
      double linpred = alpha_flat[i];
      size_t xoff = (size_t)i;
      for (int k = 0; k < l2; ++k) {
        // X is column-major: element (i,k) at X_flat[k*l1 + i]
        linpred += X_flat[(size_t)k * l1 + xoff] * bj[k];
      }
      // logistic
      xb_flat[(size_t)j * l1 + i] = 1.0 / (1.0 + std::exp(-linpred));
    }
  }
  
  // Kernel outputs are:
  //   qf_flat  (length = m1)
  //   xb_flat  (length = l1 * m1)
  
  //--------------------------------------------------------------------------
  // 3) OUTPUT MARSHALLING: flat arrays → Rcpp types
  //--------------------------------------------------------------------------
  NumericVector qf(m1);
  for (int j = 0; j < m1; ++j) {
    qf[j] = qf_flat[j];
  }
  
  NumericMatrix xb(l1, m1);
  for (int j = 0; j < m1; ++j) {
    for (int i = 0; i < l1; ++i) {
      xb(i, j) = xb_flat[(size_t)j * l1 + i];
    }
  }
  
  return List::create(
    Named("xb") = xb,   // l1 × m1 logistic probabilities
    Named("qf") = qf    // m1 prior quadratic forms
  );
}





// [[Rcpp::export]]
List f2_binomial_logit_prep_v3(
    NumericMatrix b,        // (l2 × m1) grid of β
    NumericVector y,        // (l1) responses (unused in prep)
    NumericMatrix x,        // (l1 × l2) design matrix
    NumericMatrix mu,       // (l2 × 1) prior means
    NumericMatrix P,        // (l2 × l2) precision matrix
    NumericVector alpha,    // (l1) offsets
    NumericVector wt,       // (l1) weights (unused)
    int progbar = 0
) {
  // 1) INPUT MARSHALLING
  int l1 = x.nrow();
  int l2 = x.ncol();
  int m1 = b.ncol();
  
  std::vector<double> X_flat     = flattenMatrix(x);   // length = l1*l2
  std::vector<double> B_flat     = flattenMatrix(b);   // length = l2*m1
  std::vector<double> mu_flat    = flattenMatrix(mu);  // length = l2
  std::vector<double> P_flat     = flattenMatrix(P);   // length = l2*l2
  std::vector<double> alpha_flat = copyVector(alpha);  // length = l1
  
  // 2) INVOKE PROXY KERNEL
  std::vector<double> qf_flat(m1);
  std::vector<double> xb_flat((size_t)l1 * m1);
  
  
  // —— DEBUG INPUTS ——  
  Rcpp::Rcout << "[DEBUG] l1=" << l1 
              << "  l2=" << l2 
              << "  m1=" << m1 << "\n";
  
  // show first 3 entries of each flat buffer  
  int K = std::min(3, m1);  
  Rcpp::Rcout << "X_flat[0..2]: ";  
  for (int i = 0; i < std::min(3, (int)X_flat.size()); ++i)  
    Rcpp::Rcout << X_flat[i] << " ";  
  Rcpp::Rcout << "\nB_flat[0..2]: ";  
  for (int i = 0; i < std::min(3, (int)B_flat.size()); ++i)  
    Rcpp::Rcout << B_flat[i] << " ";  
  Rcpp::Rcout << "\nmu_flat[0..2]: ";  
  for (int i = 0; i < std::min(3, (int)mu_flat.size()); ++i)  
    Rcpp::Rcout << mu_flat[i] << " ";  
  Rcpp::Rcout << "\nalpha_flat[0..2]: ";  
  for (int i = 0; i < std::min(3, (int)alpha_flat.size()); ++i)  
    Rcpp::Rcout << alpha_flat[i] << " ";  
  Rcpp::Rcout << "\n\n";  
  R_FlushConsole();
  // —— end DEBUG ——  
  
  
  f2_binomial_logit_prep_kernel_proxy(
    X_flat, B_flat, mu_flat, P_flat, alpha_flat,
    l1, l2, m1,
    qf_flat, xb_flat,
    progbar
  );
  
  

  
  
  // —— DEBUG OUTPUTS ——  
  Rcpp::Rcout << "[DEBUG] first " << K << " qf_flat: ";  
  for (int j = 0; j < K; ++j)  
    Rcpp::Rcout << qf_flat[j] << " ";  
  
  Rcpp::Rcout << "\n[DEBUG] first xb_flat(·,0): ";  
  for (int j = 0; j < K; ++j)  
    Rcpp::Rcout << xb_flat[(size_t)j * l1 + 0] << " ";  
  
  Rcpp::Rcout << "\n\n";  
  R_FlushConsole();  
  // —— end DEBUG ——  
  
  
  // 3) OUTPUT MARSHALLING
  NumericVector qf(m1);
  for (int j = 0; j < m1; ++j) {
    qf[j] = qf_flat[j];
  }
  
  NumericMatrix xb(l1, m1);
  for (int j = 0; j < m1; ++j) {
    for (int i = 0; i < l1; ++i) {
      xb(i, j) = xb_flat[(size_t)j * l1 + i];
    }
  }
  
  return List::create(
    Named("xb") = xb,
    Named("qf") = qf
  );
}


NumericVector f2_binomial_logit_accum(
    NumericMatrix xb,        // n × m matrix of π = P(y=1)
    NumericVector qf,        // length m: 0.5*(b-μ)'P(b-μ)
    NumericVector y,         // length n observed {0,1}
    NumericVector wt,        // length n weights
    int progbar          // 0 = no bar, 1 = show bar
) {
  int n = xb.nrow();
  int m = xb.ncol();
  NumericVector res(m);
  
  
  
  
  for (int i = 0; i < m; ++i) {
    Rcpp::checkUserInterrupt();
    //if (progbar == 1) {
    //  progress_bar2(i, m - 1);
    //  if (i == m - 1) Rcpp::Rcout << std::endl;
    //}
    
    // extract column i of xb
    NumericVector xbi = xb(_, i);
    
    // per-observation log-likelihoods
    NumericVector ll = dbinom_glmb(y, wt, xbi, true);
    
    // sum of log-likelihoods
    double sumll = std::accumulate(ll.begin(), ll.end(), 0.0);
    
    // total negative log-lik = quadratic form + (− sum log-lik)
    res[i] = qf[i] - sumll;
  }
  
  
  
  return res;
}