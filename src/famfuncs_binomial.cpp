// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <cmath>
#include "RcppArmadillo.h"
#include <limits>
#include <RcppParallel.h>
#include <Rmath.h>   // libR Mathlib: Rf_dbinom_raw, Rf_pnorm5, Rf_dnorm4, Rf_dbinom
#include "openclPort.h"
#include "famfuncs.h"
#include "progress_utils.h"

// Namespaces

using namespace Rcpp;
using namespace RcppParallel;
using namespace openclPort;
using namespace glmbayes::fam;
using namespace glmbayes::progress;



/////////////////////////////////////////////////////////////


namespace glmbayes{

namespace fam {
// Neg log binomial likelihood, vectorized
void neg_dbinom_glmb_rmat(const RVector<double>& x,          // success proportion
                          const RVector<double>& N,          // trial count
                          const std::vector<double>& p_vec,  // success probabilities
                          const std::vector<double>& q_vec,  // failure probabilities
                          std::vector<double>& res,          // output buffer
                          const int lg)                      // log=TRUE?
{
  std::size_t n = x.size();
  if (res.size() != n)
    res.resize(n);
  
  for (std::size_t i = 0; i < n; ++i) {
    double trials  = std::round(N[i]);
    double success = std::round(x[i] * N[i]);
    double p       = p_vec[i];
    double q       = q_vec[i];
    
    res[i] = -::Rf_dbinom_raw(success, trials, p, q, lg);
    
    // Optional diagnostics
    // if (!std::isfinite(res[i])) {
    //     Rcpp::Rcout << "[WARN] neg_dbinom_glmb_rmat: non-finite res[" << i << "]\n"
    //                 << "  trials  = " << trials  << "\n"
    //                 << "  success = " << success << "\n"
    //                 << "  p       = " << p       << "\n"
    //                 << "  q       = " << q       << "\n"
    //                 << "  raw x   = " << x[i]    << "\n"
    //                 << "  raw N   = " << N[i]    << "\n";
    // }
  }
}




NumericVector dbinom_glmb(NumericVector x, NumericVector N, NumericVector means, int lg) {
  int n = x.size();
  NumericVector res(n);
  
  for (int i = 0; i < n; i++) {
    // Round to nearest integer for trial count and success count
    int trials  = static_cast<int>(std::round(N[i]));
    int success = static_cast<int>(std::round(x[i] * N[i]));
    
    // Clamp probabilities to avoid extreme values
    double p = std::min(1.0, std::max(0.0, means[i]));
    
    res[i] = ::Rf_dbinom(success, trials, p, lg);

    
      }
  
  return res;
}


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
        progress_bar(i, m1-1);
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
    //yy=-dbinom_glmb_raw(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), res1);

    }
    
    return res;      
}



arma::vec f2_binomial_logit_rmat(
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    const int progbar /*=0*/
) {
  std::size_t l1 = x.nrow();
  std::size_t l2 = x.ncol();
  std::size_t m1 = b.ncol();
  
  // Armadillo views over RMatrix memory
  arma::mat b2full(const_cast<double*>(&*b.begin()), l2, m1, false);
  arma::mat x2(const_cast<double*>(&*x.begin()), l1, l2, false);
  arma::mat mu2(const_cast<double*>(&*mu.begin()), l2, 1, false);
  arma::mat P2(const_cast<double*>(&*P.begin()), l2, l2, false);
  arma::mat alpha2(const_cast<double*>(&*alpha.begin()), l1, 1, false);
  
  arma::vec res(m1, arma::fill::none);
  arma::mat bmu(l2, 1, arma::fill::none);
  
  // Buffers (declare once, reuse per i)
  std::vector<double> p_vec(l1), q_vec(l1), yy_temp(l1);
  
  for (std::size_t i = 0; i < m1; ++i) {
    arma::mat b_i(b2full.colptr(i), l2, 1, false);
    
    // Mahalanobis penalty
    bmu = b_i - mu2;
    double mahal = 0.5 * arma::as_scalar(bmu.t() * P2 * bmu);
    
    // Compute p and q stably for each observation
    for (std::size_t j = 0; j < l1; ++j) {
      double eta = alpha2(j,0) + arma::dot(x2.row(j), b_i);
      double p, q;
      if (eta >= 0) {
        double e = std::exp(-eta);
        p = 1.0 / (1.0 + e);
        q = e / (1.0 + e);
      } else {
        double e = std::exp(eta);
        p = e / (1.0 + e);
        q = 1.0 / (1.0 + e);
      }
      p_vec[j] = p;
      q_vec[j] = q;
    }
    
    // Call refactored backend with both p and q
    neg_dbinom_glmb_rmat(y, wt, p_vec, q_vec, yy_temp, /*lg=*/1);
    
    // Accumulate with penalty
    res(i) = std::accumulate(yy_temp.begin(), yy_temp.end(), mahal);
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
    progress_bar(i, m1-1);
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




// Combined f2/f3 for binomial–logit (CPU, non‑OpenCL)
// Matches the existing f2_binomial_logit and f3_binomial_logit numerically,
// but computes both in a single pass over the data.

Rcpp::List f2_f3_binomial_logit(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar
) {
  // Dimensions
  int l1 = x.nrow();   // observations
  int l2 = x.ncol();   // coefficients
  int m1 = b.ncol();   // grid points
  
  // Temporary views / buffers
  Rcpp::NumericMatrix b2temp(l2, 1);
  arma::mat x2    (x.begin(),     l1, l2, false);
  arma::mat alpha2(alpha.begin(), l1, 1,  false);
  
  Rcpp::NumericVector xb(l1);
  arma::colvec xb2(xb.begin(), l1, false);  // shares memory with xb
  
  Rcpp::NumericVector yy(l1);
  Rcpp::NumericVector qf(m1);               // f2 output
  
  Rcpp::NumericMatrix bmu(l2, 1);
  arma::mat mu2 (mu.begin(),  l2, 1, false);
  arma::mat bmu2(bmu.begin(), l2, 1, false);
  arma::mat P2  (P.begin(),   l2, l2, false);
  
  // Gradient buffer: l2 × m1, then transposed to m1 × l2 on return
  Rcpp::NumericMatrix out(l2, m1);
  arma::mat out2(out.begin(), l2, m1, false);
  
  for (int i = 0; i < m1; ++i) {
    Rcpp::checkUserInterrupt();
    if (progbar == 1) {
      progress_bar(i, m1 - 1);
      if (i == m1 - 1) Rcpp::Rcout << "" << std::endl;
    }
    
    // Current grid point b[, i]
    b2temp = b(Range(0, l2 - 1), Range(i, i));
    arma::mat b2(b2temp.begin(), l2, 1, false);
    
    // Prior term: 0.5 * (b - mu)' P (b - mu)
    bmu2 = b2 - mu2;
    double prior = 0.5 * arma::as_scalar(bmu2.t() * P2 * bmu2);
    
    // Linear predictor and logistic transform
    // xb2(j) = exp(-alpha - x b), xb(j) = p_j after the loop below
    xb2 = exp(-alpha2 - x2 * b2);
    
    for (int j = 0; j < l1; ++j) {
      xb(j) = 1.0 / (1.0 + xb(j));  // p_j = logistic(η_j)
    }
    
    // f2: negative binomial log-likelihood (vectorized)
    yy = -dbinom_glmb(y, wt, xb, /*lg=*/1);
    qf(i) = std::accumulate(yy.begin(), yy.end(), prior);
    
    // f3 residuals: (p - y) * wt, reusing xb / xb2 buffer
    for (int j = 0; j < l1; ++j) {
      xb(j) = (xb(j) - y(j)) * wt(j);
    }
    
    // Gradient: P (b - mu) + X' * residual
    Rcpp::NumericMatrix::Column outcol = out(_, i);
    arma::mat outtemp2(outcol.begin(), 1, l2, false);
    outtemp2 = P2 * bmu2 + x2.t() * xb2;
  }
  
  arma::mat grad = trans(out2);  // m1 × l2, matching f3_binomial_logit
  
  return Rcpp::List::create(
    Rcpp::Named("qf")   = qf,    // same as f2_binomial_logit(...)
    Rcpp::Named("grad") = grad   // same as f3_binomial_logit(...)
  );
}



///////////////////////// Probit Functions ///////////////////////////////////////
//
// Vectorized pnorm(xb, 0, 1) / dnorm(xb, 0, 1) on NumericVector below are NOT
// custom glmbayes code: they are Rcpp's element-wise stats API (Rcpp/stats/norm.h),
// which ships with Rcpp and ultimately calls libR (Rf_pnorm5 / Rf_dnorm4) per element.
// Left unchanged for now; scalar parallel paths use ::Rf_pnorm5 / ::Rf_dbinom_raw directly.

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
        progress_bar(i, m1-1);
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


arma::vec f2_binomial_probit_rmat(
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    const int progbar
) {
  std::size_t l1 = x.nrow();
  std::size_t l2 = x.ncol();
  std::size_t m1 = b.ncol();
  
  // Armadillo views over RMatrix memory
  arma::mat b2full(const_cast<double*>(&*b.begin()), l2, m1, false);
  arma::mat x2(const_cast<double*>(&*x.begin()), l1, l2, false);
  arma::mat mu2(const_cast<double*>(&*mu.begin()), l2, 1, false);
  arma::mat P2(const_cast<double*>(&*P.begin()), l2, l2, false);
  arma::mat alpha2(const_cast<double*>(&*alpha.begin()), l1, 1, false);
  
  arma::vec res(m1, arma::fill::none);
  arma::mat bmu(l2, 1, arma::fill::none);
  
  // Buffers
  std::vector<double> p_vec(l1), q_vec(l1), yy_temp(l1);
  
  for (std::size_t i = 0; i < m1; ++i) {
    arma::mat b_i(b2full.colptr(i), l2, 1, false);
    
    bmu = b_i - mu2;
    double mahal = 0.5 * arma::as_scalar(bmu.t() * P2 * bmu);
    
    // Compute linear predictor
    arma::colvec eta = alpha2 + x2 * b_i;
    
    // Compute p and q stably via probit link
    for (std::size_t j = 0; j < l1; j++) {
      double z = eta(j);
      double p = ::Rf_pnorm5(z, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/0);
      double q = ::Rf_pnorm5(z, 0.0, 1.0, /*lower_tail=*/0, /*log_p=*/0);
      p_vec[j] = p;
      q_vec[j] = q;
    }
    
    // Call refactored backend
    neg_dbinom_glmb_rmat(y, wt, p_vec, q_vec, yy_temp, 1);
    
    res(i) = std::accumulate(yy_temp.begin(), yy_temp.end(), mahal);
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
        progress_bar(i, m1-1);
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


// Rcpp::List f2_f3_binomial_probit(
//     Rcpp::NumericMatrix  b,
//     Rcpp::NumericVector  y,
//     Rcpp::NumericMatrix  x,
//     Rcpp::NumericMatrix  mu,
//     Rcpp::NumericMatrix  P,
//     Rcpp::NumericVector  alpha,
//     Rcpp::NumericVector  wt,
//     int                  progbar
// ) {
//   // Dimensions
//   int l1 = x.nrow();
//   int l2 = x.ncol();
//   int m1 = b.ncol();
//   
//   // Outputs (same storage order as original f3: we’ll transpose at the end)
//   Rcpp::NumericVector qf(m1);
//   Rcpp::NumericMatrix grad(l2, m1);
//   arma::mat grad2(grad.begin(), l2, m1, false);
//   
//   // Common Armadillo views
//   arma::mat x2(x.begin(),     l1, l2, false);
//   arma::mat mu2(mu.begin(),   l2, 1,  false);
//   arma::mat P2(P.begin(),     l2, l2, false);
//   arma::mat alpha2(alpha.begin(), l1, 1, false);
//   
//   // Temporaries matching original f2/f3
//   Rcpp::NumericMatrix b2temp(l2, 1);
//   Rcpp::NumericVector xb(l1);
//   arma::colvec xb2(xb.begin(), l1, false);
//   
//   Rcpp::NumericMatrix bmu(l2, 1);
//   arma::mat bmu2(bmu.begin(), l2, 1, false);
//   
//   Rcpp::NumericVector yy(l1);
//   Rcpp::NumericVector p1(l1);
//   Rcpp::NumericVector p2(l1);
//   Rcpp::NumericVector d1(l1);
//   
//   for (int i = 0; i < m1; i++) {
//     Rcpp::checkUserInterrupt();
//     if (progbar == 1) {
//       progress_bar(i, m1 - 1);
//       if (i == m1 - 1) Rcpp::Rcout << "" << std::endl;
//     }
//     
//     // --- extract b_i exactly as in originals ---
//     b2temp = b(Rcpp::Range(0, l2 - 1), Rcpp::Range(i, i));
//     arma::mat b2(b2temp.begin(), l2, 1, false);
//     
//     // --- prior term (shared by f2 and f3) ---
//     bmu2 = b2 - mu2;
//     double mahal = 0.5 * arma::as_scalar(bmu2.t() * P2 * bmu2);
//     
//     // =========================
//     // f2 part (matches f2_binomial_probit)
//     // =========================
//     xb2 = alpha2 + x2 * b2;          // eta
//     xb  = pnorm(xb, 0.0, 1.0);       // p = Φ(eta)
//     
//     yy   = -dbinom_glmb(y, wt, xb, true);
//     qf[i] = std::accumulate(yy.begin(), yy.end(), mahal);
//     
//     // =========================
//     // f3 part (matches f3_binomial_probit)
//     // =========================
//     xb2 = alpha2 + x2 * b2;          // eta again
//     p1  = pnorm(xb, 0.0, 1.0);       // Φ(eta)
//     p2  = pnorm(-xb, 0.0, 1.0);      // Φ(-eta)
//     d1  = dnorm(xb, 0.0, 1.0);       // φ(eta)
//     
//     for (int j = 0; j < l1; j++) {
//       xb(j) = (y(j) * d1(j) / p1(j) - (1 - y(j)) * d1(j) / p2(j)) * wt(j);
//     }
//     
//     // xb2 shares memory with xb, as in original f3
//     // This fills columns in the same order as original out2
//     grad2.col(i) = P2 * bmu2 - x2.t() * xb2;
//   }
//   
//   // Match legacy f3 return orientation (so EnvelopeEval sees 9 x 2, not 2 x 9)
//   arma::mat grad_out = grad2.t();
//   
//   return Rcpp::List::create(
//     Rcpp::Named("qf")   = qf,
//     Rcpp::Named("grad") = Rcpp::wrap(grad_out)
//   );
// }


Rcpp::List f2_f3_binomial_probit(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar
) {
  // Dimensions
  int l1 = x.nrow();
  int l2 = x.ncol();
  int m1 = b.ncol();
  
  // Outputs (same storage order as original f3: we’ll transpose at the end)
  Rcpp::NumericVector qf(m1);
  Rcpp::NumericMatrix grad(l2, m1);
  arma::mat grad2(grad.begin(), l2, m1, false);
  
  // Common Armadillo views
  arma::mat x2(x.begin(),         l1, l2, false);
  arma::mat mu2(mu.begin(),       l2, 1,  false);
  arma::mat P2(P.begin(),         l2, l2, false);
  arma::mat alpha2(alpha.begin(), l1, 1, false);
  
  // Temporaries matching original f2/f3
  Rcpp::NumericMatrix b2temp(l2, 1);
  Rcpp::NumericVector xb(l1);                 // mutable buffer
  arma::colvec xb2(xb.begin(), l1, false);    // view on xb
  
  Rcpp::NumericMatrix bmu(l2, 1);
  arma::mat bmu2(bmu.begin(), l2, 1, false);
  
  Rcpp::NumericVector yy(l1);
  Rcpp::NumericVector p1(l1);
  Rcpp::NumericVector p2(l1);
  Rcpp::NumericVector d1(l1);
  
  // NEW: separate, immutable buffer for eta = alpha + X b
  Rcpp::NumericVector eta(l1);
  arma::colvec eta2(eta.begin(), l1, false);
  
  for (int i = 0; i < m1; i++) {
    Rcpp::checkUserInterrupt();
    if (progbar == 1) {
      progress_bar(i, m1 - 1);
      if (i == m1 - 1) Rcpp::Rcout << "" << std::endl;
    }
    
    // --- extract b_i exactly as in originals ---
    b2temp = b(Rcpp::Range(0, l2 - 1), Rcpp::Range(i, i));
    arma::mat b2(b2temp.begin(), l2, 1, false);
    
    // --- prior term (shared by f2 and f3) ---
    bmu2 = b2 - mu2;
    double mahal = 0.5 * arma::as_scalar(bmu2.t() * P2 * bmu2);
    
    // =====================================================
    // Compute eta once and keep it in a separate buffer
    // =====================================================
    eta2 = alpha2 + x2 * b2;      // eta = alpha + X b
    
    // =========================
    // f2 part (matches f2_binomial_probit)
    // =========================
    for (int j = 0; j < l1; j++) {
      xb[j] = R::pnorm(eta[j], 0.0, 1.0, 1, 0);   // p = Φ(eta)
    }
    
    yy   = -dbinom_glmb(y, wt, xb, true);
    qf[i] = std::accumulate(yy.begin(), yy.end(), mahal);
    
    // =========================
    // f3 part (matches f3_binomial_probit)
    // =========================
    for (int j = 0; j < l1; j++) {
      double eta_j = eta[j];
      double p     = ::Rf_pnorm5(eta_j, 0.0, 1.0, 1, 0);   // Φ(eta)
      double q     = ::Rf_pnorm5(-eta_j, 0.0, 1.0, 1, 0);  // Φ(-eta)
      double d     = ::Rf_dnorm4(eta_j, 0.0, 1.0, 0);      // φ(eta)
      
      p1[j] = p;
      p2[j] = q;
      d1[j] = d;
      
      xb[j] = (y[j] * d / p - (1.0 - y[j]) * d / q) * wt[j];
    }
    
    // xb2 shares memory with xb, as in original f3
    grad2.col(i) = P2 * bmu2 - x2.t() * xb2;
  }
  
  // Match legacy f3 return orientation (so EnvelopeEval sees m1 × l2)
  arma::mat grad_out = grad2.t();
  
  return Rcpp::List::create(
    Rcpp::Named("qf")   = qf,
    Rcpp::Named("grad") = Rcpp::wrap(grad_out)
  );
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
  
    for(int j=0;j<l1;j++){
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
        progress_bar(i, m1-1);
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

arma::vec f2_binomial_cloglog_rmat(
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    const int progbar /*=0*/
) {
  std::size_t l1 = x.nrow();
  std::size_t l2 = x.ncol();
  std::size_t m1 = b.ncol();
  
  // Armadillo views over RMatrix memory
  arma::mat b2full(const_cast<double*>(&*b.begin()), l2, m1, false);
  arma::mat x2(const_cast<double*>(&*x.begin()), l1, l2, false);
  arma::mat mu2(const_cast<double*>(&*mu.begin()), l2, 1, false);
  arma::mat P2(const_cast<double*>(&*P.begin()), l2, l2, false);
  arma::mat alpha2(const_cast<double*>(&*alpha.begin()), l1, 1, false);
  
  arma::vec res(m1, arma::fill::none);
  arma::mat bmu(l2, 1, arma::fill::none);
  
  // Buffers reused per coefficient index i
  std::vector<double> p_vec(l1), q_vec(l1), yy_temp(l1);
  
  for (std::size_t i = 0; i < m1; ++i) {
    arma::mat b_i(b2full.colptr(i), l2, 1, false);
    
    // Mahalanobis penalty
    bmu = b_i - mu2;
    double mahal = 0.5 * arma::as_scalar(bmu.t() * P2 * bmu);
    
    // Linear predictor: cloglog uses η = α + x b
    arma::colvec eta = alpha2 + x2 * b_i;
    
    // p = 1 - exp(-exp(η)), q = exp(-exp(η))
    for (std::size_t j = 0; j < l1; ++j) {
      double t = std::exp(eta(j));              // t = exp(η), may be Inf for large η
      double q = std::exp(-t);                  // q = exp(-t), safe even if t = Inf -> q = 0
      double p;
      // For tiny t, use expm1 to avoid cancellation: 1 - exp(-t) ≈ t
      if (t < 1e-6) {
        p = -std::expm1(-t);                  // p ≈ t for small t
      } else {
        p = 1.0 - q;
      }
      p_vec[j] = p;
      q_vec[j] = q;
    }
    
    // Tail-stable backend with both p and q
    neg_dbinom_glmb_rmat(y, wt, p_vec, q_vec, yy_temp, /*lg=*/1);
    
    // Accumulate with penalty
    res(i) = std::accumulate(yy_temp.begin(), yy_temp.end(), mahal);
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
        progress_bar(i, m1-1);
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

// Rcpp::List f2_f3_binomial_cloglog(
//     Rcpp::NumericMatrix  b,
//     Rcpp::NumericVector  y,
//     Rcpp::NumericMatrix  x,
//     Rcpp::NumericMatrix  mu,
//     Rcpp::NumericMatrix  P,
//     Rcpp::NumericVector  alpha,
//     Rcpp::NumericVector  wt,
//     int                  progbar
// ) {
//   // Dimensions
//   int l1 = x.nrow();
//   int l2 = x.ncol();
//   int m1 = b.ncol();
//   
//   // Outputs (same storage order as original f3: we’ll transpose at the end)
//   Rcpp::NumericVector qf(m1);
//   Rcpp::NumericMatrix grad(l2, m1);
//   arma::mat grad2(grad.begin(), l2, m1, false);
//   
//   // Common Armadillo views
//   arma::mat x2(x.begin(),     l1, l2, false);
//   arma::mat mu2(mu.begin(),   l2, 1,  false);
//   arma::mat P2(P.begin(),     l2, l2, false);
//   arma::mat alpha2(alpha.begin(), l1, 1, false);
//   
//   // Temporaries matching original f2/f3
//   Rcpp::NumericMatrix b2temp(l2, 1);
//   Rcpp::NumericVector xb(l1);
//   arma::colvec xb2(xb.begin(), l1, false);
//   
//   Rcpp::NumericMatrix bmu(l2, 1);
//   arma::mat bmu2(bmu.begin(), l2, 1, false);
//   
//   Rcpp::NumericVector yy(l1);
//   Rcpp::NumericVector p1(l1);
//   Rcpp::NumericVector p2(l1);
//   Rcpp::NumericVector atemp(l1);
//   
//   for (int i = 0; i < m1; i++) {
//     Rcpp::checkUserInterrupt();
//     if (progbar == 1) {
//       progress_bar(i, m1 - 1);
//       if (i == m1 - 1) Rcpp::Rcout << "" << std::endl;
//     }
//     
//     // --- extract b_i exactly as in originals ---
//     b2temp = b(Rcpp::Range(0, l2 - 1), Rcpp::Range(i, i));
//     arma::mat b2(b2temp.begin(), l2, 1, false);
//     
//     // --- prior term (shared by f2 and f3) ---
//     bmu2 = b2 - mu2;
//     double mahal = 0.5 * arma::as_scalar(bmu2.t() * P2 * bmu2);
//     
//     // =========================
//     // f2 part (matches f2_binomial_cloglog)
//     // =========================
//     xb2 = alpha2 + x2 * b2;          // eta
//     
//     for (int j = 0; j < l1; j++) {
//       xb(j) = 1.0 - std::exp(-std::exp(xb(j)));
//     }
//     
//     yy   = -dbinom_glmb(y, wt, xb, true);
//     qf[i] = std::accumulate(yy.begin(), yy.end(), mahal);
//     
//     // =========================
//     // f3 part (matches f3_binomial_cloglog)
//     // =========================
//     xb2 = alpha2 + x2 * b2;          // eta again
//     
//     for (int j = 0; j < l1; j++) {
//       p1(j)    = 1.0 - std::exp(-std::exp(xb(j)));
//       p2(j)    = std::exp(-std::exp(xb(j)));
//       atemp(j) = std::exp(xb(j) - std::exp(xb(j)));
//       xb(j)    = ( (y(j) * atemp(j) / p1(j)) -
//         ((1.0 - y(j)) * atemp(j) / p2(j)) ) * wt(j);
//     }
//     
//     // xb2 shares memory with xb, as in original f3
//     grad2.col(i) = P2 * bmu2 - x2.t() * xb2;
//   }
//   
//   // Match legacy f3 return orientation
//   arma::mat grad_out = grad2.t();
//   
//   return Rcpp::List::create(
//     Rcpp::Named("qf")   = qf,
//     Rcpp::Named("grad") = Rcpp::wrap(grad_out)
//   );
// }

Rcpp::List f2_f3_binomial_cloglog(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar
) {
  // Dimensions
  int l1 = x.nrow();
  int l2 = x.ncol();
  int m1 = b.ncol();
  
  // Outputs (same storage order as original f3: we’ll transpose at the end)
  Rcpp::NumericVector qf(m1);
  Rcpp::NumericMatrix grad(l2, m1);
  arma::mat grad2(grad.begin(), l2, m1, false);
  
  // Common Armadillo views
  arma::mat x2(x.begin(),         l1, l2, false);
  arma::mat mu2(mu.begin(),       l2, 1,  false);
  arma::mat P2(P.begin(),         l2, l2, false);
  arma::mat alpha2(alpha.begin(), l1, 1,  false);
  
  // Temporaries matching original f2/f3
  Rcpp::NumericMatrix b2temp(l2, 1);
  Rcpp::NumericVector xb(l1);                 // mutable buffer
  arma::colvec xb2(xb.begin(), l1, false);    // view on xb
  
  Rcpp::NumericMatrix bmu(l2, 1);
  arma::mat bmu2(bmu.begin(), l2, 1, false);
  
  Rcpp::NumericVector yy(l1);
  Rcpp::NumericVector p1(l1);
  Rcpp::NumericVector p2(l1);
  Rcpp::NumericVector atemp(l1);
  
  // NEW: separate, immutable buffer for eta = alpha + X b
  Rcpp::NumericVector eta(l1);
  arma::colvec eta2(eta.begin(), l1, false);
  
  for (int i = 0; i < m1; i++) {
    Rcpp::checkUserInterrupt();
    if (progbar == 1) {
      progress_bar(i, m1 - 1);
      if (i == m1 - 1) Rcpp::Rcout << "" << std::endl;
    }
    
    // --- extract b_i exactly as in originals ---
    b2temp = b(Rcpp::Range(0, l2 - 1), Rcpp::Range(i, i));
    arma::mat b2(b2temp.begin(), l2, 1, false);
    
    // --- prior term (shared by f2 and f3) ---
    bmu2 = b2 - mu2;
    double mahal = 0.5 * arma::as_scalar(bmu2.t() * P2 * bmu2);
    
    // =====================================================
    // Compute eta once and keep it in a separate buffer
    // =====================================================
    eta2 = alpha2 + x2 * b2;      // eta = alpha + X b
    
    // =========================
    // f2 part (matches f2_binomial_cloglog)
    // =========================
    for (int j = 0; j < l1; j++) {
      double e_eta = std::exp(eta[j]);
      xb[j] = 1.0 - std::exp(-e_eta);   // p = 1 - exp(-exp(eta))
    }
    
    yy   = -dbinom_glmb(y, wt, xb, true);
    qf[i] = std::accumulate(yy.begin(), yy.end(), mahal);
    
    // =========================
    // f3 part (matches f3_binomial_cloglog)
    // =========================
    for (int j = 0; j < l1; j++) {
      double e_eta = std::exp(eta[j]);                 // exp(eta)
      double p     = 1.0 - std::exp(-e_eta);           // p1
      double q     = std::exp(-e_eta);                 // p2
      double a     = std::exp(eta[j] - e_eta);         // atemp
      
      p1[j]    = p;
      p2[j]    = q;
      atemp[j] = a;
      
      xb[j] = ( (y[j] * a / p) - ((1.0 - y[j]) * a / q) ) * wt[j];
    }
    
    // xb2 shares memory with xb, as in original f3
    grad2.col(i) = P2 * bmu2 - x2.t() * xb2;
  }
  
  // Match legacy f3 return orientation
  arma::mat grad_out = grad2.t();
  
  return Rcpp::List::create(
    Rcpp::Named("qf")   = qf,
    Rcpp::Named("grad") = Rcpp::wrap(grad_out)
  );
}


}  //famfuncs
}  //glmbayes
