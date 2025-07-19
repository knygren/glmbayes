// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

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
                          const RVector<double>& p_vec,     // standard-scale probabilities
                          RVector<double>& res,             // output buffer
                          const int lg) {
  std::size_t n = x.length();
  
  for (std::size_t i = 0; i < n; ++i) {
    // Input transformation: round values like dbinom_glmb()
    double trials  = std::round(N[i]);
    double success = std::round(x[i] * N[i]);
    
    double p = p_vec[i];         // success probability
    double q = 1.0 - p;          // failure probability
  
//  Rcpp::Rcout << "success: " << success <<  "\n";
//  Rcpp::Rcout << "trials: " << trials <<  "\n";
//  Rcpp::Rcout << "p: " << p <<  "\n";
//  Rcpp::Rcout << "q: " << q <<  "\n";
  
    
    // Evaluate negative log-likelihood using R-accurate approximation
//    res[i] = -dbinom_raw_log(success, trials, p, q);
    
//         Rcpp::Rcout << "c++: " << res[i] <<  "\n";
    
    res[i] = - dbinom_raw(success, trials, p, q, 1);  // log=TRUE

    //Rcpp::Rcout << "dbinom_raw: " << res[i] <<  "\n";
    
    //res[i] = - dbinom_raw_local(success, trials, p, q, 1);  // log=TRUE
    
//    Rcpp::Rcout << "dbinom_raw_local: " << res[i] <<  "\n";
    
  //      res[i] = - R::dbinom(success, trials, p, true);
  //  Rcpp::Rcout << "R: " << res[i] <<  "\n";
    
      }
}



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

arma::vec f2_binomial_logit_arma(NumericMatrix b, NumericVector y,
                                NumericMatrix x, NumericMatrix mu,
                                NumericMatrix P, NumericVector alpha,
                                NumericVector wt, int progbar = 0) {
  int l1 = x.nrow(), l2 = x.ncol();
  int m1 = b.ncol();
  
  arma::mat b2full(b.begin(), l2, m1, false);  // Shallow view over b
  
  Rcpp::NumericMatrix b2temp(l2, 1);
  arma::mat x2(x.begin(), l1, l2, false);
  arma::mat alpha2(alpha.begin(), l1, 1, false);
  
  Rcpp::NumericVector xb(l1);
  arma::colvec xb2(xb.begin(), l1, false);
  
  NumericVector yy(l1);
//  NumericVector res(m1);
  NumericMatrix bmu(l2, 1);
  arma::mat mu2(mu.begin(), l2, 1, false);
  arma::mat bmu2(bmu.begin(), l2, 1, false);

  arma::vec res2(m1);  // Owning allocation
  
  
  for (int i = 0; i < m1; i++) {
    Rcpp::checkUserInterrupt();
    if (progbar == 1) {
      progress_bar2(i, m1 - 1);
      if (i == m1 - 1) { Rcpp::Rcout << "" << std::endl; }
    }
    
    b2temp = b(Range(0, l2 - 1), Range(i, i));
    //arma::mat b2(b2temp.begin(), l2, 1, false);
    arma::mat b2(b2full.colptr(i), l2, 1, false);  // View of column i, size l2 × 1
    
    arma::mat P2(P.begin(), l2, l2, false);
    
    bmu2 = b2 - mu2;
    double res1 = 0.5 * arma::as_scalar(bmu2.t() * P2 * bmu2);
    
    xb2 = exp(-alpha2 - x2 * b2);  // eta = -alpha - Xb
 
    for (int j = 0; j < l1; j++) {
      xb(j) = 1.0 / (1.0 + xb(j));  // logistic link
    }

    
    // Wrap existing R memory buffers
    RcppParallel::RVector<double> y_view(y);
    RcppParallel::RVector<double> wt_view(wt);
    RcppParallel::RVector<double> xb_view(xb);
    RcppParallel::RVector<double> yy_view(yy); //must be preallocated to match y.size()
    
    // In-place evaluation using your log-scale accurate backend
    neg_dbinom_glmb_rmat(y_view, wt_view, xb_view, yy_view,1.0);
        
//    yy=-dbinom_glmb(y,wt,xb,true);
    
    
    //res(i) =std::accumulate(yy.begin(), yy.end(), res1);
    res2(i) = std::accumulate(yy.begin(), yy.end(), res1);

      }
  
  return res2;
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



arma::vec  f2_binomial_probit_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
  
  // Get dimensions of x - Note: should match dimensions of
  //  y, b, alpha, and wt (may add error checking)
  
  // May want to add method for dealing with alpha and wt when 
  // constants instead of vectors
  
  int l1 = x.nrow(), l2 = x.ncol();
  int m1 = b.ncol();
  
  //    int lalpha=alpha.nrow();
  //    int lwt=wt.nrow();
  
  arma::mat b2full(b.begin(), l2, m1, false);  // Shallow view over b
  Rcpp::NumericMatrix b2temp(l2,1);
  
  arma::mat x2(x.begin(), l1, l2, false); 
  arma::mat alpha2(alpha.begin(), l1, 1, false); 
  
  Rcpp::NumericVector xb(l1);
  arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
  
  //   Note: Does not seem to be used-- Editing out        
  //    NumericVector invwt=1/sqrt(wt);
  
  // Moving Loop inside the function is key for speed
  
  NumericVector yy(l1);
//  NumericVector res(m1);
  NumericMatrix bmu(l2,1);
  
  arma::mat mu2(mu.begin(), l2, 1, false); 
  arma::mat bmu2(bmu.begin(), l2, 1, false); 
  arma::vec res2(m1);  // Owning allocation
  
  double res1=0;
  
//  Rcpp::Rcout << "m1=" << m1 <<  "\n";
  
    
  for(int i=0;i<m1;i++){
    Rcpp::checkUserInterrupt();
    if(progbar==1){ 
      progress_bar2(i, m1-1);
      if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
    };  
    
    
    
    b2temp=b(Range(0,l2-1),Range(i,i));
    //arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat b2(b2full.colptr(i), l2, 1, false);  // View of column i, size l2 × 1
    
    arma::mat P2(P.begin(), l2, l2, false); 
    
    bmu2=b2-mu2;
    
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
    
    xb2=alpha2+ x2 * b2;   
    xb=pnorm(xb,0.0,1.0);

    // Wrap existing R memory buffers
    RcppParallel::RVector<double> y_view(y);
    RcppParallel::RVector<double> wt_view(wt);
    RcppParallel::RVector<double> xb_view(xb);
    RcppParallel::RVector<double> yy_view(yy); //must be preallocated to match y.size()
    
    // In-place evaluation using your log-scale accurate backend
    
//    Rcpp::Rcout << " Entering neg_dbin_glmb_rmat  " << "\n";


//  Rcout << "y: " << y << std::endl;
//  Rcout << "wt: " << wt << std::endl;
//  Rcout << "xb: " << xb << std::endl;
  
    neg_dbinom_glmb_rmat(y_view, wt_view, xb_view, yy_view,1.0);

//    Rcout << "yy: " << yy << std::endl;
    
        
//    Rcpp::Rcout << " Exiting neg_dbin_glmb_rmat  " << "\n";
    
    //yy=-dbinom_glmb(y,wt,xb,true);
    
    
    //res(i) =std::accumulate(yy.begin(), yy.end(), res1);
    res2(i) = std::accumulate(yy.begin(), yy.end(), res1);
    
//    Rcpp::Rcout << "i=" << i << " draws=" << draws[i] << "\n";
  
//    Rcpp::Rcout << "res2=" << res2(i) <<  "\n";
  
  }

  //    Rcpp::Rcout << " Returning res2  " << "\n";
  //  Rcpp::Rcout << "m1=" << m1 <<  "\n";
  
  return res2;      
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




arma::vec  f2_binomial_cloglog_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
  
  // Get dimensions of x - Note: should match dimensions of
  //  y, b, alpha, and wt (may add error checking)
  
  // May want to add method for dealing with alpha and wt when 
  // constants instead of vectors
  
  int l1 = x.nrow(), l2 = x.ncol();
  int m1 = b.ncol();
  
  //    int lalpha=alpha.nrow();
  //    int lwt=wt.nrow();
  
  arma::mat b2full(b.begin(), l2, m1, false);  // Shallow view over b
  Rcpp::NumericMatrix b2temp(l2,1);
  
  arma::mat x2(x.begin(), l1, l2, false); 
  arma::mat alpha2(alpha.begin(), l1, 1, false); 
  
  Rcpp::NumericVector xb(l1);
  arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
  
  //   Note: Does not seem to be used-- Editing out      
  //    NumericVector invwt=1/sqrt(wt);
  
  // Moving Loop inside the function is key for speed
  
  NumericVector yy(l1);
//  NumericVector res(m1);
  NumericMatrix bmu(l2,1);
  
  arma::mat mu2(mu.begin(), l2, 1, false); 
  arma::mat bmu2(bmu.begin(), l2, 1, false); 
  arma::vec res2(m1);  // Owning allocation
  
  double res1=0;
  
  
  for(int i=0;i<m1;i++){
    Rcpp::checkUserInterrupt();
    
    if(progbar==1){ 
      progress_bar2(i, m1-1);
      if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
    };  
    
    
    b2temp=b(Range(0,l2-1),Range(i,i));
    //arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat b2(b2full.colptr(i), l2, 1, false);  // View of column i, size l2 × 1
    
    arma::mat P2(P.begin(), l2, l2, false); 
    
    bmu2=b2-mu2;
    
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
    xb2=alpha2+ x2 * b2;   
    
    for(int j=0;j<l1;j++){
      xb(j)=1-exp(-exp(xb(j)));
    }
    
    
    
    // Wrap existing R memory buffers
    RcppParallel::RVector<double> y_view(y);
    RcppParallel::RVector<double> wt_view(wt);
    RcppParallel::RVector<double> xb_view(xb);
    RcppParallel::RVector<double> yy_view(yy); //must be preallocated to match y.size()
    
    // In-place evaluation using your log-scale accurate backend
    neg_dbinom_glmb_rmat(y_view, wt_view, xb_view, yy_view,1.0);
    
    
    //yy=-dbinom_glmb(y,wt,xb,true);
    
    
    //res(i) =std::accumulate(yy.begin(), yy.end(), res1);
    res2(i) = std::accumulate(yy.begin(), yy.end(), res1);
    
  }
  
  return res2;      
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


