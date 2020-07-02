#include "RcppArmadillo.h"
#include <Rcpp.h>
using namespace Rcpp;


NumericVector dnorm_glmb( NumericVector x, NumericVector means, NumericVector sds,int lg)
{
  
  
  int n = x.size() ;
  NumericVector res(n) ;
  for( int i=0; i<n; i++) res[i] = R::dnorm( x[i], means[i], sds[i],lg ) ;

  return res ;
}


NumericVector RSS(NumericVector y, NumericMatrix x,NumericMatrix b,NumericVector alpha,NumericVector wt)
{
  // Step 1: Set up dimensions
  
  int l1 = x.nrow(), l2 = x.ncol(); // Dimensions of x matrix (dims for y,alpha, and wt needs to be consistent) 
  int m1 = b.ncol();                // Number of columns for which output is needed
  
  // Step 2: Initialize b2temp and other Rcpp and arma objects used in calculations
  
  Rcpp::NumericMatrix b2temp(l2,1);
  Rcpp::NumericMatrix restemp(1,1);
  arma::mat y2(y.begin(), l1, 1, false);
  arma::mat x2(x.begin(), l1, l2, false); 
  arma::mat alpha2(alpha.begin(), l1, 1, false); 
  
  Rcpp::NumericVector xb(l1);
  arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below

  NumericVector sqrt_wt=sqrt(wt);
  arma::mat sqrt_wt2(sqrt_wt.begin(), l1, 1, false); 
  
  //  NumericVector invwt=1/sqrt(wt);
  
  // Moving Loop inside the function is key for speed
  
  NumericVector yy(l1);
  NumericVector res(m1);
  arma::colvec res2(res.begin(),m1,false); // Reuse memory - update both below
  
  for(int i=0;i<m1;i++){
    
  // Grab one column at a time from b and one row at a time from res
  
    b2temp=b(Range(0,l2-1),Range(i,i));

  // Point b2 to memory for that column
  
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat restemp(res.begin()+i, 1, 1, false); 
    
  // calculate weighted residuals (element by element multiplication with weights)
  
    xb2=(y2-alpha2- x2 * b2)%sqrt_wt2;
  
  // This is where RSS should be calculated
  // Not sure if this will complain about type differences
    
    restemp=trans(xb2)*xb2;

  }
  
  return res;      

  }



NumericVector  f1_gaussian(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt)
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
  
  
  NumericVector invwt=1/sqrt(wt);
  
  // Moving Loop inside the function is key for speed
  
  NumericVector yy(l1);
  NumericVector res(m1);
  
  
  for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    
    
    xb2=alpha2+ x2 * b2;
    
    yy=-dnorm_glmb(y,xb,invwt,true);
    
    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);
    
  }
  
  return res;      
}


// [[Rcpp::export(".f2_gaussian_vector")]]


NumericVector  f2_gaussian(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt)
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
  
  
  NumericVector invwt=1/sqrt(wt);
  
  // Moving Loop inside the function is key for speed
  
  NumericVector yy(l1);
  NumericVector res(m1);
  NumericMatrix bmu(l2,1);
  
  arma::mat mu2(mu.begin(), l2, 1, false); 
  arma::mat bmu2(bmu.begin(), l2, 1, false); 


  double res1=0;

  
  for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 
    
    bmu2=b2-mu2;
    
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
    xb2=alpha2+ x2 * b2;
    
    yy=-dnorm_glmb(y,xb,invwt,true);
  
    res(i) =std::accumulate(yy.begin(), yy.end(), res1);
    
  }
  
  return res;      
}




arma::mat  f3_gaussian(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt)
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
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);
    
    bmu2=b2-mu2;
    xb2=alpha2+ x2 * b2-y2;
    outtemp2= P2 * bmu2+x2.t() * Ptemp2 * xb2;
  }
  
  // return  b;
  
  return trans(out2);      
}

// [[Rcpp::export(".Inv_f3_gaussian")]]


arma::mat  Inv_f3_gaussian(NumericMatrix cbars,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt)
{
  
  // Get dimensions of x - Note: should match dimensions of
  //  y, b, alpha, and wt (may add error checking)
  
  // May want to add method for dealing with alpha and wt when 
  // constants instead of vectors
  
  int l1 = x.nrow(), l2 = x.ncol();
  int m1 = cbars.ncol();  // Check if this makes sense or if it is better to pass the rows
  
  //    int lalpha=alpha.nrow();
  //    int lwt=wt.nrow();
  
  Rcpp::NumericMatrix cbars2temp(l2,1);
  
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
  NumericMatrix out(l2,m1);
  
  arma::mat cbars2(cbars.begin(), l2, 1, false); 
  
  arma::mat mu2(mu.begin(), l2, 1, false); 
  arma::mat P2(P.begin(), l2, l2, false); 
  arma::mat Ptemp2(Ptemp.begin(), l1, l1, false);
  arma::mat out2(out.begin(), l2, m1, false);
  
  NumericMatrix::Column outtemp=out(_,0);

  arma::mat outtemp2(outtemp.begin(),1,l2,false);

  for(int i=0;i<m1;i++){
    cbars2temp=cbars(Range(0,l2-1),Range(i,i));
    arma::mat cbars2(cbars2temp.begin(), l2, 1, false); 
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);

    xb2=alpha2 -y2;
    
    outtemp2=-inv_sympd(P2+x2.t() * Ptemp2*x2)*(-cbars2+x2.t() * Ptemp2 *xb2+P2*mu2);
    
  }
  

  return trans(out2);      
}
