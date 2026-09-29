// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppParallel.h>
#include <Rmath.h>   // libR Mathlib: Rf_dgamma
#include "famfuncs.h"
#include "progress_utils.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace glmbayes::fam;
using namespace glmbayes::progress;




namespace glmbayes{

namespace fam {
void neg_dgamma_glmb_rmat(const RVector<double>& x,           // observations
                          const RVector<double>& shape,       // shape parameters
//                          const RVector<double>& scale,       // scale parameters
                          const std::vector<double>& scale, // scale parameters
                          std::vector<double>& res,         // output buffer
                          const int lg)                       // log=TRUE?
{
  std::size_t n = x.length();
  
  for (std::size_t i = 0; i < n; ++i) {
    double value = x[i];
    double k     = shape[i];
    double theta = scale[i];
    
    res[i] = -::Rf_dgamma(value, k, theta, lg);

      }
}



NumericVector dgamma_glmb( NumericVector x, NumericVector shape, NumericVector scale, int lg){
    int n = x.size() ;
    NumericVector res(n) ;
  for( int i=0; i<n; i++) res[i] = ::Rf_dgamma( x[i], shape[i], scale[i], lg );
    return res ;
}


////////////////////////////////////////////////////////////////
// See if it is possible to avoid having some or all of these functions exported


NumericVector  f1_gamma(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt)
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
    arma::mat wt2(wt.begin(), l1, 1, false);
    
    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);


    for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
  

//  mu<-t(exp(alpha+x%*%b))
//  disp2<-1/wt

//  -sum(dgamma(y,shape=1/disp2,scale=mu*disp2,log=TRUE))


    xb2=exp(alpha2+ x2 * b2);
    
    for(int j=0;j<l1;j++){
      
    xb[j]=xb[j]/wt[j];  
    }

    yy=-dgamma_glmb(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);

    }
    
    return res;      
}




NumericVector  f2_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
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
    arma::mat wt2(wt.begin(), l1, 1, false);

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     

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
    
    xb2=exp(alpha2+ x2 * b2);
    
    for(int j=0;j<l1;j++){
      
    xb[j]=xb[j]/wt[j];  
    }

    yy=-dgamma_glmb(y,wt,xb,true);


    res(i) =std::accumulate(yy.begin(), yy.end(), res1);

    }
    
    return res;      
}



arma::vec f2_gamma_rmat(
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
    
    
    //xb_temp2 = alpha2+  x2 * b_i;
    
//    for (std::size_t j = 0; j < l1; j++) {
//      xb_temp[j] =1-  exp(-exp(xb_temp[j]));
//    }
    
    xb_temp2=exp(alpha2+ x2 * b_i);
    
    for (std::size_t  j = 0; j < l1; j++) {      
      xb_temp[j]=xb_temp[j]/wt[j];  
    }
    
    
    //    for(int j=0;j<l1;j++){
    //      xb(j)=1-exp(-exp(xb(j)));
    //    }
    
    
    // In-place evaluation using your log-scale accurate backend
    neg_dgamma_glmb_rmat(y, wt, xb_temp, yy_temp,1.0);
    
    
    res(i) =std::accumulate(yy_temp.begin(), yy_temp.end(), mahal);
    
  }
  
  return res;
}







arma::mat  f3_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
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
    arma::mat wt2(wt.begin(), l1, 1, false);

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
       
// Note: Seem not to be used - Editing out
//    NumericVector invwt=1/wt;

       
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

//  		mu2<-t(exp(alpha+x%*%b))
//		t(x)%*%(t(1-y/mu2)*wt)+P%*%(b-mu)


    bmu2=b2-mu2;
    
    
    xb2=exp(alpha2+ x2 * b2);
    
    for(int j=0;j<l1;j++){
      xb[j]=(1-y[j]/xb[j])*wt[j];
      
    }

    outtemp2= P2 * bmu2+x2.t() * xb2;
    }
    
   // return  b;
  
    return trans(out2);      
}

// Rcpp::List f2_f3_gamma(
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
//   // Outputs (same layout as original f3: we transpose at the end)
//   Rcpp::NumericVector qf(m1);
//   Rcpp::NumericMatrix grad(l2, m1);
//   arma::mat grad2(grad.begin(), l2, m1, false);
//   
//   // Armadillo views
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
//   Rcpp::NumericVector p1(l1);     // not used but kept for structural symmetry
//   Rcpp::NumericVector p2(l1);     // not used but kept for structural symmetry
//   Rcpp::NumericVector atemp(l1);  // not used but kept for structural symmetry
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
//     // f2 part (matches f2_gamma)
//     // =========================
//     xb2 = arma::exp(alpha2 + x2 * b2);   // xb = exp(alpha + X b)
//     
//     for (int j = 0; j < l1; j++) {
//       xb[j] = xb[j] / wt[j];            // EXACTLY as in original f2_gamma
//     }
//     
//     yy = -dgamma_glmb(y, wt, xb, true);
//     
//     qf[i] = std::accumulate(yy.begin(), yy.end(), mahal);
//     
//     // =========================
//     // f3 part (matches f3_gamma)
//     // =========================
//     xb2 = arma::exp(alpha2 + x2 * b2);  // mu = exp(alpha + Xb)
//     
//     for (int j = 0; j < l1; j++) {
//       xb[j] = (1.0 - y[j] / xb[j]) * wt[j];   // EXACT gradient term
//     }
//     
//     // xb2 shares memory with xb, as in original f3
//     grad2.col(i) = P2 * bmu2 + x2.t() * xb2;
//   }
//   
//   // Match legacy f3 return orientation (trans(out2))
//   arma::mat grad_out = grad2.t();
//   
//   return Rcpp::List::create(
//     Rcpp::Named("qf")   = qf,
//     Rcpp::Named("grad") = Rcpp::wrap(grad_out)
//   );
// }


Rcpp::List f2_f3_gamma(
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
  
  // Outputs (same layout as original f3: we transpose at the end)
  Rcpp::NumericVector qf(m1);
  Rcpp::NumericMatrix grad(l2, m1);
  arma::mat grad2(grad.begin(), l2, m1, false);
  
  // Armadillo views
  arma::mat x2(x.begin(),          l1, l2, false);
  arma::mat mu2(mu.begin(),        l2, 1,  false);
  arma::mat P2(P.begin(),          l2, l2, false);
  arma::mat alpha2(alpha.begin(),  l1, 1,  false);
  
  // Temporaries matching original f2/f3
  Rcpp::NumericMatrix b2temp(l2, 1);
  Rcpp::NumericVector xb(l1);                 // mutable buffer
  arma::colvec xb2(xb.begin(), l1, false);    // view on xb
  
  Rcpp::NumericMatrix bmu(l2, 1);
  arma::mat bmu2(bmu.begin(), l2, 1, false);
  
  Rcpp::NumericVector yy(l1);
  
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
    // f2 part (matches f2_gamma)
    // =========================
    // mu = exp(eta); we copy into xb (mutable) for f2
    for (int j = 0; j < l1; j++) {
      xb[j] = std::exp(eta[j]) / wt[j];   // EXACTLY as original: xb = exp(eta)/wt
    }
    
    yy = -dgamma_glmb(y, wt, xb, true);
    qf[i] = std::accumulate(yy.begin(), yy.end(), mahal);
    
    // =========================
    // f3 part (matches f3_gamma)
    // =========================
    // Reuse eta to get mu again, without recomputing X*b
    for (int j = 0; j < l1; j++) {
      double mu_j = std::exp(eta[j]);                 // mu = exp(eta)
      xb[j] = (1.0 - y[j] / mu_j) * wt[j];            // EXACT gradient term
    }
    
    // xb2 shares memory with xb, as in original f3
    grad2.col(i) = P2 * bmu2 + x2.t() * xb2;
  }
  
  // Match legacy f3 return orientation (trans(out2))
  arma::mat grad_out = grad2.t();
  
  return Rcpp::List::create(
    Rcpp::Named("qf")   = qf,
    Rcpp::Named("grad") = Rcpp::wrap(grad_out)
  );
}

} //famfuncs

} //glmbayes

///////////////////////////////////////////////////////////////////