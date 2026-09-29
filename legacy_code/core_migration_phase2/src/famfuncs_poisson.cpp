// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppParallel.h>
#include "famfuncs.h"
#include "progress_utils.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace glmbayes::fam;
using namespace glmbayes::progress;



// f1 is negative log-likelihood
// f2 is negative log-posterior
// f3 is gradient for log-posterior


// Log Poisson PMF via lgamma(x+1): allows non-integer x; not Rf_dpois (which forceints counts).
double dpois2(double x,double lambda,int lg){
  
  //test=max(abs(round(x)-x))
  
  //if(test>0){
  //  warning("Non-Integer Values to Poisson Density - Switching to Gamma Function to Evaluate Factorial")
    return(-lambda+x*log(lambda)-lgamma(x+1));
  
  //} 
  
  //return(dpois(x,lambda,log=TRUE))
}


namespace glmbayes{

namespace fam {


void neg_dpois_glmb_rmat(const RVector<double>& x,     // observed counts
                       const std::vector<double>& means, // Poisson rates
                       std::vector<double>& res,         // output buffer (preallocated)
                       const int lg)                     // log=TRUE?
{
  std::size_t n = x.size();
  if (res.size() != n)
    res.resize(n);  // optional: ensure res is sized correctly
  
  for (std::size_t i = 0; i < n; ++i) {
    double count  = std::round(x[i]);     // match integer behavior
    double lambda = means[i];             // rate parameter
    
    res[i] = -dpois2(count, lambda, lg);  // thread-safe Poisson backend
  }
}


NumericVector dpois_glmb( NumericVector x, NumericVector means, int lg){
    int n = x.size() ;
    NumericVector res(n) ;

//    for( int i=0; i<n; i++) res[i] = R::dpois( x[i], means[i], lg ) ;
    for( int i=0; i<n; i++) res[i] = dpois2( x[i], means[i], lg ) ;
    return res ;
}

NumericVector  f1_poisson(NumericMatrix b,NumericVector y,NumericMatrix x,NumericMatrix alpha,NumericVector wt)
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
  
    
    xb2=exp(alpha2+ x2 * b2);

    yy=-dpois_glmb(y,xb,true);
    
    for(int j=0;j<l1;j++){
    yy[j]=yy[j]*wt[j];  
    }

    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);

    }
    
    return res;      
}





NumericVector  f2_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
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
     
// Note: Seem not to be used - Editing out
//    NumericVector invwt=1/sqrt(wt);

    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector yy_alt(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);

    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 

    double res1=0;


    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      

      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 

    bmu2=b2-mu2;
        
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
        xb2=exp(alpha2+ x2 * b2);

        yy=-dpois_glmb(y,xb,true);


    for(int j=0;j<l1;j++){
    yy[j]=yy[j]*wt[j];  
    }


    
    res(i) =std::accumulate(yy.begin(), yy.end(), res1);
    

    }
    
    return res;      
}




// Thread-safe Poisson likelihood using fully wrapped views

arma::vec f2_poisson_rmat(const RMatrix<double>& b,
                             const RVector<double>& y,
                             const RMatrix<double>& x,
                             const RMatrix<double>& mu,
                             const RMatrix<double>& P,
                             const RVector<double>& alpha,
                             const RVector<double>& wt,
                             const int progbar=0)
{
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
  
    

//  Rcpp::NumericVector xb_vec(l1);   // temp buffer for mean
//  Rcpp::NumericVector yy_vec(l1);   // temp buffer for log-likelihood
  
//  RVector<double> xb(xb_vec);       // thread-safe wrapper
//  RVector<double> yy(yy_vec);       // thread-safe wrapper
  
//  arma::colvec xb2(xb.begin(), l1, false);  // view for exp(alpha + x * b)
  
  
 //       Rcpp::Rcout << "i=" << i  << "\n";
  

  for (std::size_t i = 0; i < m1; ++i) {
    arma::mat b_i(b2full.colptr(i), l2, 1, false);
    
    bmu = b_i - mu2;
    
//    Rcpp::Rcout << "b_i=" << b_i  << "\n";
//    Rcpp::Rcout << "mu2=" << mu2  << "\n";
//    Rcpp::Rcout << "bmu=" << bmu  << "\n";
    
    double mahal = 0.5 * arma::as_scalar(bmu.t() * P2 * bmu);

//    Rcpp::Rcout << "mahal=" << mahal  << "\n";
    
    // Compute exp(alpha + X * b)
    
    xb_temp2 = arma::exp(alpha2 + x2 * b_i);
//    xb2 = arma::exp(alpha2 + x2 * b_i);

//    Rcpp::Rcout << "xb_temp2=" << xb_temp2  << "\n";

    
    // Thread-safe Poisson log-likelihood
    neg_dpois_glmb_rmat(y, xb_temp, yy_temp, 1);
    
//    Rcpp::Rcout << "yy_temp: ";
//    for (std::size_t j = 0; j < yy_temp.size(); ++j) {
//      Rcpp::Rcout << yy_temp[j] << " ";
//    }
//    Rcpp::Rcout << "\n";
    
        // Evaluate density using raw RVector views
//    neg_dpois_glmb_rmat(y, xb, yy, 1);
    
    for (std::size_t j = 0; j < l1; ++j) {
//      yy[j] *= wt[j];
      yy_temp[j] *= wt[j];
    }
    
    res(i) = std::accumulate(yy_temp.begin(), yy_temp.end(), mahal);
//    res(i) = std::accumulate(yy.begin(), yy.end(), mahal);
  }
  
  return res;
}


////////////////////////////////////////////////////////////////////////////////////



arma::mat  f3_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
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
    xb2=alpha2+ x2 * b2;
//    xb2=y2-exp(alpha2+ x2 * b2);

    for(int j=0;j<l1;j++){
    xb(j)=(y(j)-exp(xb(j)))*wt(j);  
    }


//        -t(x)%*%((y-exp(alpha+x%*%b))*wt)+P%*%(b-mu)

    outtemp2= P2 * bmu2-x2.t() * xb2;
    }
    
   // return  b;
  
    return trans(out2);      
}

// Rcpp::List f2_f3_poisson(
//     Rcpp::NumericMatrix  b,
//     Rcpp::NumericVector  y,
//     Rcpp::NumericMatrix  x,
//     Rcpp::NumericMatrix  mu,
//     Rcpp::NumericMatrix  P,
//     Rcpp::NumericVector  alpha,
//     Rcpp::NumericVector  wt,
//     int                  progbar
// ) {
//   // Dimensions (match original f2/f3)
//   int l1 = x.nrow();
//   int l2 = x.ncol();
//   int m1 = b.ncol();
//   
//   // Outputs (same layout as original f3: we’ll transpose at the end)
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
//     // f2 part (matches f2_poisson)
//     // =========================
//     xb2 = arma::exp(alpha2 + x2 * b2);   // xb = exp(alpha + X b)
//     
//     yy = -dpois_glmb(y, xb, true);
//     
//     for (int j = 0; j < l1; j++) {
//       yy[j] *= wt[j];
//     }
//     
//     qf[i] = std::accumulate(yy.begin(), yy.end(), mahal);
//     
//     // =========================
//     // f3 part (matches f3_poisson)
//     // =========================
//     xb2 = alpha2 + x2 * b2;              // eta
//     
//     for (int j = 0; j < l1; j++) {
//       xb(j) = (y(j) - std::exp(xb(j))) * wt(j);
//     }
//     
//     // xb2 shares memory with xb, as in original f3
//     grad2.col(i) = P2 * bmu2 - x2.t() * xb2;
//   }
//   
//   // Match legacy f3_poisson return orientation (trans(out2))
//   arma::mat grad_out = grad2.t();
//   
//   return Rcpp::List::create(
//     Rcpp::Named("qf")   = qf,
//     Rcpp::Named("grad") = Rcpp::wrap(grad_out)
//   );
// }

Rcpp::List f2_f3_poisson(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar
) {
  // Dimensions (match original f2/f3)
  int l1 = x.nrow();
  int l2 = x.ncol();
  int m1 = b.ncol();
  
  // Outputs (same layout as original f3: we’ll transpose at the end)
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
    // f2 part (matches f2_poisson)
    // =========================
    // mu = exp(eta); copy into xb (mutable) for f2
    for (int j = 0; j < l1; j++) {
      xb[j] = std::exp(eta[j]);   // xb = mu
    }
    
    yy = -dpois_glmb(y, xb, true);
    
    for (int j = 0; j < l1; j++) {
      yy[j] *= wt[j];             // EXACTLY as in original f2_poisson
    }
    
    qf[i] = std::accumulate(yy.begin(), yy.end(), mahal);
    
    // =========================
    // f3 part (matches f3_poisson)
    // =========================
    // Reuse eta to get mu again, without recomputing X*b
    for (int j = 0; j < l1; j++) {
      double mu_j = std::exp(eta[j]);              // mu = exp(eta)
      xb[j] = (y[j] - mu_j) * wt[j];               // EXACT gradient term
    }
    
    // xb2 shares memory with xb, as in original f3
    grad2.col(i) = P2 * bmu2 - x2.t() * xb2;
  }
  
  // Match legacy f3_poisson return orientation (trans(out2))
  arma::mat grad_out = grad2.t();
  
  return Rcpp::List::create(
    Rcpp::Named("qf")   = qf,
    Rcpp::Named("grad") = Rcpp::wrap(grad_out)
  );
}


} //famfuncs

} //glmbayes