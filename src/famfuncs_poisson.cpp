// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;
void progress_bar2(double x, double N);


// f1 is negative log-likelihood
// f2 is negative log-posterior
// f3 is gradient for log-posterior


double dpois2(double x,double lambda,int lg){
  
  //test=max(abs(round(x)-x))
  
  //if(test>0){
  //  warning("Non-Integer Values to Poisson Density - Switching to Gamma Function to Evaluate Factorial")
    return(-lambda+x*log(lambda)-lgamma(x+1));
  
  //} 
  
  //return(dpois(x,lambda,log=TRUE))
}




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

// void neg_dpois_glmb_rmat_old(const RVector<double>& x,         // observed counts
//                          const RVector<double>& means,     // Poisson rates
//                          RVector<double>& res,             // output buffer
//                          const int lg)                     // log=TRUE?
// {
//   std::size_t n = x.length();
//   
//   for (std::size_t i = 0; i < n; ++i) {
//     double count  = std::round(x[i]);     // match integer behavior
//     double lambda = means[i];             // rate parameter
//     
//     res[i] = -dpois2(count, lambda, lg);  // thread-safe Poisson backend
//   }
// }
// 

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
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);

    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 

    double res1=0;


    for(int i=0;i<m1;i++){
//      Rcpp::checkUserInterrupt();
      
//      if(progbar==1){ 
//        progress_bar2(i, m1-1);
//        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
//      };  
      
      
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


// arma::vec  f2_poisson_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
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
//   
//   Rcpp::NumericMatrix b2temp(l2,1);
//   
//   arma::mat x2(x.begin(), l1, l2, false); 
//   arma::mat alpha2(alpha.begin(), l1, 1, false); 
//   
//   Rcpp::NumericVector xb(l1);
//   arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
//   
//   // Note: Seem not to be used - Editing out
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
//   //arma::vec res2(res.begin(), m1, false);  // Shallow view into res
//   
//   arma::vec res2(m1);  // Owning allocation
// 
//   double res1=0;
//   
//   
//   for(int i=0;i<m1;i++){
//     //      Rcpp::checkUserInterrupt();
//     
//     //      if(progbar==1){ 
//     //        progress_bar2(i, m1-1);
//     //        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
//     //      };  
//     
//     
//     b2temp=b(Range(0,l2-1),Range(i,i));
// //    arma::mat b2(b2temp.begin(), l2, 1, false);
//     arma::mat b2(b2full.colptr(i), l2, 1, false);  // View of column i, size l2 × 1
//     
//     arma::mat P2(P.begin(), l2, l2, false); 
//     
//     bmu2=b2-mu2;
//     
//     res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
//     
//     xb2=exp(alpha2+ x2 * b2);
// 
//     // Wrap existing R memory buffers
//     RcppParallel::RVector<double> y_view(y);
//     RcppParallel::RVector<double> wt_view(wt);
//     RcppParallel::RVector<double> xb_view(xb);
//     RcppParallel::RVector<double> yy_view(yy); //must be preallocated to match y.size()
//     
//     // In-place evaluation using your log-scale accurate backend
//     neg_dpois_glmb_rmat_old(y_view, xb_view, yy_view,1.0);
//     
// 
//             
//     //yy=-dpois_glmb(y,xb,true);
//     
//     for(int j=0;j<l1;j++){
//       yy[j]=yy[j]*wt[j];  
//     }
//     
//   //  res(i) =std::accumulate(yy.begin(), yy.end(), res1);
//       res2(i) = std::accumulate(yy.begin(), yy.end(), res1);  
//   }
//   
// //  return res;      
//   return res2;      
//   
// }


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
        progress_bar2(i, m1-1);
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

