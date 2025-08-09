
//#include <Rcpp.h>
#include <vector>
#include <string>
#include "kernel_loader.h"
#include "OpenCL_helper.h"
#include "kernel_wrappers.h"
#include "kernel_runners.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

using namespace OpenCLHelper;


// // forward declare the runner
// void f2_binomial_logit_prep_kernel_runner(
//     const std::string& kernel_source,
//     const char*        kernel_name,
//     int                l1,
//     int                l2,
//     int                m1,
//     const std::vector<double>& X_flat,
//     const std::vector<double>& B_flat,
//     const std::vector<double>& mu_flat,
//     const std::vector<double>& P_flat,
//     const std::vector<double>& alpha_flat,
//     std::vector<double>&       qf_flat,
//     std::vector<double>&       xb_flat,
//     int progbar
// );

// [[Rcpp::export]]
Rcpp::List f2_binomial_logit_prep_opencl(
    Rcpp::NumericMatrix b,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    int progbar 
) {
  int l1 = x.nrow(), l2 = x.ncol(), m1 = b.ncol();
  
  // prepare flat inputs
  auto X_flat     = flattenMatrix(x);
  auto B_flat     = flattenMatrix(b);
  auto mu_flat    = flattenMatrix(mu);
  auto P_flat     = flattenMatrix(P);
  auto alpha_flat = copyVector(alpha);
  

  
  // prepare outputs
  std::vector<double> qf_flat(m1);
  std::vector<double> xb_flat((size_t)l1 * m1);
  
  Rcpp::NumericVector qf(m1);
  Rcpp::NumericMatrix xb(l1, m1);
  
   
#ifdef USE_OPENCL

  
  // load kernels
  //std::string core_src = load_kernel_library("f2_binomial_logit_prep"); 
  std::string ksrc     = load_kernel_source("src/f2_binomial_logit_prep.cl");
//  std::string all_src  = core_src + "\n" + ksrc;
  std::string all_src  = ksrc;
  
  // call runner
  f2_binomial_logit_prep_kernel_runner(
    all_src,
    "f2_binomial_logit_prep",
    l1, l2, m1,
    X_flat, B_flat, mu_flat, P_flat, alpha_flat,
    qf_flat, xb_flat,
    progbar
  );
  

  
  // reconstruct R outputs
  for (int j = 0; j < m1; ++j) qf[j] = qf_flat[j];
  
  for (int j = 0; j < m1; ++j)
    for (int i = 0; i < l1; ++i)
      xb(i, j) = xb_flat[(size_t)j * l1 + i];
#else
  Rcpp::Rcout << "[INFO] OpenCL not available — returning zero vector/matrices.\n";
  
#endif 
  
  return Rcpp::List::create(
    Rcpp::Named("xb") = xb,
    Rcpp::Named("qf") = qf
  );
}




// [[Rcpp::export]]
List f2_binomial_logit_prep_grad_opencl(
    const NumericMatrix& b,
    const NumericVector& y,
    const NumericMatrix& x,
    const NumericMatrix& mu,
    const NumericMatrix& P,
    const NumericVector& alpha,
    const NumericVector& wt,
    int                  progbar
) {
  // dimensions
  int l1 = x.nrow(), l2 = x.ncol(), m1 = b.ncol();
  
  // flatten inputs
  auto X_flat     = flattenMatrix(x);
  auto B_flat     = flattenMatrix(b);
  auto mu_flat    = flattenMatrix(mu);
  auto P_flat     = flattenMatrix(P);
  auto alpha_flat = copyVector(alpha);
  auto y_flat     = copyVector(y);
  auto wt_flat    = copyVector(wt);
  
  // allocate outputs
  std::vector<double> qf_flat(m1);
  std::vector<double> xb_flat((size_t)l1 * m1);
  std::vector<double> grad_flat((size_t)m1 * l2);


  NumericMatrix xb(l1, m1);
  NumericVector qf(m1);
  
  #ifdef USE_OPENCL

    
  // load & call kernel runner
  std::string ksrc    = load_kernel_source("src/f2_binomial_logit_prep_grad.cl");
  std::string all_src = ksrc;
 
  f2_binomial_logit_prep_grad_kernel_runner(
    all_src,
    "f2_binomial_logit_prep_grad",
    l1, l2, m1,
    X_flat, B_flat, mu_flat, P_flat, alpha_flat,
    y_flat, wt_flat,
    qf_flat, xb_flat, grad_flat,
    progbar
  );
  
  // rebuild xb, qf exactly as before
  for (int j = 0; j < m1; ++j) {
    qf[j] = qf_flat[j];
    for (int i = 0; i < l1; ++i) {
      xb(i, j) = xb_flat[(size_t)j * l1 + i];
    }
  }
 
#else
 Rcpp::Rcout << "[INFO] OpenCL not available — returning zero vector/matrices.\n";
 
#endif
  
  // wrap gradient directly as arma::mat (m1 rows × l2 cols),
  // no copies, column‐major data matches armadillo & R
  arma::mat grad_arma(
      grad_flat.data(),  // pointer to your flat array
      m1,                // n_rows
      l2,                // n_cols
      false,             // don't copy memory
      false              // strict = false
  );
  
  return List::create(
    Named("xb")   = xb,
    Named("qf")   = qf,
    Named("grad") = grad_arma
  );
}


