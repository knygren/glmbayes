//#include <Rcpp.h>
#include <vector>
#include <string>
#include "openclPort.h"
#include <RcppArmadillo.h>
#include "famfuncs.h"
#include "opencl.h"

using namespace Rcpp;

using namespace openclPort;
using namespace glmbayes::opencl;

namespace glmbayes {

namespace opencl {

Rcpp::List f2_f3_opencl(
    std::string family,
    std::string link,
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar
) {
  
  [[maybe_unused]] int l1 = x.nrow();
  [[maybe_unused]] int l2 = x.ncol();
  int m1 = b.ncol();
  
//  int l1 = x.nrow(), l2 = x.ncol(), m1 = b.ncol();
  
  
  
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
  std::vector<double> grad_flat(static_cast<size_t>(m1) * l2);
  
  
  Rcpp::NumericVector qf(m1);
  
  
  // Dispatch kernel name and source
  std::string kernel_name;
  // std::string kernel_file;
  std::string all_src;
  
  
#ifdef USE_OPENCL

  if (family == "binomial"||family == "quasibinomial") {
    if (link == "logit") {
      kernel_name = "f2_f3_binomial_logit";
      // kernel_file = "src/f2_f3_binomial_logit.cl";
    } 
    else if (link == "probit") {
      kernel_name = "f2_f3_binomial_probit";
      // kernel_file = "src/f2_f3_binomial_probit.cl";
    }
    else if (link == "cloglog") {
      kernel_name = "f2_f3_binomial_cloglog";
      // kernel_file = "src/f2_f3_binomial_cloglog.cl";
    }
    else {
      Rcpp::stop("Unsupported link function for binomial family: " + link);
    }
  }
  
  else if (family =="poisson"||family =="quasipoisson"){
    kernel_name = "f2_f3_poisson";
    // kernel_file  = "src/f2_f3_poisson.cl";
  }
  
  else if (family=="Gamma"){
    kernel_name = "f2_f3_gamma";
    // kernel_file  = "src/f2_f3_gamma.cl";
  }
  
  else if (family=="gaussian"){
    kernel_name = "f2_f3_gaussian";
    // kernel_file  = "src/f2_f3_gaussian.cl";
  }
  
  else {
    Rcpp::stop("Unsupported family: " + family);
  }

  // all_src = load_likelihood_subgradient_program(family, link, "glmbayes");
  all_src = load_likelihood_subgradient_program_v2(family, link, "glmbayes");

  // Legacy glmbayes program assembly (replaced by load_likelihood_subgradient_program above)
  // std::string OPENCL_source     = load_kernel_source("OPENCL.cl");
  // std::string rmath_source     = load_kernel_library("rmath","glmbayes", false);
  // std::string nmath_source     = load_kernel_library("nmath","glmbayes", false);
  // std::string dpq_source     = load_kernel_library("dpq","glmbayes", false);
  // std::string ksrc    = load_kernel_source(kernel_file);
  //
  // all_src = OPENCL_source +
  //   "\n" +   rmath_source +
  //   "\n" + dpq_source +
  //   "\n" +nmath_source
  // + "\n" +   ksrc;
  
  // Rcpp::Rcout << "Entering f2_f3_kernel runner \n";
  
  
  f2_f3_kernel_runner(
    all_src,
    kernel_name.c_str(),  // ✅ convert to const char*
    l1, l2, m1,
    X_flat, B_flat, mu_flat, P_flat, alpha_flat,
    y_flat, wt_flat,
    qf_flat, 
    // xb_flat,
    grad_flat,
    progbar
  );
  
  // Rcpp::Rcout << "Exiting f2_f3_kernel runner \n";
  
  
  // rebuild xb, qf exactly as before
  for (int j = 0; j < m1; ++j) {
    qf[j] = qf_flat[j];
    // for (int i = 0; i < l1; ++i) {
    //   xb(i, j) = xb_flat[static_cast<size_t>(j) * l1 + i];
    // }
  }
  
#else
  
  Rcpp::Rcout << "[INFO] OpenCL not available — returning zero vector/matrices.\n";
  
#endif
  
  // wrap gradient directly as arma::mat (m1 rows × l2 cols),
  // no copies, column-major data matches Armadillo & R
  arma::mat grad_arma(
      grad_flat.data(),  // pointer to your flat array
      m1,                // n_rows - dimensions
      l2,                // n_cols - grid points
      false,             // don't copy memory
      false              // strict = false
  );
  
  return Rcpp::List::create(
    // Rcpp::Named("xb")   = xb,
    Rcpp::Named("qf")   = qf,
    Rcpp::Named("grad") = grad_arma
  );
}

}
}

