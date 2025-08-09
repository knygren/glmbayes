#include "OpenCL_proxy_kernels.h"
#include <Rcpp.h>      // Rcpp::checkUserInterrupt, progress_bar2
#include <vector>
#include <cmath>

namespace OpenCLProxy {

void f2_binomial_logit_prep_kernel_proxy(
    const std::vector<double>& X_flat,
    const std::vector<double>& B_flat,
    const std::vector<double>& mu_flat,
    const std::vector<double>& P_flat,
    const std::vector<double>& alpha_flat,
    int l1,
    int l2,
    int m1,
    std::vector<double>& qf_flat,
    std::vector<double>& xb_flat,
    int progbar
) {
  
  

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
  
  
  
  
  
  // allocate/reset outputs
  qf_flat.assign(m1,            0.0);
  xb_flat.assign((size_t)l1*m1, 0.0);
  
  
  
  
  
  // temporaries
  std::vector<double> bmu(l2), tmp(l2);
  
  for (int j = 0; j < m1; ++j) {
    Rcpp::checkUserInterrupt();
 
    
    // pointers into flattened buffers
    const double* bj  = &B_flat[(size_t)j * l2];
    const double* muj = mu_flat.data();
    
    // 2a) Quadratic form: 0.5 * (b_j – mu)' P (b_j – mu)
    for (int k = 0; k < l2; ++k) {
      bmu[k] = bj[k] - muj[k];
    }
    for (int r = 0; r < l2; ++r) {
      double acc = 0.0;
      size_t prow = (size_t)r * l2;
      for (int c = 0; c < l2; ++c) {
        acc += P_flat[prow + c] * bmu[c];
      }
      tmp[r] = acc;
    }
    double qval = 0.0;
    for (int k = 0; k < l2; ++k) {
      qval += bmu[k] * tmp[k];
    }
    qf_flat[j] = 0.5 * qval;
    
    // 2b) Logistic prep: π_{i,j} = 1/(1 + exp(α[i] + X[i,·]·b_j))
    for (int i = 0; i < l1; ++i) {
      double linpred = alpha_flat[i];
      for (int k = 0; k < l2; ++k) {
        linpred += X_flat[(size_t)k * l1 + i] * bj[k];
      }
      xb_flat[(size_t)j * l1 + i] = 1.0 / (1.0 + std::exp(-linpred));
    }
  }
}

}