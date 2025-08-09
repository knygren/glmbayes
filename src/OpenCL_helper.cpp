#include "OpenCL_helper.h"

namespace OpenCLHelper {

std::vector<double> flattenMatrix(const Rcpp::NumericMatrix& mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  std::vector<double> out;
  out.reserve(static_cast<size_t>(nrow) * ncol);
  
  // Column-major traversal
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      out.push_back(mat(i, j));
    }
  }
  
  return out;
}

std::vector<double> copyVector(const Rcpp::NumericVector& vec) {
  return std::vector<double>(vec.begin(), vec.end());
}

}