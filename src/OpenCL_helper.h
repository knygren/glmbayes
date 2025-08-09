#pragma once

//#include <Rcpp.h>
#include "RcppArmadillo.h"

#include <vector>

namespace OpenCLHelper {

// Flatten an Rcpp::NumericMatrix (column-major) into a std::vector<double>
std::vector<double> flattenMatrix(const Rcpp::NumericMatrix& mat);

// Copy an Rcpp::NumericVector into a std::vector<double>
std::vector<double> copyVector(const Rcpp::NumericVector& vec);

}