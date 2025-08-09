#pragma once

#include <Rcpp.h>

// Rcpp‐exported wrapper that loads your .cl, flattens inputs, calls the runner,
// and reconstructs the xb/qf outputs for cpp.
//
// b      : l2 × m1 grid of β values
// y      : length l1 responses (unused in prep)
// x      : l1 × l2 design matrix
// mu     : l2 × 1 prior mean
// P      : l2 × l2 prior precision
// alpha  : length l1 offsets
// wt     : length l1 weights (unused here)
// progbar: 0 = no text bar, 1 = show progress
//
// Returns a List with components:
//   xb : NumericMatrix (l1 × m1) logistic preps
//   qf : NumericVector (length m1) quadratic forms
// 
// Implementation lives in f2_binomial_logit_prep_wrapper.cpp
// and calls f2_binomial_logit_prep_kernel_runner().
Rcpp::List f2_binomial_logit_prep_opencl(
    Rcpp::NumericMatrix b,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    int                  progbar = 0
);

// Rcpp‐exported wrapper that loads your .cl, flattens inputs, calls the
// f2_binomial_logit_prep_grad_kernel_runner, and reconstructs xb/qf/grad
// for downstream C++.
//
// b      : NumericMatrix, size l2 × m1 grid of β values
// y      : NumericVector, length l1 responses
// x      : NumericMatrix, size l1 × l2 design matrix
// mu     : NumericMatrix, size l2 × 1 prior mean
// P      : NumericMatrix, size l2 × l2 prior precision
// alpha  : NumericVector, length l1 offsets
// wt     : NumericVector, length l1 observation weights
// progbar: int, 0 = no text bar, 1 = show progress
//
// Returns an Rcpp::List with components:
//   xb   : NumericMatrix (l1 × m1) logistic preps
//   qf   : NumericVector  (length = m1) prior quadratic forms
//   grad : NumericMatrix  (m1 × l2) gradient w.r.t. β matching f3 output
//
// Implementation lives in f2_binomial_logit_prep_grad_wrapper.cpp
// and calls f2_binomial_logit_prep_grad_kernel_runner().
Rcpp::List f2_binomial_logit_prep_grad_opencl(
    Rcpp::NumericMatrix b,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    int                  progbar = 0
);