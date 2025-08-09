// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;

List EnvelopeBuild_c(NumericVector bStar,NumericMatrix A,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,std::string family="binomial",std::string link="logit",int Gridtype=2, int n=1,bool sortgrid=false,
                       bool use_opencl = false,        // Enables OpenCL acceleration during envelope construction
                         bool verbose = false            // Enables diagnostic output
                       
                       );  
