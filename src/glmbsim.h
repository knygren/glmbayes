// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;

double ctrnorm_cpp(double lgrt,double lglt,double mu,double sigma);

void progress_bar2(double x, double N);

Rcpp::List  rnnorm_reg_std_cpp(int n,NumericVector y,NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,Function f2,Rcpp::List  Envelope,Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, int progbar=1);

Rcpp::List rnnorm_reg_cpp(int n,NumericVector y,NumericMatrix x,NumericVector mu,NumericMatrix P,NumericVector offset2,NumericVector wt,double dispersion,Rcpp::List famfunc, Function f1,Function f2,Function f3,NumericVector start,std::string family="binomial",std::string link="logit",int Gridtype=2);

Rcpp::List rnorm_reg_cpp(int n,NumericVector y,NumericMatrix x, NumericVector mu,NumericMatrix P,NumericVector offset2,NumericVector wt,double dispersion,Rcpp::List famfunc, Function f1,Function f2,Function f3,NumericVector start,std::string family="binomial",std::string link="logit",int Gridtype=2);


