// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;

NumericVector dbinom_glmb( NumericVector x, NumericVector N, NumericVector means, int lg);
NumericVector  f1_binomial_logit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar=0);
arma::mat  f3_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar=0);

NumericVector  f1_binomial_probit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar=0);
arma::mat  f3_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar=0);

NumericVector  f1_binomial_cloglog(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0);
arma::mat  f3_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar=0);

NumericVector dpois_glmb( NumericVector x, NumericVector means, int lg);
NumericVector  f1_poisson(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0);
arma::mat  f3_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar=0);

NumericVector dgamma_glmb( NumericVector x, NumericVector shape, NumericVector scale, int lg);
NumericVector  f1_gamma(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0);
arma::mat  f3_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar=0);

NumericVector dnorm_glmb( NumericVector x, NumericVector means, NumericVector sds,int lg);
NumericVector  f1_gaussian(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_gaussian(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
arma::mat  f3_gaussian(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);


//List Set_Grid(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint);
//List   setlogP(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3);

//Rcpp::List Set_Grid_C(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint,Rcpp::NumericMatrix Down,Rcpp::NumericMatrix Up,Rcpp::NumericMatrix lglt,Rcpp::NumericMatrix lgrt,Rcpp::NumericMatrix lgct,Rcpp::NumericMatrix logU,Rcpp::NumericMatrix logP);
//Rcpp::List   setlogP_C(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst);


//void Set_Grid_C2(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, 
//Rcpp::NumericMatrix Lint,
//Rcpp::NumericMatrix Down,
//Rcpp::NumericMatrix Up,
//Rcpp::NumericMatrix lglt,
//Rcpp::NumericMatrix lgrt,
//Rcpp::NumericMatrix lgct,
//Rcpp::NumericMatrix logU,
//Rcpp::NumericMatrix logP);

//void setlogP_C2(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst);
