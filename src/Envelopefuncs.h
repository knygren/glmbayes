// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;

List glmbenvelope_c(NumericVector bStar,NumericMatrix A,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,std::string family="binomial",std::string link="logit",int Gridtype=2, int n=1,bool sortgrid=false);  
List Set_Grid(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint);
List   setlogP(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3);
Rcpp::List Set_Grid_C(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint,Rcpp::NumericMatrix Down,Rcpp::NumericMatrix Up,Rcpp::NumericMatrix lglt,Rcpp::NumericMatrix lgrt,Rcpp::NumericMatrix lgct,Rcpp::NumericMatrix logU,Rcpp::NumericMatrix logP);
Rcpp::List   setlogP_C(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst);


void Set_Grid_C2(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint,Rcpp::NumericMatrix Down,Rcpp::NumericMatrix Up,Rcpp::NumericMatrix lglt,Rcpp::NumericMatrix lgrt,Rcpp::NumericMatrix lgct,Rcpp::NumericMatrix logU,Rcpp::NumericMatrix logP);
void setlogP_C2(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst);
