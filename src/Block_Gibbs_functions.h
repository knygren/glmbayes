// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;

double Initialize_bstar(const NumericVector& y, arma::mat& x2,arma::mat& mu2, const arma::mat& P2,const arma::mat& alpha2,const NumericVector& wt,arma::mat& b2,NumericVector& xb,arma::colvec& xb2,arma::mat& Ptemp2,arma::mat& bmu2,arma::vec& bstar2,NumericVector& yy,std::string family="binomial",std::string link="logit");

double Find_Value(const NumericVector& y,arma::mat& x2,arma::mat& mu2,const arma::mat& P2,const arma::mat& alpha2, const NumericVector& wt, const arma::vec& b2, NumericVector& xb,NumericVector& yy,arma::vec& grad2,arma::mat& Pout2,arma::mat& Varout,arma::colvec& xb2,arma::mat& bmu2,arma::colvec& xbtemp2,std::string family="binomial",std::string link="logit");

double set_candidate(const arma::vec& b2, const double& stepsize,const arma::mat& Pout2,const arma::mat& Varout,const arma::mat& P2,const arma::mat& bmu2,const arma::mat& alpha2,const arma::mat& x2,const arma::mat& xb2,const arma::mat& mu2,arma::vec& btemp2,arma::vec& bmutemp2,arma::colvec& xbtemp2,const NumericVector& y,const NumericVector& wt,NumericVector& xbtemp,NumericVector& yy, const double& res2,std::string family="binomial",std::string link="logit");

void set_Pout(const arma::vec& b2,const NumericVector& y, const arma::mat& alpha2,const int& l1,const arma::mat& P2, const arma::mat& x2,const NumericVector& wt,const NumericVector& xbtemp,arma::colvec& xbtemp2,arma::mat& xrow2,arma::mat& Pout2,std::string family="binomial",std::string link="logit");

Rcpp::List optPostMode(NumericVector y,NumericMatrix x,NumericVector mu,NumericMatrix P, NumericVector alpha,NumericVector wt,NumericVector b,NumericVector bstar,std::string family="binomial",std::string link="logit");

Rcpp::List glmbsim_NGauss2_cpp(int n,NumericVector y,NumericMatrix x, 
NumericVector mu,NumericMatrix P,NumericVector offset2,NumericVector wt,double dispersion,Rcpp::List famfunc, Function f1,Function f2,Function f3,NumericVector start,std::string family="binomial",std::string link="logit",int Gridtype=2);

double get_epsilon1(double rstar,double epsilonstar,double U_out,double alpha_out,double nstar, double gammastar,double tstar,double mu_constant);

double golden_r(double upper_bound,double lower_bound,double epsilonstar, double U_out,double alpha_out,double nstar, double gammastar, double tstar, double mu_constant);

double find_nstar(double upper_bound,double lower_bound,double rstar2,double epsilon,double U_out,double alpha_out,double gammastar,double t_star,double mu_const,double epsilon_converge);

double get_n(double gammastar,double trace_const, double lambda_star,double epsilon1, double epsilon_converge,double gammastar_lower,double mu_const,double beta_const, int type=0 );

Rcpp::List golden_n(double trace_const, double lambda_star,double epsilon1, double epsilon_converge,double gamma_star_lower,double mu_const, double beta_const);

arma::mat Mat_pow(arma::mat A, double k);

Rcpp::List set_nstar(NumericMatrix x, NumericMatrix P, NumericMatrix P_0,arma::vec mu_0, arma::vec mu_star,arma::vec beta_star,arma::vec beta_star2,double epsilon_converge,NumericMatrix PD);

Rcpp::List set_beta_const(arma::mat x2,arma::vec offset2b,arma::vec mu_star2,int n,NumericVector y,NumericMatrix xtemp, NumericVector mutemp,NumericMatrix P,NumericVector offset2,NumericVector wt,double dispersion,Rcpp::List famfunc, Function f1,Function f2,Function f3,NumericMatrix betatemp,NumericMatrix x,NumericVector mu,NumericMatrix P_0,NumericVector offset3,NumericVector wt3,std::string family="binomial",std::string link="logit",int Gridtype=2);

void Set_PD(const arma::vec& b2,const NumericVector& y, const arma::mat& alpha2,const int& l1,const arma::mat& P2, const arma::mat& x2,const NumericVector& wt,const NumericVector& xbtemp,arma::colvec& xbtemp2,arma::mat& Pout2,std::string family="binomial",std::string link="logit");