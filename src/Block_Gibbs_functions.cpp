// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "famfuncs.h"
#include "glmbsim.h"
#include "Envelopefuncs.h"

using namespace Rcpp;


// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar



double Initialize_bstar(const NumericVector& y, arma::mat& x2,
                        arma::mat& mu2, 
                        const arma::mat& P2,
                        const arma::mat& alpha2,const NumericVector& wt,
                        arma::mat& b2,NumericVector& xb,
                        arma::colvec& xb2,
                        arma::mat& Ptemp2,
                        arma::mat& bmu2,
                        arma::vec& bstar2,
                        NumericVector& yy,
                        std::string family="binomial",
                        std::string link="logit"){
  
  int l1 = x2.n_rows;
  
  arma::mat xrow2=x2.row(0);
  
  int j;
  double res2=0;
  
  /////////////////// binomial - logit /////////////////////////////
  
  if(family=="binomial" && link=="logit")
  {
    xb2=exp(-alpha2- x2 * b2);
    for(j=0;j<l1;j++){
      xb(j)=1/(1+xb(j));  
      xrow2=x2.row(j);
      
      Ptemp2=Ptemp2+wt(j)*xb(j)*(1-xb(j))*trans(xrow2)*xrow2;
      
      //Ptemp2(j,j)=Ptemp2(j,j)+wt(j)*xb(j)*(1-xb(j))*trans(xrow2)*xrow2;
    }
    
    
    if(arma::is_finite(bstar2)){
      b2=inv_sympd(Ptemp2)*((Ptemp2-P2)*bstar2+P2*mu2); 
      xb2=exp(-alpha2- x2 * b2);
      for(int j=0;j<l1;j++)
      {   xb(j)=1/(1+xb(j));  } 
    }
    
    bmu2=b2-mu2;
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dbinom_glmb(y,wt,xb,true);    
    res2=std::accumulate(yy.begin(), yy.end(), res1);
  }
  
  /////////////////// binomial - probit /////////////////////////////
  
  if(family=="binomial" && link=="probit")
  {
    
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector d1(l1);
    
    xb2=alpha2+ x2 * b2;
    
    p1=pnorm(xb,0.0,1.0);
    p2=pnorm(-xb,0.0,1.0);
    d1=dnorm(xb,0.0,1.0);
    
    // This part must be edited (using Hessian)
    
    for(int j=0;j<l1;j++){
      xrow2=x2.row(j);
      Ptemp2=Ptemp2
        +wt(j)*d1(j)*(y(j)*(d1(j)+xb(j)*p1(j))/(p1(j)*p1(j))
                        +(1-y(j))*(d1(j)-xb(j)*p2(j))/(p2(j)*p2(j)))*trans(xrow2)*xrow2;
    }
    
    xb=pnorm(xb,0.0,1.0);
    
    
    ///////////////////////////
    
    
    if(arma::is_finite(bstar2)){
      b2=inv_sympd(Ptemp2)*((Ptemp2-P2)*bstar2+P2*mu2); 
      xb2=alpha2+ x2 * b2;
      xb=pnorm(xb,0.0,1.0);
    }
    
    bmu2=b2-mu2;
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dbinom_glmb(y,wt,xb,true);    
    res2=std::accumulate(yy.begin(), yy.end(), res1);
  }
  
  
  
  /////////////////// binomial - cloglog /////////////////////////////
  
  if(family=="binomial" && link=="cloglog")
  {
    
    NumericVector exb(l1);
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector atemp(l1);
    
    xb2=alpha2+ x2 * b2;
    exb=exp(xb);
    
    
    // This part must be edited (using Hessian)
    
    
    for(j=0;j<l1;j++){
      //      p1(j)=1-exp(-exb(j));
      //      p2(j)=exp(-exb(j));
      
      atemp(j)=exp(exb(j))-1;
      
      xrow2=x2.row(j);
      Ptemp2=Ptemp2+wt(j)*(1-y(j))*exb(j)*trans(xrow2)*xrow2
        +wt(j)*y(j)*(exb(j)*exb(j)*exp(exb(j))/(atemp(j)*atemp(j) ))*trans(xrow2)*xrow2
        -wt(j)*y(j)*(exb(j)/atemp(j))*trans(xrow2)*xrow2;
        
        
    }
    
    
    
    ///////////////////////////
    
    
    if(arma::is_finite(bstar2)){
      b2=inv_sympd(Ptemp2)*((Ptemp2-P2)*bstar2+P2*mu2); 
      xb2=alpha2+ x2 * b2;
      exb=exp(xb);
      for(j=0;j<1;j++){
        xb(j)=1-exp(-exb(j));
      }
      
    }
    
    bmu2=b2-mu2;
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dbinom_glmb(y,wt,xb,true);    
    res2=std::accumulate(yy.begin(), yy.end(), res1);
  }
  
  
  /////////////////// quasibinomial - logit /////////////////////////////
  
  if(family=="quasibinomial"  && link=="logit")
  {
    xb2=exp(-alpha2- x2 * b2);
    for(j=0;j<l1;j++){
      xb(j)=1/(1+xb(j));  
      xrow2=x2.row(j);
      Ptemp2=Ptemp2+wt(j)*xb(j)*(1-xb(j))*trans(xrow2)*xrow2;}
    
    if(arma::is_finite(bstar2)){
      b2=inv_sympd(Ptemp2)*((Ptemp2-P2)*bstar2+P2*mu2); 
      xb2=exp(-alpha2- x2 * b2);
      for(int j=0;j<l1;j++)
      {   xb(j)=1/(1+xb(j));  } 
    }
    
    bmu2=b2-mu2;
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dbinom_glmb(y,wt,xb,true);    
    res2=std::accumulate(yy.begin(), yy.end(), res1);
  }
  
  
  
  /////////////////// poisson /////////////////////////////
  
  if(family=="poisson")
  {
    xb2=exp(alpha2+ x2 * b2);
    for(j=0;j<l1;j++){
      xrow2=x2.row(j);
      Ptemp2=Ptemp2+wt(j)*xb(j)*trans(xrow2)*xrow2;
    }
    
    if(arma::is_finite(bstar2)){
      b2=inv_sympd(Ptemp2)*((Ptemp2-P2)*bstar2+P2*mu2); 
      xb2=exp(alpha2+ x2 * b2);
      
      bmu2=b2-mu2;
      double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
      yy=-dpois_glmb(y,xb,true);    
      
      for(int j=0;j<l1;j++){
        yy[j]=yy[j]*wt[j];  
      }
      
      res2=std::accumulate(yy.begin(), yy.end(), res1);
    }
  }
  
  
  /////////////////// quasipoisson /////////////////////////////
  
  if(family=="quasipoisson")
  {
    xb2=exp(alpha2+ x2 * b2);
    for(j=0;j<l1;j++){
      xrow2=x2.row(j);
      Ptemp2=Ptemp2+wt(j)*xb(j)*trans(xrow2)*xrow2;
    }
    
    if(arma::is_finite(bstar2)){
      b2=inv_sympd(Ptemp2)*((Ptemp2-P2)*bstar2+P2*mu2); 
      xb2=exp(alpha2+ x2 * b2);
      
      bmu2=b2-mu2;
      double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
      yy=-dpois_glmb(y,xb,true);    
      
      for(int j=0;j<l1;j++){
        yy[j]=yy[j]*wt[j];  
      }
      
      res2=std::accumulate(yy.begin(), yy.end(), res1);
      
    }
  }
  
  /////////////////// Gamma /////////////////////////////
  
  if(family=="Gamma")
  {
    xb2=exp(-alpha2- x2 * b2);
    for(j=0;j<l1;j++){
      xrow2=x2.row(j);
      Ptemp2=Ptemp2+wt(j)*y(j)*xb(j)*trans(xrow2)*xrow2;
    }
    
    if(arma::is_finite(bstar2)){
      b2=inv_sympd(Ptemp2)*((Ptemp2-P2)*bstar2+P2*mu2); 
      xb2=exp(alpha2+ x2 * b2);
      
      bmu2=b2-mu2;
      
      for(int j=0;j<l1;j++){
        
        xb[j]=xb[j]/wt[j];  
      }
      
      double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
      yy=-dgamma_glmb(y,wt,xb,true);
      
      
      res2=std::accumulate(yy.begin(), yy.end(), res1);
    }
  }
  
  
  
  return(res2);
  
  
}



double Find_Value(const NumericVector& y,arma::mat& x2,arma::mat& mu2,
                  const arma::mat& P2,const arma::mat& alpha2, const NumericVector& wt, 
                  const arma::vec& b2, 
                  NumericVector& xb,
                  NumericVector& yy,
                  arma::vec& grad2,
                  arma::mat& Pout2,
                  arma::mat& Varout,
                  arma::colvec& xb2,
                  arma::mat& bmu2,
                  arma::colvec& xbtemp2,
                  std::string family="binomial",
                  std::string link="logit"
){
  
  int l1 = x2.n_rows;
  double res2=0;
  //  int l2 = x2.n_cols;
  
  arma::mat xrow2=x2.row(0);
  
  /////////////////// binomial - logit /////////////////////////////
  
  if(family=="binomial" && link=="logit")
  {
    
    xb2=exp(-alpha2- x2 * b2);
    bmu2=b2-mu2;
    
    for(int j=0;j<l1;j++){
      xb(j)=1/(1+xb(j));  
      xrow2=x2.row(j);
      Pout2=Pout2+wt(j)*xb(j)*(1-xb(j))*trans(xrow2)*xrow2;
    }
    
    
    Varout=inv_sympd(Pout2);
    
    
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dbinom_glmb(y,wt,xb,true);    
    res2=std::accumulate(yy.begin(), yy.end(), res1);
    
    for(int j=0;j<l1;j++){
      xb(j)=(xb(j)-y(j))*wt(j);
    }
    
    grad2=(P2 * bmu2+x2.t() * xb2);
    
  }
  
  /////////////////// binomial - probit /////////////////////////////
  
  if(family=="binomial" && link=="probit")
  {
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector d1(l1);
    
    xb2=alpha2+ x2 * b2;
    
    bmu2=b2-mu2;
    
    p1=pnorm(xb,0.0,1.0);
    p2=pnorm(-xb,0.0,1.0);
    d1=dnorm(xb,0.0,1.0);
    
    
    // Edit (Hessian)
    for(int j=0;j<l1;j++){
      xrow2=x2.row(j);
      Pout2=Pout2
        +wt(j)*d1(j)*(y(j)*(d1(j)+xb(j)*p1(j))/(p1(j)*p1(j))
                        +(1-y(j))*(d1(j)-xb(j)*p2(j))/(p2(j)*p2(j)))*trans(xrow2)*xrow2;
    }
    
    
    // This part should be good 
    xb=pnorm(xb,0.0,1.0);
    
    
    Varout=inv_sympd(Pout2);
    
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dbinom_glmb(y,wt,xb,true);    
    res2=std::accumulate(yy.begin(), yy.end(), res1);
    
    
    for(int j=0;j<l1;j++){
      xb(j)=(y(j)*d1(j)/p1(j)-(1-y(j))*d1(j)/p2(j))*wt(j);    
    }
    
    
    grad2=(P2 * bmu2-x2.t() * xb2);
    
    
  }
  
  /////////////////// binomial - cloglog /////////////////////////////
  
  if(family=="binomial" && link=="cloglog")
  {
    
    NumericVector exb(l1);
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector atemp(l1);
    
    xb2=alpha2+ x2 * b2;
    exb=exp(xb);
    
    bmu2=b2-mu2;
    
    for(int j=0;j<l1;j++){
      p1(j)=1-exp(-exb(j));
      xb(j)=1-exp(-exb(j));
      p2(j)=exp(-exb(j));
      
      atemp(j)=exp(exb(j))-1;
      
      xrow2=x2.row(j);
      Pout2=Pout2+wt(j)*(1-y(j))*exb(j)*trans(xrow2)*xrow2
        +wt(j)*y(j)*(exb(j)*exb(j)*exp(exb(j))/(atemp(j)*atemp(j) ))*trans(xrow2)*xrow2
        -wt(j)*y(j)*(exb(j)/atemp(j))*trans(xrow2)*xrow2;
        
    }
    
    
    Varout=inv_sympd(Pout2);
    
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dbinom_glmb(y,wt,xb,true);    
    res2=std::accumulate(yy.begin(), yy.end(), res1);
    
    xb2=alpha2+ x2 * b2;
    
    for(int j=0;j<l1;j++){
      atemp(j)=exp(xb(j)-exb(j));
      
      xb(j)=((y(j)*atemp(j)/p1(j))-((1-y(j))*atemp(j)/p2(j)))*wt(j);    
    }
    
    
    grad2=(P2 * bmu2-x2.t() * xb2);
    
    
  }
  
  
  
  
  /////////////////// quasi-binomial - logit /////////////////////////////
  
  if(family=="quasibinomial" && link=="logit")
  {
    
    xb2=exp(-alpha2- x2 * b2);
    bmu2=b2-mu2;
    
    for(int j=0;j<l1;j++){
      xb(j)=1/(1+xb(j));  
      xrow2=x2.row(j);
      Pout2=Pout2+wt(j)*xb(j)*(1-xb(j))*trans(xrow2)*xrow2;
    }
    
    
    Varout=inv_sympd(Pout2);
    
    
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dbinom_glmb(y,wt,xb,true);    
    res2=std::accumulate(yy.begin(), yy.end(), res1);
    
    for(int j=0;j<l1;j++){
      xb(j)=(xb(j)-y(j))*wt(j);
    }
    
    grad2=(P2 * bmu2+x2.t() * xb2);
    
  }    
  
  /////////////////// poisson /////////////////////////////
  
  if(family=="poisson" )
  {
    
    xb2=exp(alpha2+ x2 * b2);
    bmu2=b2-mu2;
    
    for(int j=0;j<l1;j++){  
      xrow2=x2.row(j);
      Pout2=Pout2+wt(j)*xb(j)*trans(xrow2)*xrow2;
    }
    
    
    Varout=inv_sympd(Pout2);
    
    
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dpois_glmb(y,xb,true);    
    for(int j=0;j<l1;j++){
      yy[j]=yy[j]*wt[j];  }
    res2=std::accumulate(yy.begin(), yy.end(), res1);
    
    
    for(int j=0;j<l1;j++){
      xb(j)=(y(j)-xb(j))*wt(j);  
    }
    
    
    grad2=(P2 * bmu2-x2.t() * xb2);
    
  }
  
  /////////////////// quasipoisson /////////////////////////////
  
  if(family=="quasipoisson" )
  {
    
    xb2=exp(alpha2+ x2 * b2);
    bmu2=b2-mu2;
    
    for(int j=0;j<l1;j++){  
      xrow2=x2.row(j);
      Pout2=Pout2+wt(j)*xb(j)*trans(xrow2)*xrow2;
    }
    
    
    Varout=inv_sympd(Pout2);
    
    
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dpois_glmb(y,xb,true);    
    for(int j=0;j<l1;j++){
      yy[j]=yy[j]*wt[j];  }
    res2=std::accumulate(yy.begin(), yy.end(), res1);
    
    
    for(int j=0;j<l1;j++){
      xb(j)=(y(j)-xb(j))*wt(j);  
    }
    
    
    grad2=(P2 * bmu2-x2.t() * xb2);
    
  }
  
  if(family=="Gamma" )
  {
    
    xb2=exp(alpha2+ x2 * b2);
    bmu2=b2-mu2;
    
    for(int j=0;j<l1;j++){
      xrow2=x2.row(j);
      Pout2=Pout2+(wt(j)*y(j)/xb(j))*trans(xrow2)*xrow2;
    }
    
    Varout=inv_sympd(Pout2);
    
    
    for(int j=0;j<l1;j++){
      xb[j]=xb[j]/wt[j];  
    }
    
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dgamma_glmb(y,wt,xb,true);
    
    res2=std::accumulate(yy.begin(), yy.end(), res1);
    
    
    xb2=exp(alpha2+ x2 * b2);
    
    for(int j=0;j<l1;j++){
      xb[j]=(1-y[j]/xb[j])*wt[j];
    }
    
    grad2= P2 * bmu2+x2.t() * xb2;
    
  }
  
  
  
  
  Rcpp::List out=Rcpp::List::create(Rcpp::Named("grad2")=grad2,
                                    Rcpp::Named("Pout2")=Pout2);  
  
  
  return(res2);
}




double set_candidate(const arma::vec& b2, const double& stepsize,
                     const arma::mat& Pout2,
                     const arma::mat& Varout,
                     const arma::mat& P2,const arma::mat& bmu2,
                     const arma::mat& alpha2,
                     const arma::mat& x2,const arma::mat& xb2,const arma::mat& mu2,
                     arma::vec& btemp2,arma::vec& bmutemp2,arma::colvec& xbtemp2,
                     const NumericVector& y,const NumericVector& wt,
                     NumericVector& xbtemp,NumericVector& yy, const double& res2
                       ,
                       std::string family="binomial",
                       std::string link="logit"){
  
  int l1 = x2.n_rows;
  double res3=0;
  
  //////////////   Set candidate point and check function value 
  
  /////////////////// binomial - logit /////////////////////////////
  
  if(family=="binomial" && link=="logit")
  {
    btemp2=b2-stepsize*Varout*(P2 * bmu2+x2.t() * xb2);    
    bmutemp2=btemp2-mu2;
    
    xbtemp2=exp(-alpha2- x2 * btemp2);
    
    for(int j=0;j<l1;j++){
      xbtemp(j)=1/(1+xbtemp(j));      
    } 
    
    
    double res1=0.5*arma::as_scalar(bmutemp2.t() * P2 *  bmutemp2);
    yy=-dbinom_glmb(y,wt,xbtemp,true);    
    res3=std::accumulate(yy.begin(), yy.end(), res1);
  }
  
  
  /////////////////// binomial - probit /////////////////////////////
  
  if(family=="binomial" && link=="probit")
  {
    
    btemp2=b2-stepsize*Varout*(P2 * bmu2-x2.t() * xb2);    
    
    bmutemp2=btemp2-mu2;
    
    xbtemp2=alpha2+ x2 * btemp2;
    xbtemp=pnorm(xbtemp,0.0,1.0);
    
    double res1=0.5*arma::as_scalar(bmutemp2.t() * P2 *  bmutemp2);
    yy=-dbinom_glmb(y,wt,xbtemp,true);    
    res3=std::accumulate(yy.begin(), yy.end(), res1);
    
    
  }
  
  /////////////////// binomial - probit /////////////////////////////
  
  if(family=="binomial" && link=="cloglog")
  {
    NumericVector exb(l1);
    
    
    btemp2=b2-stepsize*Varout*(P2 * bmu2-x2.t() * xb2);    
    bmutemp2=btemp2-mu2;
    
    xbtemp2=alpha2+ x2 * btemp2;
    
    exb=exp(xbtemp);
    
    for(int j=0;j<l1;j++){
      xbtemp(j)=1-exp(-exb(j));
    }
    
    
    double res1=0.5*arma::as_scalar(bmutemp2.t() * P2 *  bmutemp2);
    yy=-dbinom_glmb(y,wt,xbtemp,true);    
    res3=std::accumulate(yy.begin(), yy.end(), res1);
    
    
  }
  
  /////////////////// quasi-binomial - logit /////////////////////////////
  
  if(family=="quasibinomial" && link=="logit")
  {
    btemp2=b2-stepsize*Varout*(P2 * bmu2+x2.t() * xb2);    
    bmutemp2=btemp2-mu2;
    
    xbtemp2=exp(-alpha2- x2 * btemp2);
    
    for(int j=0;j<l1;j++){
      xbtemp(j)=1/(1+xbtemp(j));      
    } 
    
    
    double res1=0.5*arma::as_scalar(bmutemp2.t() * P2 *  bmutemp2);
    yy=-dbinom_glmb(y,wt,xbtemp,true);    
    res3=std::accumulate(yy.begin(), yy.end(), res1);
  }
  
  /////////////////// poisson /////////////////////////////
  
  if(family=="poisson" )
  {
    btemp2=b2-stepsize*Varout*(P2 * bmu2-x2.t() * xb2);    
    bmutemp2=btemp2-mu2;
    
    xbtemp2=exp(alpha2+ x2 * btemp2);
    
    
    double res1=0.5*arma::as_scalar(bmutemp2.t() * P2 *  bmutemp2);
    yy=-dpois_glmb(y,xbtemp,true); 
    for(int j=0;j<l1;j++){
      yy[j]=yy[j]*wt[j];  }
    res3=std::accumulate(yy.begin(), yy.end(), res1);
  }
  
  /////////////////// quasipoisson /////////////////////////////
  
  if(family=="quasipoisson" )
  {
    btemp2=b2-stepsize*Varout*(P2 * bmu2-x2.t() * xb2);    
    bmutemp2=btemp2-mu2;
    
    xbtemp2=exp(alpha2+ x2 * btemp2);
    
    
    double res1=0.5*arma::as_scalar(bmutemp2.t() * P2 *  bmutemp2);
    yy=-dpois_glmb(y,xbtemp,true); 
    for(int j=0;j<l1;j++){
      yy[j]=yy[j]*wt[j];  }
    res3=std::accumulate(yy.begin(), yy.end(), res1);
  }
  
  /////////////////// Gamma /////////////////////////////
  
  if(family=="Gamma" )
  {
    btemp2=b2-stepsize*Varout*(P2 * bmu2+x2.t() * xb2);    
    bmutemp2=btemp2-mu2;
    
    xbtemp2=exp(alpha2+ x2 * btemp2);
    
    for(int j=0;j<l1;j++){
      xbtemp2[j]=xbtemp2[j]/wt[j];  
    }
    
    double res1=0.5*arma::as_scalar(bmutemp2.t() * P2 *  bmutemp2);
    yy=-dgamma_glmb(y,wt,xbtemp,true);
    
    
    res3=std::accumulate(yy.begin(), yy.end(), res1);
  }
  
  
  
  return(res3);
  
}




void set_Pout(const arma::vec& b2,const NumericVector& y, 
              const arma::mat& alpha2,
              const int& l1,const arma::mat& P2, const arma::mat& x2,
              const NumericVector& wt,const NumericVector& xbtemp,arma::colvec& xbtemp2,
              arma::mat& xrow2,
              arma::mat& Pout2,
              std::string family="binomial",
              std::string link="logit"
){
  
  Pout2=P2;
  
  /////////////////// binomial - logit /////////////////////////////
  
  if(family=="binomial" && link=="logit")
  {
    for(int j=0;j<l1;j++){
      xrow2=x2.row(j);
      Pout2=Pout2+wt(j)*xbtemp(j)*(1-xbtemp(j))*trans(xrow2)*xrow2;
    }
  }
  
  
  /////////////////// binomial - probit /////////////////////////////
  
  if(family=="binomial" && link=="probit")
  {
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector d1(l1);
    
    // Moved to before the evaluation of p1,p2, d1 (after likely problematic - Could be cause of weird earlier behavior)
    
    xbtemp2=alpha2+ x2 * b2;
    
    p1=pnorm(xbtemp,0.0,1.0);
    p2=pnorm(-xbtemp,0.0,1.0);
    d1=dnorm(xbtemp,0.0,1.0);
    
    //        xbtemp2=alpha2+ x2 * b2;
    
    // Edit (Hessian)
    for(int j=0;j<l1;j++){
      xrow2=x2.row(j);
      Pout2=Pout2
        +wt(j)*d1(j)*(y(j)*(d1(j)+xbtemp(j)*p1(j))/(p1(j)*p1(j))
                        +(1-y(j))*(d1(j)-xbtemp(j)*p2(j))/(p2(j)*p2(j)))*trans(xrow2)*xrow2;
                        
    }
  }
  
  
  
  /////////////////// quasi-binomial - logit /////////////////////////////
  
  if(family=="binomial" && link=="cloglog")
  {
    
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector atemp(l1);
    NumericVector exb(l1);
    
    xbtemp2=alpha2+ x2 * b2;
    
    exb=exp(xbtemp);
    
    for(int j=0;j<l1;j++){
      
      
      atemp(j)=exp(exb(j))-1;
      
      xrow2=x2.row(j);
      Pout2=Pout2+wt(j)*(1-y(j))*exb(j)*trans(xrow2)*xrow2
        +wt(j)*y(j)*(exb(j)*exb(j)*exp(exb(j))/(atemp(j)*atemp(j) ))*trans(xrow2)*xrow2
        -wt(j)*y(j)*(exb(j)/atemp(j))*trans(xrow2)*xrow2;
        
    }
    
    
  }  
  
  /////////////////// binomial - logit /////////////////////////////
  
  if(family=="quasibinomial" && link=="logit")
  {
    for(int j=0;j<l1;j++){
      xrow2=x2.row(j);
      Pout2=Pout2+wt(j)*xbtemp(j)*(1-xbtemp(j))*trans(xrow2)*xrow2;
    }
  }
  
  
  /////////////////// poisson /////////////////////////////
  
  if(family=="poisson")
  {
    for(int j=0;j<l1;j++){
      xrow2=x2.row(j);
      Pout2=Pout2+wt(j)*xbtemp(j)*trans(xrow2)*xrow2;
    }
  }
  
  /////////////////// quasipoisson /////////////////////////////
  
  if(family=="quasipoisson")
  {
    for(int j=0;j<l1;j++){
      xrow2=x2.row(j);
      Pout2=Pout2+wt(j)*xbtemp(j)*trans(xrow2)*xrow2;
    }
  }
  
  
  /////////////////// Gamma /////////////////////////////
  
  if(family=="Gamma")
  {
    for(int j=0;j<l1;j++){
      xrow2=x2.row(j);
      Pout2=Pout2+(wt(j)*y(j)/xbtemp(j))*trans(xrow2)*xrow2;
    }
  }
  
  
  
}








Rcpp::List optPostMode(NumericVector y,NumericMatrix x,
                       NumericVector mu,NumericMatrix P, NumericVector alpha,
                       NumericVector wt,NumericVector b,NumericVector bstar,std::string family="binomial",
                       std::string link="logit"
){
  int l1 = x.nrow(), l2 = x.ncol();
  
  NumericMatrix Pout(clone(P));
  NumericMatrix Varout1(clone(P));
  NumericMatrix bmu(l2,1);
  Rcpp::NumericVector xb(l1);
  NumericVector xrow = x( 0, _);
  NumericMatrix Ptemp(clone(P));
  NumericVector grad(l2);
  NumericVector gradb(l2);
  NumericVector btemp(l2);
  NumericVector bmutemp(l2);
  Rcpp::NumericVector xbtemp(l1);
  double maxgrad;
  double stepsize=1;
  double res3;
  int k;
  
  arma::mat y2(y.begin(), l1, 1, false);
  arma::mat x2(x.begin(), l1, l2, false); 
  arma::mat mu2(mu.begin(), l2, 1, false); 
  arma::mat P2(P.begin(), l2, l2, false); 
  arma::mat alpha2(alpha.begin(), l1, 1, false); 
  arma::mat b2(b.begin(), l2, 1, false);
  arma::vec bstar2(bstar.begin(), l2, false);
  
  arma::mat Pout2(Pout.begin(), l2, l2, false); 
  arma::mat Varout(Varout1.begin(), l2, l2, false); 
  arma::mat bmu2(bmu.begin(), l2, 1, false); 
  arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below  
  //    arma::vec xrow2(xrow.begin(),l2);
  arma::mat xrow2=x2.row(0);
  
  arma::mat Ptemp2(Ptemp.begin(), l2, l2, false);  
  arma::vec grad2(grad.begin(),l2);
  arma::vec gradb2(grad.begin(),l2);
  arma::vec btemp2(btemp.begin(), l2, 1, false);
  arma::vec bmutemp2(bmutemp.begin(), l2, 1, false);
  arma::colvec xbtemp2(xbtemp.begin(),l1,false); // Reuse memory - update both below  
  Rcpp::List valuet;
  
  double res_final;
  NumericVector yy(l1);    
  
  // Initialize bstar
  
  //    Rcpp::Rcout << "Enter Initialize_bstar:" << std::endl << res_final << std::endl;
  
  double res2=Initialize_bstar(y, x2, mu2,P2, alpha2,wt,b2, xb,xb2,
                               Ptemp2,bmu2,bstar2,yy,family,link);
  
  //    Rcpp::Rcout << "Exit Initialize_bstar:" << std::endl << res2 << std::endl;
  //    b2.print("b2 - Initial value");
  
  ///////////////////////////////////////////////////////
  
  
  // Initialize while loop
  
  int i=0;
  int reset=0;
  int check=0;
  int check2=0;
  
  while(i<5 && check==0){
    
    /////////////////////////////////////////////////////////////
    
    // Calculate Function Value and gradient At Latest Iteration 
    
    
    Pout2=P2;
    
    //    Rcpp::Rcout << "Enter Find Value:" << std::endl << res_final << std::endl;
    
    res2=Find_Value(y,x2, mu2, P2,alpha2,  wt,  b2, xb,yy,grad2,Pout2,Varout,
                    xb2,bmu2,xbtemp2,family,link);
    //    Rcpp::Rcout << "Exit Find Value:" << std::endl << res_final << std::endl;
    //      grad2.print("Value for gradient:");
    //        b2.print("b2 value after find value");      
    
    res_final=res2;
    
    reset=0;
    
    if(arma::any(grad2)==false){
      check=1;
    }
    // Why is this using P2 instead of Pout2?
    
    //    grad2.print("Value for gradient");
    
    //    gradb2=inv_sympd(P2)*grad2;
    
    gradb2=inv_sympd(Pout2)*grad2;
    
    //    gradb2.print("Value for gradb2 - Unstandardized");
    
    gradb2=abs(2*gradb2/(2*b2+gradb2));
    //    gradb2.print("Value for gradb2 - standardized");
    maxgrad=max(gradb2);
    if(maxgrad<0.0001){
      check=1;
    }
    
    
    
    
    // Update b2 (Newton-Rhapson update)
    // Reduce stepsize if function value increases
    
    stepsize=1;
    
    check2=0;
    k=0;
    
    while(check2==0&& k<5 && check==0){
      
      //////////////   Set candidate point and check function value 
      
      //    Rcpp::Rcout << "Enter set_candidate:" << std::endl << res2 << std::endl;
      //    b2.print("b2 entering set_candidate:");
      
      
      res3=set_candidate( b2,  stepsize, Pout2, Varout,P2, bmu2, alpha2, 
                          x2,xb2, mu2, btemp2, bmutemp2, xbtemp2, y, wt, xbtemp, yy,res2,family,link);
      //    btemp2.print("btemp2 exiting set_candidate:");
      
      //    b2.print("b2 exiting set candidate");
      //      btemp2.print("btemp2 exiting set candidate");
      //      Rcpp::Rcout << "Value for res2 exiting set candidate:" << std::endl << res2 << std::endl;
      //      Rcpp::Rcout << "Value for res3 exiting set candidate:" << std::endl << res3 << std::endl;
      //      Rcpp::Rcout << "Proposed change in value:" << std::endl << res3-res2 << std::endl;
      
      if(res3<res2){
        b2= btemp2;
        reset=1;
        res_final=res3;
        
        check2=1;
        
      }
      
      else{
        stepsize=stepsize/2.0;
      } 
      
      k++;
    }
    
    
    i++;
    
    
  }
  
  // If needed - recalculate Pout - Only if end of loop without full convergence
  
  if(reset==1){
    //            Rcpp::Rcout << "Enter set_Pout:" << std::endl << res2 << std::endl;
    
    set_Pout(b2,y,alpha2,l1,P2,x2,wt,xbtemp,xbtemp2,xrow2,Pout2,family,link);
    //          Rcpp::Rcout << "Exit set_Pout:" << std::endl << res2 << std::endl;
  }
  
  //    b2.print("b2 final value-New optimization:");
  
  Rcpp::List opt=Rcpp::List::create(Rcpp::Named("bstar")=b,
                                    Rcpp::Named("Pout")=Pout,Rcpp::Named("minval")=res_final);  
  
  return(opt);
}



// glmbsim_NGauss2_cpp may have modifications used by Two-Block Gibbs samplers


// // [[Rcpp::export]]

Rcpp::List glmbsim_NGauss2_cpp(int n,NumericVector y,NumericMatrix x, 
                               NumericVector mu,NumericMatrix P,NumericVector offset2,NumericVector wt,
                               double dispersion,Rcpp::List
                                 famfunc, Function f1,Function f2,Function f3,NumericVector start,
                                 std::string family="binomial",
                                 std::string link="logit",
                                 int Gridtype=2      
) {
  
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function asVec("as.vector");
  int l1=x.ncol();
  int l2=x.nrow();
  
  int l1b=mu.length();
  int l1c=P.ncol();
  int l1d=P.nrow();
  
  if(l1b!=l1) Rcpp::stop("Number of rows in mu not consistent with number of columns in matrix x");;
  if(l1c!=l1) Rcpp::stop("Number of columns in matrix P not consistent with number of columns in matrix x");;
  if(l1d!=l1) Rcpp::stop("Number of rows in matrix P not consistent with number of columns in matrix x");;
  
  int l2b=y.length();
  int l2c=offset2.length();
  int l2d=wt.length();
  
  if(l2b!=l2) Rcpp::stop("Number of rows in y not consistent with number of rows in matrix x");;
  if(l2c!=l2) Rcpp::stop("Number of rows in offset2 vector not consistent with number of rows in matrix x");;
  if(l2d!=l2) Rcpp::stop("Number of rows in wt vector not consistent with number of rows in matrix x");;
  
  
  double dispersion2;
  NumericVector alpha(l2);
  NumericMatrix mu2a=asMat(mu);
  
  arma::mat x2(x.begin(), l2, l1, false);
  arma::vec alpha2(alpha.begin(),l2,false);  
  arma::vec offset2b(offset2.begin(),l2,false);  
  arma::mat mu2(mu2a.begin(), mu2a.nrow(), mu2a.ncol(), false);
  
  arma::mat P2(P.begin(), P.nrow(), P.ncol(), false);
  
  if(family=="poisson"||family=="binomial") dispersion2=1;
  else dispersion2=dispersion;
  
  NumericVector  wt2=wt/dispersion2;
  
  alpha2=x2*mu2+offset2b;
  
  int i;
  
  
  NumericVector parin=start-mu;
  NumericVector mu1=mu-mu;
  Rcpp::Function optfun("optim");
  NumericVector bstar(l2);
  
  NumericMatrix b2a(l1);
  NumericVector parin2(clone(parin));
  
  
  if(family=="binomial" && link=="logit"){bstar=log(y/(1-y))-alpha;}  
  if(family=="quasibinomial" && link=="logit"){bstar=log(y/(1-y))-alpha;}  
  if(family=="binomial" && link=="probit"){bstar=qnorm(y,0.0,1.0)-alpha;}  
  if(family=="quasibinomial" && link=="probit"){bstar=qnorm(y,0.0,1.0)-alpha;}  
  if(family=="binomial" && link=="cloglog"){bstar=log(-log(1-y))-alpha;}  
  if(family=="quasibinomial" && link=="cloglog"){bstar=log(-log(1-y))-alpha;}  
  if(family=="poisson"){bstar=log(y)-alpha;}  
  if(family=="quasipoisson"){bstar=log(y)-alpha;}  
  if(family=="Gamma"){bstar=log(y)-alpha;}  
  
  
  
  List opt1=optPostMode(y,x,mu1, P, alpha,wt2,
                        parin2,bstar,family,link);
  
  b2a=asMat(opt1(0));
  NumericMatrix A1=opt1(1);
  NumericVector min1=asVec(opt1[2]);
  int conver1=0;
  
  arma::mat b2(b2a.begin(), b2a.nrow(), b2a.ncol(), false);
  arma::mat A1_b(A1.begin(), l1, l1, false); 
  
  if(conver1>0){Rcpp::stop("Posterior Optimization failed");}
  
  arma::vec mu_0(mu.begin(), l1, false);
  
  arma::vec eigval_1;
  arma::mat eigvec_1;
  
  eig_sym(eigval_1, eigvec_1, A1_b);
  
  NumericMatrix L2Inv_1(l1, l1);
  arma::mat L2Inv(L2Inv_1.begin(), L2Inv_1.nrow(), L2Inv_1.ncol(), false);
  
  
  // Standardize Model to Have Diagonal Variance-Covariance Matrix at Posterior Mode
  
  arma::mat D1=arma::diagmat(eigval_1);
  arma::mat L2= arma::sqrt(D1)*trans(eigvec_1);
  L2Inv=eigvec_1*sqrt(inv_sympd(D1));
  arma::mat b3=L2*b2; 
  arma::mat mu3=L2*mu2;
  arma::mat x3=x2*L2Inv;
  arma::mat P3=trans(L2Inv)*P2*L2Inv;
  
  //   Find diagonal matrix that has "smaller" precision than prior  
  //   Follows Definition 3, and procedure on p. 1150 in Nygren 
  //   Puts model into standard form 
  
  
  arma::mat P3Diag=arma::diagmat(arma::diagvec(P3));
  arma::mat epsilon=P3Diag;
  arma::mat P4=P3Diag;   
  
  double scale=1;
  int check=0;
  arma::vec eigval_2;
  arma::mat eigvec_2;
  double eigval_temp;
  
  while(check==0){
    epsilon=scale*P3Diag;
    P4=P3-epsilon;				
    eig_sym(eigval_2, eigvec_2, P4);
    eigval_temp=arma::min(eigval_2);
    if(eigval_temp>0){check=1;}
    else{scale=scale/2;}
  }
  
  
  int check2=0;
  double scale2=scale;
  arma::mat epsilon_temp=P3Diag;   
  arma::mat P4_temp=P3Diag;   
  
  while(check2==0){
    
    scale=scale+(scale2/10);
    epsilon_temp=scale*P3Diag;
    P4_temp=P3-epsilon_temp;
    eig_sym(eigval_2, eigvec_2, P4_temp);
    eigval_temp=arma::min(eigval_2);
    
    if(eigval_temp>0){
      epsilon=epsilon_temp;
      P4=P4_temp;	
    }		
    else{check2=1;}
    
  }    
  
  arma::mat ident=arma::mat (l1,l1,arma::fill::eye);
  arma::mat A3=ident-epsilon;	
  
  //   Put into Standard form
  
  eig_sym(eigval_2, eigvec_2, epsilon);
  
  arma::mat D2=arma::diagmat(eigval_2);
  
  
  NumericMatrix b4_1(l1,1);
  NumericMatrix mu4_1(l1,1);
  NumericMatrix x4_1(x.nrow(), x.ncol());
  NumericMatrix A4_1(l1, l1);
  NumericMatrix P5_1(l1, l1);
  NumericMatrix P6Temp_1(l1, l1);
  NumericMatrix L3Inv_1(l1, l1);
  arma::mat b4(b4_1.begin(), b4_1.nrow(), b4_1.ncol(), false);
  arma::mat mu4(mu4_1.begin(), mu4_1.nrow(), mu4_1.ncol(), false);
  arma::mat x4(x4_1.begin(), x4_1.nrow(), x4_1.ncol(), false);
  arma::mat A4(A4_1.begin(), A4_1.nrow(), A4_1.ncol(), false);
  arma::mat P5(P5_1.begin(), P5_1.nrow(), P5_1.ncol(), false);
  arma::mat P6Temp(P6Temp_1.begin(), P6Temp_1.nrow(), P6Temp_1.ncol(), false);
  arma::mat L3Inv(L3Inv_1.begin(), L3Inv_1.nrow(), L3Inv_1.ncol(), false);
  
  
  arma::mat L3= arma::sqrt(D2)*trans(eigvec_2);
  L3Inv=eigvec_2*sqrt(inv_sympd(D2));
  b4=L3*b3; 
  mu4=L3*mu3; 
  x4=x3*L3Inv;
  A4=trans(L3Inv)*A3*L3Inv;
  P5=trans(L3Inv)*P4*L3Inv;
  P6Temp=P5+ident;  
  NumericVector b5=asVec(b4_1);
  Rcpp::List Envelope;
  
  NumericMatrix mu5_1=0*mu4_1; // Does this modify mu4_1?
  
  
  if(n==1){
    Envelope=glmbenvelope_c(b5, A4_1,y, x4_1,mu5_1,
                            P5_1,alpha,wt2,family,link,Gridtype, n,false);
  }
  if(n>1){
    Envelope=glmbenvelope_c(b5, A4_1,y, x4_1,mu5_1,P5_1,alpha,wt2,family,link,Gridtype, n,true);
  }
  
  Rcpp::List sim=glmbsim_cpp(n,y,x4_1,mu5_1,P5_1,alpha,wt2,f2,Envelope,family,link);
  
  
  NumericMatrix sim2=sim[0];
  arma::mat sim2b(sim2.begin(), sim2.nrow(), sim2.ncol(), false);
  NumericMatrix out(l1,n);
  arma::mat out2(out.begin(), out.nrow(), out.ncol(), false);
  
  out2=L2Inv*L3Inv*trans(sim2b);
  NumericVector LL(n);
  
  
  for(i=0;i<n;i++){
    out(_,i)=out(_,i)+mu;
    
    // How much does log-likelihood evaluation slows this down?
    LL[i]=as<double>(f1(_["b"]=out(_,i),_["y"]=y,_["x"]=x,offset2,wt2));
  }
  
  
  
  Rcpp::List outlist=Rcpp::List::create(
    Rcpp::Named("coefficients")=trans(out2),
    //  Rcpp::Named("Envelope")=Envelope,
    Rcpp::Named("loglike")=LL
  );  
  
  return(outlist);
}




double get_epsilon1(double rstar,double epsilonstar,double U_out,
                    double alpha_out,double nstar, double gammastar,double tstar,double mu_constant){
  
  double d1=pow(1-epsilonstar,rstar);
  
  double d2=pow(U_out,rstar)/pow(alpha_out,1-rstar);
  
  
  //  ropt<-function(rstar,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant){
  //    d1<-(1-epsilonstar)^rstar
  //    d2<-(U_out^rstar)/(alpha_out^(1-rstar))
  //    d1^nstar+(d2^nstar)*(1+gammastar*(1+tstar)+mu_constant*(1+tstar))
  
  double epsilonout=pow(d1,nstar)+pow(d2,nstar)*(1+gammastar*(1+tstar)+mu_constant*(1+tstar));
  
  return(epsilonout);
  
}



double golden_r(double upper_bound,double lower_bound,double epsilonstar, double U_out,double alpha_out,
                double nstar, double gammastar, double tstar, double mu_constant
){
  
  // Inititialize Variables 
  
  double golden_ratio = 2/(sqrt(5) + 1);
  //  double  upper_bound=rstar;
  
  //  double lower_bound=0;
  double tolerance=0.00001;
  
  // Use the golden ratio to set the initial test points
  
  double r1 = upper_bound - golden_ratio*(upper_bound-lower_bound);
  double r2 = lower_bound + golden_ratio*(upper_bound - lower_bound);
  
  
  double  f1=get_epsilon1(r1,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant);
  double  f2=get_epsilon1(r2,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant);
  
  
  
  int iteration = 0;
  
  double gap=upper_bound - lower_bound;
  
  //      f2<-ropt(r2,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant)
  
  
  while (gap > tolerance && iteration<30){
    
    
    iteration = iteration + 1;
    
    
    if (f2 > f1){
      // then the minimum is to the left of r2
      // let r2 be the new upper bound
      // let r1 be the new upper test point
      
      //### Set the new upper bound
      //### Set the new upper bound
      upper_bound= r2;  
      
      
      // Set the new upper test point
      // Use the special result of the golden ratio
      
      r2 = r1;
      
      f2 = f1;
      
      //### Set the new lower test point
      
      r1 = upper_bound - golden_ratio*(upper_bound - lower_bound);
      f1=  get_epsilon1(r1,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant);
      
      gap=  upper_bound - lower_bound;    
    }
    
    else{
      
      //# the minimum is to the right of x1
      //# let r1 be the new lower bound
      //# let r2 be the new lower test point
      
      //### Set the new lower bound
      lower_bound = r1;
      
      //### Set the new lower test point
      r1 = r2;
      
      f1 = f2;
      
      //### Set the new upper test point
      r2 = lower_bound + golden_ratio*(upper_bound - lower_bound);
      
      f2=get_epsilon1(r2,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant);
      gap=  upper_bound - lower_bound;    
      
      
    }
    
  }
  
  //### Use the mid-point of the final interval as the estimate of the optimzer
  
  double rstar_out = (lower_bound + upper_bound)/2;
  
  //  double rstar=1;
  
  return(rstar_out);
  
  
}




double find_nstar(double upper_bound,double lower_bound,double rstar2,
                  double epsilon,double U_out,double alpha_out,
                  double gammastar,
                  double t_star,
                  double mu_const,
                  double epsilon_converge){
  
  
  double golden_ratio = 2/(sqrt(5) + 1);
  
  // Use the golden ratio to set the initial test points
  
  double  n1 = upper_bound - golden_ratio*(upper_bound - lower_bound);
  double  n2 = lower_bound + golden_ratio*(upper_bound - lower_bound);
  
  
  // ### Evaluate the function at the test points
  
  double epsilon_temp2=get_epsilon1( rstar2,epsilon, U_out,alpha_out,n1,  gammastar,t_star, mu_const);
  double  f1 = (epsilon_temp2-epsilon_converge)*(epsilon_temp2-epsilon_converge);
  
  
  
  epsilon_temp2=get_epsilon1( rstar2,epsilon, U_out,alpha_out,n2,  gammastar,t_star, mu_const);
  double  f2 = (epsilon_temp2-epsilon_converge)*(epsilon_temp2-epsilon_converge);
  
  
  
  double gap=upper_bound-lower_bound;
  
  int iter1=0;
  
  
  while(gap>1 && iter1<20){
    
    iter1=iter1+1;
    
    
    if(f2>f1){
      
      // then the minimum is to the left of n2
      // let n2 be the new upper bound
      // let n1 be the new upper test point
      
      
      upper_bound = n2;
      
      n2=n1;
      
      f2=f1;
      
      n1 = upper_bound - golden_ratio*(upper_bound - lower_bound);
      
      epsilon_temp2=get_epsilon1( rstar2,epsilon, U_out,alpha_out,n1,  gammastar,t_star, mu_const);
      f1 = (epsilon_temp2-epsilon_converge)*(epsilon_temp2-epsilon_converge);
      
      gap=upper_bound-lower_bound;
    }
    
    else{
      
      
      //      Rcpp::Rcout << "f1>f2:  "  << std::endl;
      
      //      Rcpp::Rcout << "n1:                                " << std::flush << n1 << std::endl;
      //      Rcpp::Rcout << "n2:                                " << std::flush << n2 << std::endl;
      //      Rcpp::Rcout << "f1:                                " << std::flush << f1 << std::endl;
      //      Rcpp::Rcout << "f2:                                " << std::flush << f2 << std::endl;
      
      lower_bound = n1;
      
      n1=n2;
      
      f1=f2;
      
      n2 = lower_bound + golden_ratio*(upper_bound - lower_bound);
      
      epsilon_temp2=get_epsilon1( rstar2,epsilon, U_out,alpha_out,n2,  gammastar,t_star, mu_const);
      f2 = (epsilon_temp2-epsilon_converge)*(epsilon_temp2-epsilon_converge);
      
      gap=upper_bound-lower_bound;
      
    }
    
    //      Rcpp::Rcout << "lower_bound:  " << std::flush << lower_bound << std::endl;
    //    Rcpp::Rcout << "upper_bound:  " << std::flush << upper_bound << std::endl;
    
    
  }
  
  //  Rcpp::Rcout << "rstar_first:                                " << std::flush << rstar1 << std::endl;
  //  Rcpp::Rcout << "nstar_first:                                " << std::flush << nstar << std::endl;
  
  
  //  rstar1=rstar2;
  
  double nstar=(n1+n2)/2;
  
  return(nstar);
}




double get_n(double gammastar,double trace_const, double lambda_star,
             double epsilon1, double epsilon_converge,double gammastar_lower,double mu_const,
             double beta_const, int type=0 ){
  
  
  double t_star=sqrt(beta_const/gammastar);
  double gammastar2=(1+t_star)*gammastar;
  double trace_const2=(1+t_star)*trace_const;
  double mu_const2=(1+t_star)*mu_const;
  
  double beta_const2=0;
  if(t_star>0) beta_const2=beta_const/t_star;
  
  double gammastar_lower2=(1+t_star)*gammastar_lower;
  
  double epsilon=exp(log(epsilon1)-beta_const2-gammastar2);
  
  
  
  double alpha_out=(1+gammastar2)/(1+trace_const2+lambda_star*gammastar2);
  double U_out=(1+trace_const2+2*lambda_star*gammastar2);
  //    double A1_out=(1-exp(log(epsilon1)-beta_const2-gammastar2));
  double lg_A1_out=log(1-exp(log(epsilon1)-beta_const2-gammastar2));
  
  double x1=exp(log(epsilon1)-beta_const2-gammastar2);
  
  if(x1==0){
    
    Rcpp::Rcout << "gammastar"  << std::endl << gammastar  << std::endl ;
    Rcpp::Rcout << "trace_const"  << std::endl << trace_const  << std::endl ;
    Rcpp::Rcout << "lambda_star"  << std::endl << lambda_star  << std::endl ;
    
    
    Rcpp::stop("Unable to Calculate Finite Iteration Limit - Likely because Gammastar (trace_cont/(1-lambda_star)) is too large."    )  ;
    
    
  }
  
  double qc_lg_A1_out=-x1-0.5*x1*x1;
  
  
  
  if(lg_A1_out>-2.47036e-012)
  {
    lg_A1_out=-2.47036e-012;
    if(qc_lg_A1_out>-2.47036e-012)  lg_A1_out=qc_lg_A1_out;
    
  }
  
  
  
  
  
  //    double qc_lg_alpha_out=log(1+gammastar2)-log(1+trace_const2+lambda_star*gammastar2);
  double lg_alpha_out=log(1+(1/gammastar2))-log(lambda_star+( (1+trace_const2)/gammastar2));
  double x2=((lambda_star-1)+ (trace_const2)/gammastar2);
  
  double qc_lg_alpha_out2= -x2+0.5*x2*x2;
  
  
  
  
  if (lg_alpha_out<5.00028e-014){
    
    lg_alpha_out=5.00028e-014;
    if(qc_lg_alpha_out2<5.00028e-014) lg_alpha_out=qc_lg_alpha_out2;
    
  } 
  
  //    double rstar1=log(alpha_out)/(log(U_out)+log(alpha_out)-log(A1_out));
  double rstar1=lg_alpha_out/(log(U_out)+lg_alpha_out-lg_A1_out);
  //    double A3=pow(A1_out,rstar1);
  //    double log_A3_temp=rstar1*log(A1_out);
  double log_A3=rstar1*lg_A1_out;
  //    double log_A3=log(A3);
  double log_A3_2=-rstar1*exp(log(epsilon1)-beta_const2-gammastar2);
  
  
  
  
  if(log_A3==0){
    log_A3=log_A3_2;
  }
  
  
  
  // Initialize nstar     
  
  double nstar=((log(epsilon_converge))-log(2+gammastar_lower2+mu_const2))/log_A3;
  
  double rstar2=rstar1;  
  
  //    Rcpp::Rcout << "rstar_in:                                " << std::flush << rstar1 << std::endl;
  //    Rcpp::Rcout << "nstar_in:                                " << std::flush << nstar << std::endl;
  
  
  int i=0;
  
  while(i<10){    
    
    //    Rcpp::Rcout << "get_n_Iter:                                " << std::flush << i << std::endl;
    
    rstar2=golden_r(rstar1,0,epsilon, U_out,alpha_out,nstar, gammastar,t_star, mu_const);
    
    //    Rcpp::Rcout << "rstar2:                                " << std::flush << rstar2 << std::endl;
    
    
    double nstar_temp=nstar;
    
    //    get_epsilon1
    double epsilon_temp=get_epsilon1( rstar2,epsilon, U_out,alpha_out,nstar,  gammastar,t_star, mu_const);
    
    double epsilon_temp2=epsilon_temp;
    
    while(epsilon_temp2<epsilon_converge){
      
      nstar_temp=nstar_temp/2;
      
      epsilon_temp2=get_epsilon1( rstar2,epsilon, U_out,alpha_out,nstar_temp,  gammastar,t_star, mu_const);
      
    }
    
    
    double upper_bound=nstar;
    double lower_bound=nstar_temp;
    
    
    
    nstar=find_nstar(upper_bound,lower_bound,rstar2,
                     epsilon,U_out, alpha_out,gammastar,t_star,mu_const,epsilon_converge);
    
    //    Rcpp::Rcout << "nstar:                                " << std::flush << nstar << std::endl;
    
    
    rstar1=rstar2;
    
    
    
    //    Rcpp::Rcout << "rstar1:                                " << std::flush << rstar1 << std::endl;
    
    i=i+1;
    
    
    
  }
  
  
  
  //    double rstar3=golden_r(rstar1,0,epsilon, U_out,alpha_out,nstar, gammastar,t_star, mu_const);
  
  
  
  //    Rcpp::Rcout << "rstar_third:                                " << std::flush << rstar3 << std::endl;
  
  
  //    Rcpp::Rcout << "nstar_Out:                                " << std::flush << nstar << std::endl;
  //    Rcpp::Rcout << "rstar_Out:                                " << std::flush << rstar1 << std::endl;
  
  double out=nstar;
  
  if(type==1){out=rstar1;}
  
  
  return out;
  
}





Rcpp::List golden_n(double trace_const, double lambda_star,
                    double epsilon1, double epsilon_converge,
                    double gamma_star_lower,double mu_const, double beta_const){
  
  // Initialize Golden Section Search
  
  double gammastar=gamma_star_lower+0.1;
  double min=0;
  //    double low=0;
  //    double high=0;
  double val=0;
  double gamma_opt=gammastar;
  int upper_set=0;
  
  // Find Initial Value  
  
  
  val=get_n(gammastar,trace_const, lambda_star, epsilon1,  epsilon_converge,gamma_star_lower,mu_const, beta_const);
  
  double lower_bound=gamma_star_lower;
  min=val;
  
  //    Rcpp::Rcout << "############  Testing 2.1.2" << std::endl;
  
  
  
  // Rcpp::Rcout << "Gamma Opt: Initial Gammastar"  << std::endl << gammastar  << std::endl ;
  // Rcpp::Rcout << "Gamma Opt: Initial n"  << std::endl << val  << std::endl ;
  
  
  min=val;
  
  if(min==R_PosInf){
    //        Rcpp::stop("Min Value After Upper Is Set is Positive Infinity");
    Rcpp::warning("Initial Min Value is Positive Infinity");
    
  }
  
  //    double gamma_star_upper=0;
  
  //    low=val;
  //    high=val;
  
  // Find Upper Bound
  
  
  while(upper_set==0){
    
    Rcpp::checkUserInterrupt();
    
    //    gammastar=2*gammastar;
    
    gammastar=gammastar+1;
    
    
    //    gamma_star_upper=gammastar;
    
    //    Rcpp::Rcout << "check 1.1   "  << std::endl ;
    
    //    Rcpp::Rcout << "gammastar:Proposed"  << std::endl << gammastar  << std::endl ;
    
    // Temporarily edit this out        
    val=get_n(gammastar,trace_const, lambda_star, epsilon1,  epsilon_converge,gamma_star_lower,mu_const,beta_const);
    
    //    Rcpp::Rcout << "min-value"  << std::endl << min  << std::endl ;
    //      Rcpp::Rcout << "proposed val at upper"  << std::endl << val  << std::endl ;
    
    
    //    Rcpp::Rcout << "check 1.2   "  << std::endl ;
    
    if (val==min) {
      //      high=val;
      
      // Edit - Set min=val regardless
      min=val;
      gamma_opt=gammastar;
      
      upper_set=1;  
      
    }
    
    if (val>min) {
      //      high=val;
      
      // Edit - Set min=val regardless
      //    min=val;
      gamma_opt=gammastar;
      
      upper_set=1;  
      
    }
    
    if (val<min) {
      //      high=val;
      min=val;
      lower_bound=gamma_opt;
      gamma_opt=gammastar;
      // Was this causing an issue ?
      //    lower_set=1;  
      
    }
    
  }
  double upper_bound=gammastar;
  
  
  //    Rcpp::Rcout << "lower_bound  "  <<  std::endl << lower_bound <<std::endl;
  //    Rcpp::Rcout << "upper_bound  "  <<  std::endl << upper_bound <<std::endl;
  //    Rcpp::Rcout << "Value at current min  "  <<  std::endl << min <<std::endl;
  //    Rcpp::Rcout << "Value at upper_bound  "  <<  std::endl << val <<std::endl;
  
  
  if(min==R_PosInf){
    Rcpp::stop("Min Value After Upper Is Set is Positive Infinity");
    //      Rcpp::warning("Min Value After Upper Is Set is Positive Infinity");
  }
  
  
  
  
  double golden_ratio = 2/(sqrt(5) + 1);
  
  // Use the golden ratio to set the initial test points
  
  double  g1 = upper_bound - golden_ratio*(upper_bound - lower_bound);
  double  g2 = lower_bound + golden_ratio*(upper_bound - lower_bound);
  
  
  // ### Evaluate the function at the test points
  
  
  double  f1 = get_n(g1,trace_const, lambda_star, epsilon1,  epsilon_converge,gamma_star_lower,mu_const,beta_const);
  double  f2 = get_n(g2,trace_const, lambda_star, epsilon1,  epsilon_converge,gamma_star_lower,mu_const,beta_const);
  
  double gap=upper_bound-lower_bound;
  
  int iter1=0;
  
  while(gap>0.001 && iter1<20){
    
    iter1=iter1+1;
    
    //      Rcpp::Rcout << "iter  "  <<  std::endl << iter1 <<std::endl;
    
    if(f2>f1){
      
      //        Rcpp::Rcout << "f2>f1:  "  << std::endl;
      
      //        Rcpp::Rcout << "g1:                                " << std::flush << g1 << std::endl;
      //        Rcpp::Rcout << "g2:                                " << std::flush << g2 << std::endl;
      //        Rcpp::Rcout << "f1:                                " << std::flush << f1 << std::endl;
      //        Rcpp::Rcout << "f2:                                " << std::flush << f2 << std::endl;
      
      
      
      // then the minimum is to the left of n2
      // let n2 be the new upper bound
      // let n1 be the new upper test point
      
      
      upper_bound = g2;
      
      g2=g1;
      
      f2=f1;
      
      g1 = upper_bound - golden_ratio*(upper_bound - lower_bound);
      
      f1 = get_n(g1,trace_const, lambda_star, epsilon1,  epsilon_converge,gamma_star_lower,mu_const,beta_const);
      
      gap=upper_bound-lower_bound;
    }
    
    
    else{
      
      
      //              Rcpp::Rcout << "f1>f2:  "  << std::endl;
      
      //              Rcpp::Rcout << "g1:                                " << std::flush << g1 << std::endl;
      //              Rcpp::Rcout << "g2:                                " << std::flush << g2 << std::endl;
      //              Rcpp::Rcout << "f1:                                " << std::flush << f1 << std::endl;
      //              Rcpp::Rcout << "f2:                                " << std::flush << f2 << std::endl;
      
      lower_bound = g1;
      
      g1=g2;
      
      f1=f2;
      
      g2 = lower_bound + golden_ratio*(upper_bound - lower_bound);
      
      f2 = get_n(g2,trace_const, lambda_star, epsilon1,  epsilon_converge,gamma_star_lower,mu_const,beta_const);
      
      gap=upper_bound-lower_bound;
      
    }
    
    
    
  }
  
  
  gammastar=  (g1+g2)/2;
  
  min = get_n(gammastar,trace_const, lambda_star, epsilon1,  epsilon_converge,gamma_star_lower,mu_const,beta_const);
  double rstar = get_n(gammastar,trace_const, lambda_star, epsilon1,  epsilon_converge,gamma_star_lower,mu_const,beta_const,1);
  
  // Return Final Estimate of gammastar
  
  Rcpp::List outlist=Rcpp::List::create(
    Rcpp::Named("gammastar")=gammastar,
    Rcpp::Named("nstar")=min,
    Rcpp::Named("rstar")=rstar);  
  
  
  return(outlist);
  
  
}



arma::mat Mat_pow(arma::mat A, double k){
  
  int l1=A.n_rows;
  //  int l2=A.n_cols;
  
  arma::vec eigval;
  arma::mat eigvec;
  
  eig_sym(eigval, eigvec, A);
  arma::mat eigvecb=eigvec.t();
  
  for(int i=0;i<l1;i++)    eigvecb.row(i)=eigvecb.row(i)*pow(eigval(i),k);
  arma::mat B=eigvec*eigvecb;
  
  
  return(B);
}




Rcpp::List set_nstar(NumericMatrix x, NumericMatrix P, NumericMatrix P_0,
                     arma::vec mu_0, arma::vec mu_star,arma::vec beta_star,arma::vec beta_star2,double epsilon_converge,
                     NumericMatrix PD){
  
  
  // Temporary
  
  //    double m_D=1; // Scale to apply to P2 in order to temporarily generate P_D
  
  int l1=x.ncol();
  
  arma::mat P2(P.begin(),P.nrow(),P.ncol(),false);
  arma::mat x2(x.begin(), x.nrow(), x.ncol(), false);
  arma::mat P_0b(P_0.begin(), P_0.nrow(), P_0.ncol(), false);
  arma::mat PX=P2*x2;
  arma::mat XTPX=x2.t()*PX;
  
  
  // Temporary
  
  arma::mat P_D(PD.begin(),PD.nrow(),PD.ncol(),false);
  //    arma::mat P_D=m_D*P2;   // Temporarily set P_D equal to m_D*P2
  arma::mat V_temp=inv_sympd(P2+P_D);
  
  arma::mat W=x2.t()*P2*inv_sympd(P2+P_D)*PX;
  
  // take square root of P
  
  //    arma::mat P_1_2=Mat_pow(P2,0.5);
  
  
  
  arma::mat P_Inner=P_0b+XTPX;
  arma::vec mu_diff=mu_0-mu_star;
  //    arma::vec mu_const=0.5*mu_diff.t()*XTPX*mu_diff; // Old Calculation
  arma::vec mu_const=0.5*mu_diff.t()*W*mu_diff; // Revised Calculation
  
  arma::vec beta_diff=beta_star2-beta_star;
  
  //    arma::vec beta_const=0.5*beta_diff.t()*PX*inv_sympd(XTPX)*PX.t()*beta_diff; // Old Calculation
  arma::vec beta_const=0.5*beta_diff.t()*PX*inv_sympd(W)*PX.t()*beta_diff; // New Calculation
  
  
  arma::vec eigval;
  arma::mat eigvec;
  
  
  
  eig_sym(eigval, eigvec, XTPX);
  
  arma::mat eigvec2=eigvec.t();
  
  for(int i=0;i<l1;i++)    eigvec2.row(i)=eigvec2.row(i)/sqrt(eigval(i));
  
  arma::mat InvXTPX_1_2=eigvec*eigvec2;
  
  
  arma::mat InvW_1_2=Mat_pow(W,-0.5);
  
  
  // Calculations of lambda_star should change if Data Precision is bounded from below
  //  Need Matrix P_D so that 
  //    arma::mat P_BB=P2+P_D;
  
  
  // Old Calculations
  
  //    arma::mat P_AA=InvXTPX_1_2*P_Inner*InvXTPX_1_2;
  //    arma::mat P_AB=-InvXTPX_1_2*PX.t();
  //    arma::mat P_BA=P_AB.t();
  //    arma::mat P_BB=P2;
  //    arma::mat P2_AB=inv_sympd(P_AA)*P_AB*inv_sympd(P_BB)*P_BA;
  
  //    Revised Calculations
  
  
  arma::mat P_AA=InvW_1_2*P_Inner*InvW_1_2;
  arma::mat P_AB=-InvW_1_2*PX.t();
  arma::mat P_BA=P_AB.t();
  arma::mat P_BB=P2+P_D;
  arma::mat P2_AB=inv_sympd(P_AA)*P_AB*inv_sympd(P_BB)*P_BA;
  
  
  //      W.print("W:");
  //      InvW_1_2.print("InvW_1_2:");
  //      P_AA.print("P_AA:");
  //      P_D.print("P_D:");
  //      P2_AB.print("P2_AB:");
  
  arma::vec  eigen_out=eig_sym(P2_AB.t()*P2_AB) ;
  
  double lambda_star=sqrt(eigen_out(l1-1));
  
  arma::mat P_Initial=P_0b+XTPX;    
  
  // P_Lower calculations should change if Data Precision is bounded from below
  //  Need Matrix P_D so that 
  //      arma::mat P_Lower=P_Initial-PX.t()*inv_sympd(P2+P_D+PX*inv_sympd(P_Initial)*PX.t())*PX;
  
  
  arma::mat P_Upper=P_0b+XTPX;    
  
  // Old Calculation
  //    arma::mat P_Lower=P_Initial-PX.t()*inv_sympd(P2+PX*inv_sympd(P_Initial)*PX.t())*PX;
  
  //   New Calculation
  
  arma::mat P_Lower=P_Initial-PX.t()*inv_sympd(P2+P_D+PX*inv_sympd(P_Initial)*PX.t())*PX;
  
  double  det_P_Upper=det(P_Upper);
  double  det_P_Lower=det(P_Lower);
  double  epsilon1=sqrt(det_P_Lower/det_P_Upper);
  
  
  //    double  det_XTPX=det(XTPX);
  //    double  qc1=det_P_Upper/det_XTPX;
  //    double  qc2=det_P_Lower/det_XTPX;
  
  //    Rcpp::Rcout << "det_XTPX:     " << std::flush << det_XTPX << std::endl;
  //    Rcpp::Rcout << "qc1:     " << std::flush << qc1 << std::endl;
  //    Rcpp::Rcout << "qc2:     " << std::flush << qc2 << std::endl;
  //    Rcpp::Rcout << "det_P_Upper:     " << std::flush << det_P_Upper << std::endl;
  //    Rcpp::Rcout << "det_P_Lower:     " << std::flush << det_P_Lower << std::endl;
  //    Rcpp::Rcout << "epsilon1:     " << std::flush << epsilon1 << std::endl;
  
  
  // Verify Calculation for W_1_2 
  // needs to be edited if P_D is added to calculation 
  
  
  // Old calculation
  
  //    eigvec2=eigvec.t();
  
  //    for(int i=0;i<l1;i++)    eigvec2.row(i)=eigvec2.row(i)*sqrt(eigval(i));
  
  //    arma::mat W_1_2=eigvec*eigvec2;
  
  arma::mat W_1_2=Mat_pow(W,0.5);
  
  
  double trace_const = trace(W_1_2*inv_sympd(P_Lower)*W_1_2);
  double gamma_star_lower=trace_const/(1-lambda_star);
  
  
  //    Rcpp::Rcout << "############  Testing 2.1" << std::endl;
  
  // Initialize Golden Section Search
  
  //    double epsilon_converge=0.01;
  
  
  Rcpp::List golden_out=golden_n(trace_const, lambda_star,
                                 epsilon1,  epsilon_converge,gamma_star_lower,mu_const(0),beta_const(0));
  
  //    Rcpp::Rcout << "############  Testing 2.2" << std::endl;
  
  
  NumericVector temp=golden_out(0);
  double gammastar=temp(0);
  temp=golden_out(1);
  double nstar=temp(0);
  temp=golden_out(2);
  double rstar1=temp(0);
  
  
  
  double t_star=sqrt(beta_const(0)/gammastar);
  double gammastar2=(1+t_star)*gammastar;
  double trace_const2=(1+t_star)*trace_const;
  //    double mu_const2=(1+t_star)*mu_const(0);
  double beta_const2=0;
  if(t_star>0) beta_const2=beta_const(0)/t_star;
  
  
  double alpha_out=(1+gammastar2)/(1+trace_const2+lambda_star*gammastar2);
  double U_out=(1+trace_const2+2*lambda_star*gammastar2);
  double epsilon=exp(log(epsilon1)-beta_const2-gammastar2);
  
  
  
  double tau=1+2*((1/(1-sqrt(lambda_star))-1));
  
  //    Rcpp::Rcout << " "  << std::scientific   << std::endl ;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "trace_constant (Drift Condition):     " << std::flush << trace_const << std::endl;
  Rcpp::Rcout << "lambdastar (Drift Condition):         " << std::flush << lambda_star << std::endl;
  Rcpp::Rcout << "epsilon1 (Minorization Condition):    " << std::flush << epsilon1 << std::endl;
  Rcpp::Rcout << "gammastar_lower:                      " << std::flush << gamma_star_lower << std::endl;
  Rcpp::Rcout << "mu_constant:                          " << std::flush << mu_const(0) << std::endl;
  Rcpp::Rcout << "beta_constant:                        " << std::flush << beta_const(0) << std::endl;
  Rcpp::Rcout << "gammastar:                            " << std::flush << gammastar << std::endl;
  Rcpp::Rcout << "tstar:                                " << std::flush << t_star << std::endl;
  Rcpp::Rcout << ""  << std::scientific    ;
  Rcpp::Rcout << "epsilonstar:                          " << std::flush << epsilon << std::endl;
  Rcpp::Rcout << ""  << std::fixed    ;
  Rcpp::Rcout << "alpha:                                " << std::flush << alpha_out << std::endl;
  Rcpp::Rcout << "U:                                    " << std::flush << U_out << std::endl;
  Rcpp::Rcout << "rstar:                                " << std::flush << rstar1 << std::endl;
  Rcpp::Rcout << "nstar:                                " << std::flush << nstar << std::endl;
  Rcpp::Rcout << "sample size multiplier:               " << std::flush << tau << std::endl;
  Rcpp::Rcout << " " << std::endl;
  
  
  Rcpp::List simconstants=Rcpp::List::create(
    Rcpp::Named("trace_const")=trace_const,
    Rcpp::Named("lambda_star")=lambda_star,
    Rcpp::Named("epsilon1")=epsilon1,
    Rcpp::Named("gamma_star_lower")=gamma_star_lower, 
    Rcpp::Named("mu_constant")=mu_const(0),
    Rcpp::Named("beta_constant")=beta_const(0),
    Rcpp::Named("gammastar")=gammastar,
    Rcpp::Named("tstar")=t_star,
    Rcpp::Named("epsilonstar")=epsilon,
    Rcpp::Named("rstar")=rstar1,
    Rcpp::Named("nstar")=nstar
  );  
  
  
  Rcpp::List outlist=Rcpp::List::create(
    Rcpp::Named("nstar")=nstar,
    Rcpp::Named("tau")=tau,
    Rcpp::Named("simconstants")=simconstants
  );
  
  
  
  
  return(outlist);
  
  
  //    return(nstar);
  
}



Rcpp::List set_beta_const(arma::mat x2,arma::vec offset2b,arma::vec mu_star2,
                          int n,NumericVector y,NumericMatrix xtemp, 
                          NumericVector mutemp,NumericMatrix P,NumericVector offset2,NumericVector wt,
                          double dispersion,Rcpp::List
                            famfunc, Function f1,Function f2,Function f3,NumericMatrix betatemp,
                            NumericMatrix x,
                            NumericVector mu,
                            NumericMatrix P_0,
                            NumericVector offset3,
                            NumericVector wt3,
                            std::string family="binomial",
                            std::string link="logit",
                            int Gridtype=2){
  
  int l1=x2.n_cols;
  int l2=x2.n_rows;
  
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function asVec("as.vector");
  Rcpp::List out1;
  Rcpp::List out2;
  NumericMatrix betastarout(1000,l2);
  NumericVector temp(1);
  NumericVector b3(l2);
  arma::vec mu_star3(b3.begin(), l1);
  NumericMatrix Ptilde=clone(P);
  Ptilde(0,0)=2*P(0,0);
  
  
  arma::vec mutemp2=x2*mu_star2;
  arma::vec alpha2=mutemp2+offset2b;
  
  
  //      mutemp2.print("beta_star - inside function");
  
  //   Rcpp::Rcout << "mutemp2: Input         " << std::flush << mutemp2 << std::endl;
  
  
  
  for(int i=0;i<1000;i++){
    
    Rcpp::checkUserInterrupt();
    
    progress_bar2(0, 999);
    
    if(i==0) {    Rcpp::Rcout << "Running simulation for betastar:" << std::endl;}
    
    progress_bar2(i, 999);
    
    if(i==999) {Rcpp::Rcout << "" << std::endl;}
    
    for(int j=0;j<l2;j++){
      
      out1=glmbsim_NGauss2_cpp(1,asVec(y[j]),xtemp,
                               asMat(mutemp[j]),P,
                               asVec(offset2[j]),
                               asVec(wt[j]),
                               dispersion,
                               famfunc,f1,f2,f3,asMat(betatemp(j-1,0)),
                               family=family,
                               link=link,
                               Gridtype=Gridtype);
      temp=out1(0);
      betatemp(j,0)=temp(0);
      betastarout(i,j)=temp(0); 
      
    }
  }
  
  arma::vec beta_star=mutemp2;
  
  for(int j=0;j<l2;j++){  beta_star(j)=mean(betastarout(_,j));
  }
  
  //beta_star.print("beta_star - inside function");
  
  Rcpp::Rcout.precision(5);
  
  out2=glmbsim_Gauss_cpp(1,betatemp,x,
                         mu,P_0,offset3
                           ,wt3,
                           1/P(0,0),
                           famfunc,f1,f2,f3,
                           mu);   
  
  b3=out2(1);    
  
  // Simulate for bstar2
  
  arma::vec mu_star4=(mu_star2+mu_star3)/2;    
  mutemp2=x2*mu_star4;
  alpha2=mutemp2+offset2b;
  
  
  for(int i=0;i<1000;i++){
    
    Rcpp::checkUserInterrupt();
    
    progress_bar2(0, 999);
    
    if(i==0) {    Rcpp::Rcout << "Running simulation for betastar2:" << std::endl;}
    
    progress_bar2(i, 999);
    
    if(i==999) {Rcpp::Rcout << "" << std::endl;}
    
    for(int j=0;j<l2;j++){
      
      out1=glmbsim_NGauss2_cpp(1,asVec(y[j]),xtemp,
                               asMat(mutemp[j]),Ptilde,
                               asVec(offset2[j]),
                               asVec(wt[j]),
                               dispersion,
                               famfunc,f1,f2,f3,asMat(betatemp(j-1,0)),
                               family=family,
                               link=link,
                               Gridtype=Gridtype);
      temp=out1(0);
      betastarout(i,j)=temp(0); 
    }
    
  }
  
  arma::vec beta_star2=1*beta_star+0;
  for(int j=0;j<l2;j++){  beta_star2(j)=mean(betastarout(_,j));}
  
  //beta_star2.print("beta_star2 - Inside function");
  
  
  Rcpp::List beta_stars=Rcpp::List::create(Rcpp::Named("beta_star")=beta_star,
                                           Rcpp::Named("beta_star2")=beta_star2);
  
  //      double beta_const=0;
  
  return(beta_stars);
  
}




void Set_PD(const arma::vec& b2,const NumericVector& y, 
            const arma::mat& alpha2,
            const int& l1,const arma::mat& P2, const arma::mat& x2,
            const NumericVector& wt,const NumericVector& xbtemp,arma::colvec& xbtemp2,
            arma::mat& Pout2,
            std::string family="binomial",
            std::string link="logit"){
  
  // Initial Setup - Adjust to control for offset and X
  
  //  Pout2.print("Pout2: Beginning of Set_PD FUnction");
  
  
  int j;
  
  /////////////////// binomial - logit /////////////////////////////
  
  
  if(family=="binomial"&& link=="logit")
  {
    
    
    for(j=0;j<l1;j++){
      
      //      xbtemp2(j)=exp(-alpha2(j)-x2(j,1)*b2(j));
      xbtemp2(j)=exp(-alpha2(j)-b2(j));
      xbtemp2(j)=1/(1+xbtemp2(j));      
      Pout2(j,j)=wt(j)*xbtemp2(j)*(1-xbtemp2(j));
    } 
    xbtemp2.print("xbtemp2: End of Logit");
    
    //        Pout2.print("Pout2: End of Logit");
    
  }
  
  /////////////////// binomial - probit /////////////////////////////
  
  double d1;
  double p1;
  double p2;
  
  
  if(family=="binomial"&& link=="probit")
  {
    
    for(j=0;j<l1;j++){  
      //      xbtemp2(j)=alpha2(j)+x2(j,1)*b2(j);
      xbtemp2(j)=alpha2(j)+b2(j);
      d1=R::dnorm(xbtemp2(j),0,1,FALSE);
      p1=R::pnorm(xbtemp2(j),0,1,TRUE,FALSE);
      p2=R::pnorm(xbtemp2(j),0,1,FALSE,FALSE);
      Pout2(j,j)=wt(j)*d1*((y[j]*(d1+xbtemp2(j)*p1)/(p1*p1))+((1-y[j])*(d1-xbtemp2(j)*p2)/(p2*p2))   );
    }
  }
  
  
  /////////////////// binomial - cloglog /////////////////////////////
  
  //  Pout2.print("Pout2: End of Set_PD FUnction");
  
  double exb;
  double atemp;
  
  if(family=="binomial" && link=="cloglog")
  {
    
    
    
    for(int j=0;j<l1;j++){
      
      exb=exp(alpha2(j)+b2(j));
      //    p1(j)=1-exp(-exb(j));
      //    xb(j)=1-exp(-exb(j));
      //    p2(j)=exp(-exb(j));
      
      atemp=exp(exb)-1;
      
      //    xrow2=x2.row(j);
      Pout2(j,j)=wt(j)*(1-y(j))*exb
        +wt(j)*y(j)*(exb*exb*exp(exb)/(atemp*atemp ))
        -wt(j)*y(j)*(exb/atemp);
        
    }
  }
  
  
  /////////////////// poisson /////////////////////////////
  
  if(family=="poisson" )
  {
    
    for(int j=0;j<l1;j++){
      
      exb=exp(alpha2(j)+b2(j));
      Pout2(j,j)=wt(j)*exb;
      
    }
  }
  
  
}




