// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include "Set_Grid.h"
#include <math.h>


using namespace Rcpp;


// [[Rcpp::export(".rindep_norm_gamma_reg_std_cpp")]]

Rcpp::List  rindep_norm_gamma_reg_std_cpp(int n,NumericVector y,NumericMatrix x,
                                          NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,
                                          Function f2,Rcpp::List  Envelope,Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, int progbar=1)
{
  //RNGScope scope;

  int l1 = mu.nrow();
  //  int l2=pow(3,l1);
  
  std::string family2 = Rcpp::as<std::string>(family);
  std::string link2 = Rcpp::as<std::string>(link);  
  
  
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(n,l1);
  //double a1=0;
  double a2=0;
  double U=0;
  double test=0;
//  double test_int=0;
//  double test_data=0;
  double U2=0;
  
  //out(0,0)=1;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
  NumericMatrix cbars=Envelope["cbars"];
  

//  NumericMatrix cbars_int=Envelope["cbars_int"];
  NumericVector LLconst=Envelope["LLconst"]; 
//  NumericVector LLconst_int=Envelope["LLconst_int"]; 

  NumericVector outtemp=out(0,_);
  arma::rowvec outtemp2(outtemp.begin(),l1,false);
  NumericVector cbartemp=cbars(0,_);
  arma::rowvec cbartemp2(cbartemp.begin(),l1,false);

  // Add Corresponding prior terms here
  
//  NumericVector cbartemp_int=cbars_int(0,_);
//  arma::rowvec cbartemp_int2(cbartemp_int.begin(),l1,false);
  

    NumericMatrix testtemp(1,1);
  arma::mat testtemp2(testtemp.begin(),1,1,false);

  NumericMatrix testtemp_int(1,1);
  arma::mat testtemp_int2(testtemp_int.begin(),1,1,false);
  
  
    NumericMatrix btemp(l1,1);
  arma::mat btemp2(btemp.begin(),l1,1,false); 
  NumericVector testll(1);
  NumericVector testll_data(1);


    
  if(progbar==1){ Rcpp::Rcout << "Starting Simulation:" << std::endl;  };

    Rcpp::checkUserInterrupt();
    if(progbar==1){
//      progress_bar2(i, n-1);
//      if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
    }
    

    //a1=0;

      U=R::runif(0.0, 1.0);
      a2=0;
      J(0)=0;    
      while(a2==0){
        if(U<=PLSD(J(0))) a2=1;
        if(U>PLSD(J(0))){ 
          U=U-PLSD(J(0));
          J(0)=J(0)+1;
          
        }
        //a2=1; 
      }
      

      for(int j=0;j<l1;j++){  
        
        // Switch to using thetabars here
        
        out(0,j)=ctrnorm_cpp(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);    
        
        
      }


            
      // cbars_int2 holds the intercept part
      
      outtemp=out(0,_);   // Hopefully this doew not break link to outtemp2
      cbartemp=cbars(J(0),_);
      testtemp2=outtemp2 * trans(cbartemp2);
//      testtemp_int2=outtemp2 * trans(cbartemp_int2);
      U2=R::runif(0.0, 1.0);
      btemp2=trans(outtemp2);    
      
      // Need to modify to call correct f2 function based on family and link
      
      
      if(family2=="Gamma"){  
        testll=f2_gamma(btemp,y, x,mu,P,alpha,wt,0);
      }
      
      if(family2=="gaussian"){  
        testll=f2_gaussian(btemp,y, x,mu,P,alpha,wt);
        testll_data=f1_gaussian(btemp,y, x,alpha,wt);
      }
      

      test=LLconst(J(0))+testtemp(0,0)-testll(0);


//      test_int=LLconst_int(J(0))+testtemp_int(0,0)-(testll(0)-testll_data(0));
//      test_data=test-test_int;


  return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws,Rcpp::Named("test")=test,Rcpp::Named("J")=J,
                                        Rcpp::Named("log_U2")=log(U2)
//                                      ,  Rcpp::Named("test_int")=test_int,
//                                        Rcpp::Named("test_data")=test_data
                              );      
  
}


double p_inv_gamma(double dispersion,double shape,double rate){
  
  return(1- R::pgamma(1/dispersion,shape,1/rate,TRUE,FALSE));
}

double  q_inv_gamma(double p,double shape,double rate,double disp_upper,double disp_lower){
  double p_upp=p_inv_gamma(disp_upper,shape=shape,rate=rate);
  double p_low=p_inv_gamma(disp_lower,shape=shape,rate=rate);
  double p1=p_low+p*(p_upp-p_low);
  double p2=1-p1;
  return(1/ R::qgamma(p2,shape,1/rate,TRUE,FALSE));
}

double r_invgamma(double shape,double rate,double disp_upper,double disp_lower){
  double p= R::runif(0,1);
  return(q_inv_gamma(p=p,shape=shape,rate=rate,disp_upper=disp_upper,disp_lower));
}




// [[Rcpp::export(".rindep_norm_gamma_reg_std_V2_cpp")]]

Rcpp::List  rindep_norm_gamma_reg_std_v2_cpp(int n,NumericVector y,NumericMatrix x,
NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,
Function f2,Rcpp::List  Envelope,
Rcpp::List  gamma_list,
Rcpp::List  UB_list,
Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, int progbar=1)
{

  int l1 = mu.nrow();
  
  
  // Get various inputs frm the provided lists
  
  double shape3 =gamma_list["shape3"];
  double rate2 =gamma_list["rate2"];
  double disp_upper =gamma_list["disp_upper"];
  double disp_lower =gamma_list["disp_lower"];
  double RSS_ML =UB_list["RSS_ML"];
  double max_New_LL_UB =UB_list["max_New_LL_UB"];
  double max_LL_log_disp =UB_list["max_LL_log_disp"];
  double lm_log1 =UB_list["lm_log1"];
  double lm_log2 =UB_list["lm_log2"];
  NumericVector log_P_diff =UB_list["log_P_diff"];
  

  NumericVector iters_out(n);
  NumericVector disp_out(n);
  NumericVector weight_out(n);
  NumericMatrix beta_out(n,l1);
  double dispersion;
  
  
  int a1=0;
  //Rcpp::Function r_invgamma("r_invgamma");
  
  
  
  for(int i=0;i<n;i++){
    
   a1=0;
    
    while(a1==0){
      
//      dispersion=r_invgamma(_["shape"]=shape3,_["rate"]=rate2,_["disp_upper"]=disp_upper,
//                         _["disp_lower"]=disp_lower);
      
      a1=1;
    }  
    
  }
  
  // Temporarily just return non-sense constants equal to all 1
  
  return Rcpp::List::create(Rcpp::Named("out")=1,Rcpp::Named("draws")=1,
                            Rcpp::Named("test")=1,Rcpp::Named("J")=1,
                                        Rcpp::Named("log_U2")=1);  
  
  
  //  int l2=pow(3,l1);
  
  std::string family2 = Rcpp::as<std::string>(family);
  std::string link2 = Rcpp::as<std::string>(link);  
  
  
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(n,l1);
  //double a1=0;
  double a2=0;
  double U=0;
  double test=0;
  //  double test_int=0;
  //  double test_data=0;
  double U2=0;
  
  //out(0,0)=1;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
  NumericMatrix cbars=Envelope["cbars"];
  
  
  //  NumericMatrix cbars_int=Envelope["cbars_int"];
  NumericVector LLconst=Envelope["LLconst"]; 
  //  NumericVector LLconst_int=Envelope["LLconst_int"]; 
  
  NumericVector outtemp=out(0,_);
  arma::rowvec outtemp2(outtemp.begin(),l1,false);
  NumericVector cbartemp=cbars(0,_);
  arma::rowvec cbartemp2(cbartemp.begin(),l1,false);
  
  // Add Corresponding prior terms here
  
  //  NumericVector cbartemp_int=cbars_int(0,_);
  //  arma::rowvec cbartemp_int2(cbartemp_int.begin(),l1,false);
  
  
  NumericMatrix testtemp(1,1);
  arma::mat testtemp2(testtemp.begin(),1,1,false);
  
  NumericMatrix testtemp_int(1,1);
  arma::mat testtemp_int2(testtemp_int.begin(),1,1,false);
  
  
  NumericMatrix btemp(l1,1);
  arma::mat btemp2(btemp.begin(),l1,1,false); 
  NumericVector testll(1);
  NumericVector testll_data(1);
  
  
  
  if(progbar==1){ Rcpp::Rcout << "Starting Simulation:" << std::endl;  };
  
  Rcpp::checkUserInterrupt();
  if(progbar==1){
    //      progress_bar2(i, n-1);
    //      if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
  }
  
  
  //a1=0;
  
  U=R::runif(0.0, 1.0);
  a2=0;
  J(0)=0;    
  while(a2==0){
    if(U<=PLSD(J(0))) a2=1;
    if(U>PLSD(J(0))){ 
      U=U-PLSD(J(0));
      J(0)=J(0)+1;
      
    }
    //a2=1; 
  }
  
  
  for(int j=0;j<l1;j++){  
    
    // Switch to using thetabars here
    
    out(0,j)=ctrnorm_cpp(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);    
    
    
  }
  
  
  
  // cbars_int2 holds the intercept part
  
  outtemp=out(0,_);   // Hopefully this doew not break link to outtemp2
  cbartemp=cbars(J(0),_);
  testtemp2=outtemp2 * trans(cbartemp2);
  //      testtemp_int2=outtemp2 * trans(cbartemp_int2);
  U2=R::runif(0.0, 1.0);
  btemp2=trans(outtemp2);    
  
  // Need to modify to call correct f2 function based on family and link
  
  
  if(family2=="Gamma"){  
    testll=f2_gamma(btemp,y, x,mu,P,alpha,wt,0);
  }
  
  if(family2=="gaussian"){  
    testll=f2_gaussian(btemp,y, x,mu,P,alpha,wt);
    testll_data=f1_gaussian(btemp,y, x,alpha,wt);
  }
  
  
  test=LLconst(J(0))+testtemp(0,0)-testll(0);
  
  
  //      test_int=LLconst_int(J(0))+testtemp_int(0,0)-(testll(0)-testll_data(0));
  //      test_data=test-test_int;
  
  
  return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws,Rcpp::Named("test")=test,Rcpp::Named("J")=J,
                                        Rcpp::Named("log_U2")=log(U2)
                              //                                      ,  Rcpp::Named("test_int")=test_int,
                              //                                        Rcpp::Named("test_data")=test_data
  );      
  
}

