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


void progress_bar3(double x, double N)
{
  // how wide you want the progress meter to be
  int totaldotz=40;
  double fraction = x / N;
  // part of the progressmeter that's already "full"
  int dotz = round(fraction * totaldotz);
  
  Rcpp::Rcout.precision(3);
  Rcout << "\r                                                                 " << std::flush ;
  Rcout << "\r" << std::flush ;
  Rcout << std::fixed << fraction*100 << std::flush ;
  Rcout << "% [" << std::flush ;
  int ii=0;
  for ( ; ii < dotz;ii++) {
    Rcout << "=" << std::flush ;
  }
  // remaining part (spaces)
  for ( ; ii < totaldotz;ii++) {
    Rcout << " " << std::flush ;
  }
  // and back to line begin 
  
  Rcout << "]" << std::flush ;
  
  // and back to line begin 
  
  Rcout << "\r" << std::flush ;
  
}




// [[Rcpp::export(".rindep_norm_gamma_reg_std_cpp")]]

Rcpp::List  rindep_norm_gamma_reg_std_cpp(int n,NumericVector y,NumericMatrix x,
                                          NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,
                                          Function f2,Rcpp::List  Envelope,Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, int progbar=1)
{

  int l1 = mu.nrow();

  std::string family2 = Rcpp::as<std::string>(family);
  std::string link2 = Rcpp::as<std::string>(link);  
  
  
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(n,l1);
  double a2=0;
  double U=0;
//  double test=0;
  double U2=0;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
  NumericMatrix cbars=Envelope["cbars"];
  
      U=R::runif(0.0, 1.0);
      a2=0;
      J(0)=0;    
      while(a2==0){
        if(U<=PLSD(J(0))) a2=1;
        if(U>PLSD(J(0))){ 
          U=U-PLSD(J(0));
          J(0)=J(0)+1;
          
        }
      }
      
      for(int j=0;j<l1;j++){  out(0,j)=ctrnorm_cpp(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);          }
      U2=R::runif(0.0, 1.0);


  return Rcpp::List::create(Rcpp::Named("out")=out,
                                        Rcpp::Named("J")=J,
                                        Rcpp::Named("log_U2")=log(U2)
                              );      
  
}


double p_inv_gamma(double dispersion,double shape,double rate){
  
  return(1- R::pgamma(1/dispersion,shape,1/rate,TRUE,FALSE));
}

double  q_inv_gamma(double p,double shape,double rate,double disp_upper,double disp_lower){
  double p_upp=p_inv_gamma(disp_upper,shape,rate);
  double p_low=p_inv_gamma(disp_lower,shape,rate);
  double p1=p_low+p*(p_upp-p_low);
  double p2=1-p1;
  return(1/ R::qgamma(p2,shape,1/rate,TRUE,FALSE));
}

double r_invgamma(double shape,double rate,double disp_upper,double disp_lower){
  double p= R::runif(0,1);
  return(q_inv_gamma(p,shape,rate,disp_upper,disp_lower));
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
  int l2 = x.nrow();
  
  
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
  NumericMatrix cbars=Envelope["cbars"];
  

  NumericVector iters_out(n);
  NumericVector disp_out(n);
  NumericVector weight_out(n);
  NumericMatrix beta_out(n,l1);
  double dispersion;
  NumericVector wt2(l1);
  
  
  arma::vec wt1b(wt.begin(), x.nrow());
//  arma::vec wt2b(wt2.begin(), x.nrow());
  

  NumericMatrix cbarst(cbars.ncol(),cbars.nrow());
  NumericMatrix thetabars(cbars.nrow(),cbars.ncol());
  NumericVector New_LL(cbars.nrow());
  
  
  
  arma::mat cbarsb(cbars.begin(), cbars.nrow(), cbars.ncol(), false);
  arma::mat cbarstb(cbarst.begin(), cbarst.nrow(), cbarst.ncol(), false);

  arma::mat thetabarsb(thetabars.begin(), thetabars.nrow(), thetabars.ncol(), false);
  cbarstb=trans(cbarsb);
  
  arma::vec y2(y.begin(),l2);
  arma::vec alpha2(alpha.begin(),l2);
  arma::mat x2(x.begin(),l2,l1);
  arma::mat P2(P.begin(),l1,l1);
//  arma::vec log_P_diff2(log_P_diff.begin(),cbars.nrow());
  
  
  double UB1;
  double UB2;
  double UB3A;
  double UB3B;
  double max_New_LL;
  double New_LL_log_disp;
    
  int a1=0;
  double test=0;
  //Rcpp::Function r_invgamma("r_invgamma");
  
  
  
  for(int i=0;i<n;i++){
    
   a1=0;
    iters_out[i]=1;  
    while(a1==0){
      
      dispersion=r_invgamma(shape3,rate2,disp_upper,disp_lower);

      //Rcpp::Rcout << "wt original" << std::flush << wt << std::endl;
      
      //wt1b.print("original weight");
      //wt2b=wt1b/dispersion;
      wt2=wt/dispersion;
      
//      wt2b.print("updated weight");
//      Rcpp::Rcout << "wt2 -should match wt2b" << std::flush << wt2 << std::endl;
      
      arma::mat theta =Inv_f3_gaussian(transpose(cbars), y,x, mu, P, alpha, wt2);  
//      theta.print("new thetabars");
      thetabarsb=theta;
      
//      Rcpp::Rcout << "thetabars_new - actual " << std::flush << thetabars << std::endl;
      
      NumericVector LL_New=-f2_gaussian(transpose(thetabars),  y, x, mu, P, alpha, wt2);  
      
      Rcpp::List  sim_list=rindep_norm_gamma_reg_std_cpp(1,y,x,mu, P,alpha,wt2,
      f2,Envelope,family,link,  progbar); 
      
      NumericVector J_out=sim_list["J"];

//      Rcpp::Rcout << "J_out " << std::flush << J_out << std::endl;
      
      double log_U2=sim_list["log_U2"];
      NumericMatrix sim_out=sim_list["out"];  
      NumericVector b_out=sim_out(0,_);
      arma::rowvec b_out2(b_out.begin(),l1,false);
      NumericVector thetabars_temp=thetabars(J_out(0),_);
      arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
      NumericVector cbars_temp=cbars(J_out(0),_);
      arma::vec  cbars_temp2(cbars_temp.begin(), l1);
      
      
//            thetabars_temp2.print("Selected thetabars");

      // Passing transpose of sim_out here might work because sim_out has only one row

      NumericVector LL_Test=-f2_gaussian(transpose(sim_out),  y, x, mu, P, alpha, wt2);  
      
      // Block 1: UB1 
      //   Same form as in fixed dispersion case but thetabar is a function of the dispersion
      //   So all components that include thetabar must now be bounded as well
      
      arma::colvec betadiff=trans(b_out2)-thetabars_temp2;
      UB1=LL_New(J_out(0)) -arma::as_scalar(trans(cbars_temp2)*betadiff);
      
      //Block 2: UB2 [RSS Term bounded by shifting it to the gamma candidate]
      
      // % is element wise multiplication
      
      arma::colvec yxbeta=(y2-alpha2-x2*thetabars_temp2)%sqrt(wt1b); 
      UB2=0.5*(1/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_ML);
      
      //  UB1=arma::as_scalar(trans(cbars_temp2)*betadiff);
      //  UB1=LL_New[J_out]-Env2$cbars[J_out,1:ncol(x)]%*%(betadiff)
    
      // Block 3: UB3A (adjusts because probabilities of components in grid are different from original grid)
      // Investigate whether changing probabilities of grid componets for proposal
      // allows us to do away with this term and to thereby improve the acceptance rate
    
        for(int j=0;j<cbars.nrow();j++){
          thetabars_temp=thetabars(j,_);
          cbars_temp=cbars(j,_);
          arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
          arma::vec  cbars_temp2(cbars_temp.begin(), l1);

          New_LL(j)=arma::as_scalar(-0.5*trans(thetabars_temp2)*P2*thetabars_temp2
                                      +trans(cbars_temp2)*thetabars_temp2);
                    
        }
    
        //arma::vec New_LL2(New_LL.begin(),cbars.nrow());
        NumericVector LL_temp=log_P_diff+ New_LL;
        max_New_LL=max(LL_temp);
        UB3A=max_New_LL-LL_temp(J_out(0));
  
  
      // Block 4: UB3B  
      
      New_LL_log_disp=lm_log1+lm_log2*log(dispersion);
      UB3B=(max_New_LL_UB-max_LL_log_disp)-(max_New_LL-New_LL_log_disp);
  
  
      test= LL_Test[0]-UB1;  // Should be all negative 

//        Rcpp::Rcout << "test1 " << std::flush << test << std::endl;

      test= LL_Test[0]-(UB1+UB2);  // Should be all negative 

//      Rcpp::Rcout << "test2 " << std::flush << test << std::endl;

      test= LL_Test[0]-(UB1+UB2+UB3A);  // Should be all negative 

      //Rcpp::Rcout << "test3 " << std::flush << test << std::endl;
      
      test= LL_Test[0]-(UB1+UB2+UB3A+UB3B);  // Should be all negative 
      
      //Rcpp::Rcout << "test3 " << std::flush << test << std::endl;
      //      P4.print("P4 after step 1");  
        //      epsilon.print("epsilon after step 1");  

        // Block 1: UB1 
        //   Same form as in fixed dispersion case but thetabar is a function of the dispersion
        //   So all components that include thetabar must now be bounded as well
        
        
      test=test-log_U2;
          
      
     disp_out[i]=dispersion;  
     beta_out(i,_)=sim_out(0,_);
     weight_out[i]=max_New_LL;
       
     
     if(test>=0) a1=1;
     else{iters_out[i]=iters_out[i]+1;}        

    //  a1=1;
    }  
    
  }
  
  // Temporarily just return non-sense constants equal to all 1
  
  return Rcpp::List::create(Rcpp::Named("beta_out")=beta_out,Rcpp::Named("disp_out")=disp_out,
                            Rcpp::Named("iters_out")=iters_out,Rcpp::Named("weight_out")=weight_out);  
  
  
  //  int l2=pow(3,l1);
  
  std::string family2 = Rcpp::as<std::string>(family);
  std::string link2 = Rcpp::as<std::string>(link);  
  
  
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(n,l1);
  //double a1=0;
  double a2=0;
  double U=0;
  //double test=0;
  //  double test_int=0;
  //  double test_data=0;
  double U2=0;
  
  //out(0,0)=1;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
//  NumericMatrix cbars=Envelope["cbars"];
  
  
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


// [[Rcpp::export(".rindep_norm_gamma_reg_std_V3_cpp")]]

Rcpp::List  rindep_norm_gamma_reg_std_v3_cpp(int n,NumericVector y,NumericMatrix x,
                                             NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,
                                             Function f2,Rcpp::List  Envelope,
                                             Rcpp::List  gamma_list,
                                             Rcpp::List  UB_list,
                                             Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, int progbar=1)
{
  
  int l1 = mu.nrow();
  int l2 = x.nrow();
  
  
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
  double lmc1 =UB_list["lmc1"];
  double lmc2 =UB_list["lmc2"];
  NumericVector lg_prob_factor =UB_list["lg_prob_factor"];
  NumericMatrix cbars=Envelope["cbars"];
  
  
  NumericVector iters_out(n);
  NumericVector disp_out(n);
  NumericVector weight_out(n);
  NumericMatrix beta_out(n,l1);
  double dispersion;
  NumericVector wt2(l1);
  
  
  arma::vec wt1b(wt.begin(), x.nrow());

  
  NumericMatrix cbarst(cbars.ncol(),cbars.nrow());
  NumericMatrix thetabars(cbars.nrow(),cbars.ncol());
  NumericMatrix thetabars_new(1,cbars.ncol());
  
  NumericVector New_LL(cbars.nrow());
  
  
  
  arma::mat cbarsb(cbars.begin(), cbars.nrow(), cbars.ncol(), false);
  arma::mat cbarstb(cbarst.begin(), cbarst.nrow(), cbarst.ncol(), false);
  
  arma::mat thetabarsb(thetabars.begin(), thetabars.nrow(), thetabars.ncol(), false);
  arma::mat thetabarsb_new(thetabars_new.begin(), thetabars_new.nrow(), thetabars_new.ncol(), false);
  cbarstb=trans(cbarsb);
  
  arma::vec y2(y.begin(),l2);
  arma::vec alpha2(alpha.begin(),l2);
  arma::mat x2(x.begin(),l2,l1);
  arma::mat P2(P.begin(),l1,l1);

  double UB1;
  double UB2;
  double UB3A;
  double UB3B;
  double New_LL_log_disp;
  
  int a1=0;
  double test=0;
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(1,l1);
  double a2=0;
  double U=0;
  double U2=0;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];

  
  if(progbar==1){ Rcpp::Rcout << "Starting Simulation:" << std::endl;  };
  
  
  
  
  for(int i=0;i<n;i++){

    Rcpp::checkUserInterrupt();
    if(progbar==1){
            progress_bar3(i, n-1);
            if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
    }
    
    
        
    a1=0;
    iters_out[i]=1;  
    while(a1==0){
      
      dispersion=r_invgamma(shape3,rate2,disp_upper,disp_lower);
      
      wt2=wt/dispersion;
      
      // Simulate from discrete distribution
      
      U=R::runif(0.0, 1.0);
      a2=0;
      J(0)=0;    
      while(a2==0){
        if(U<=PLSD(J(0))) a2=1;
        if(U>PLSD(J(0))){ 
          U=U-PLSD(J(0));
          J(0)=J(0)+1;
          
        }
      }
      

      NumericMatrix cbars_small = cbars( Range(J(0),J(0)) , Range(0,cbars.ncol()-1) );
      
      arma::mat theta2 =Inv_f3_gaussian(transpose(cbars_small), y,x, mu, P, alpha, wt2);  

      thetabarsb_new=theta2;
      NumericVector LL_New2=-f2_gaussian(transpose(thetabars_new),  y, x, mu, P, alpha, wt2);  
      

      
      for(int j=0;j<l1;j++){  out(0,j)=ctrnorm_cpp(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);          }

      
      U2=R::runif(0.0, 1.0);
      
      double log_U2=log(U2);
      NumericVector J_out=J;
      NumericVector b_out=out(0,_);
      arma::rowvec b_out2(b_out.begin(),l1,false);
      NumericVector thetabars_temp=thetabars_new(0,_); // Changed

      arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
      NumericVector cbars_temp=cbars(J_out(0),_);
      arma::vec  cbars_temp2(cbars_temp.begin(), l1);
      
      

      NumericVector LL_Test=-f2_gaussian(transpose(out),  y, x, mu, P, alpha, wt2);

      // Block 1: UB1 
      //   Same form as in fixed dispersion case but thetabar is a function of the dispersion
      //   So all components that include thetabar must now be bounded as well
      
      arma::colvec betadiff=trans(b_out2)-thetabars_temp2;
      UB1=LL_New2(0) -arma::as_scalar(trans(cbars_temp2)*betadiff);
      
      //Block 2: UB2 [RSS Term bounded by shifting it to the gamma candidate]
      

      arma::colvec yxbeta=(y2-alpha2-x2*thetabars_temp2)%sqrt(wt1b); 
      UB2=0.5*(1/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_ML);
      

      // Block 3: UB3A (adjusts because probabilities of components in grid are different from original grid)
      // Investigate whether changing probabilities of grid componets for proposal
      // allows us to do away with this term and to thereby improve the acceptance rate
      
      // This is likely time consuming part
      
        for(int j=J_out(0);j<(J_out(0)+1);j++){
        thetabars_temp=thetabars_new(0,_); // Changed
        
        
        cbars_temp=cbars(j,_);
        arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
        arma::vec  cbars_temp2(cbars_temp.begin(), l1);
        
        New_LL(j)=arma::as_scalar(-0.5*trans(thetabars_temp2)*P2*thetabars_temp2
                                    +trans(cbars_temp2)*thetabars_temp2);
        
      }
      

      // Modified UB3A 
      
      UB3A= lg_prob_factor(J_out(0))+lmc1+lmc2*dispersion-New_LL(J_out(0));
      
      // Block 4: UB3B  
      
      New_LL_log_disp=lm_log1+lm_log2*log(dispersion);
      
      UB3B=(max_New_LL_UB-max_LL_log_disp+New_LL_log_disp)-(lmc1+lmc2*dispersion);
      
      //test= LL_Test[0]-UB1;  // Should be all negative 
      
      //        Rcpp::Rcout << "test1 " << std::flush << test << std::endl;
      
      //test= LL_Test[0]-(UB1+UB2);  // Should be all negative 
      
      //      Rcpp::Rcout << "test2 " << std::flush << test << std::endl;
      
      //test= LL_Test[0]-(UB1+UB2+UB3A);  // Should be all negative 
      
      // Rcpp::Rcout << "test3 " << std::flush << test << std::endl;
      
      test= LL_Test[0]-(UB1+UB2+UB3A+UB3B);  // Should be all negative 
      
      //  Rcpp::Rcout << "test4 " << std::flush << test << std::endl;
      
      //      P4.print("P4 after step 1");  
      //      epsilon.print("epsilon after step 1");  
      
      // Block 1: UB1 
      //   Same form as in fixed dispersion case but thetabar is a function of the dispersion
      //   So all components that include thetabar must now be bounded as well
      
      
      test=test-log_U2;
      
      
      disp_out[i]=dispersion;  
      beta_out(i,_)=out(0,_);

      if(test>=0) a1=1;
      else{iters_out[i]=iters_out[i]+1;}        
      
    }  
    

  }
  
  // Temporarily just return non-sense constants equal to all 1
  
  return Rcpp::List::create(Rcpp::Named("beta_out")=beta_out,Rcpp::Named("disp_out")=disp_out,
                            Rcpp::Named("iters_out")=iters_out,Rcpp::Named("weight_out")=weight_out);  
  
  
  
}



// [[Rcpp::export(".rindep_norm_gamma_reg_std_V4_cpp")]]

Rcpp::List  rindep_norm_gamma_reg_std_v4_cpp(int n,NumericVector y,NumericMatrix x,
                                             NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,
                                             Function f2,Rcpp::List  Envelope,
                                             Rcpp::List  gamma_list,
                                             Rcpp::List  UB_list,
                                             Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, int progbar=1)
{
  
  int l1 = mu.nrow();
  int l2 = x.nrow();
  
  
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
  double lmc1 =UB_list["lmc1"];
  double lmc2 =UB_list["lmc2"];
  NumericVector lg_prob_factor =UB_list["lg_prob_factor"];
  NumericMatrix cbars=Envelope["cbars"];
  
  
  NumericVector iters_out(n);
  NumericVector disp_out(n);
  NumericVector weight_out(n);
  NumericMatrix beta_out(n,l1);
  double dispersion;
  NumericVector wt2(l1);
  
  
  arma::vec wt1b(wt.begin(), x.nrow());
  
  
  NumericMatrix cbarst(cbars.ncol(),cbars.nrow());
  NumericMatrix thetabars(cbars.nrow(),cbars.ncol());
  NumericMatrix thetabars_new(1,cbars.ncol());
  
  NumericVector New_LL(cbars.nrow());
  
  
  
  arma::mat cbarsb(cbars.begin(), cbars.nrow(), cbars.ncol(), false);
  arma::mat cbarstb(cbarst.begin(), cbarst.nrow(), cbarst.ncol(), false);
  
  arma::mat thetabarsb(thetabars.begin(), thetabars.nrow(), thetabars.ncol(), false);
  arma::mat thetabarsb_new(thetabars_new.begin(), thetabars_new.nrow(), thetabars_new.ncol(), false);
  cbarstb=trans(cbarsb);
  
  arma::vec y2(y.begin(),l2);
  arma::vec alpha2(alpha.begin(),l2);
  arma::mat x2(x.begin(),l2,l1);
  arma::mat P2(P.begin(),l1,l1);
  
  double UB1;
  double UB2;
  double UB3A;
  double UB3B;
  double New_LL_log_disp;
  
  int a1=0;
  double test=0;
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(1,l1);
  double a2=0;
  double U=0;
  double U2=0;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];

  if(progbar==1){ Rcpp::Rcout << "Starting Simulation:" << std::endl;  };
  
  
  for(int i=0;i<n;i++){

    Rcpp::checkUserInterrupt();
    if(progbar==1){
      progress_bar3(i, n-1);
      if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
    }
    
    
        
    a1=0;
    iters_out[i]=1;  
    while(a1==0){

            
      
      // Simulate from discrete distribution
      
      U=R::runif(0.0, 1.0);
      a2=0;
      J(0)=0;    
      while(a2==0){
        if(U<=PLSD(J(0))) a2=1;
        if(U>PLSD(J(0))){ 
          U=U-PLSD(J(0));
          J(0)=J(0)+1;
          
        }
      }
      
      
      // Simulate for beta
      
      for(int j=0;j<l1;j++){  out(0,j)=ctrnorm_cpp(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);          }
      
      
      // Update this to make distribution contingent on component of the grid
      
      dispersion=r_invgamma(shape3,rate2,disp_upper,disp_lower);
      
      wt2=wt/dispersion;
      NumericMatrix cbars_small = cbars( Range(J(0),J(0)) , Range(0,cbars.ncol()-1) );
      arma::mat theta2 =Inv_f3_gaussian(transpose(cbars_small), y,x, mu, P, alpha, wt2);  
      
      thetabarsb_new=theta2;
      NumericVector LL_New2=-f2_gaussian(transpose(thetabars_new),  y, x, mu, P, alpha, wt2);  
      
      
      
      
      
      U2=R::runif(0.0, 1.0);
      
      double log_U2=log(U2);
      NumericVector J_out=J;
      NumericVector b_out=out(0,_);
      arma::rowvec b_out2(b_out.begin(),l1,false);
      NumericVector thetabars_temp=thetabars_new(0,_); // Changed
      
      arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
      NumericVector cbars_temp=cbars(J_out(0),_);
      arma::vec  cbars_temp2(cbars_temp.begin(), l1);
      
      
      
      NumericVector LL_Test=-f2_gaussian(transpose(out),  y, x, mu, P, alpha, wt2);
      
      // Block 1: UB1 
      //   Same form as in fixed dispersion case but thetabar is a function of the dispersion
      //   So all components that include thetabar must now be bounded as well
      
      arma::colvec betadiff=trans(b_out2)-thetabars_temp2;
      UB1=LL_New2(0) -arma::as_scalar(trans(cbars_temp2)*betadiff);
      
      //Block 2: UB2 [RSS Term bounded by shifting it to the gamma candidate]
      
      
      arma::colvec yxbeta=(y2-alpha2-x2*thetabars_temp2)%sqrt(wt1b); 
      UB2=0.5*(1/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_ML);
      
      
      // Block 3: UB3A (adjusts because probabilities of components in grid are different from original grid)
      // Investigate whether changing probabilities of grid componets for proposal
      // allows us to do away with this term and to thereby improve the acceptance rate
      
      // This is likely time consuming part
      
      for(int j=J_out(0);j<(J_out(0)+1);j++){
        thetabars_temp=thetabars_new(0,_); // Changed
        
        
        cbars_temp=cbars(j,_);
        arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
        arma::vec  cbars_temp2(cbars_temp.begin(), l1);
        
        New_LL(j)=arma::as_scalar(-0.5*trans(thetabars_temp2)*P2*thetabars_temp2
                                    +trans(cbars_temp2)*thetabars_temp2);
        
      }
      
      
      // Modified UB3A 
      
      UB3A= lg_prob_factor(J_out(0))+lmc1+lmc2*dispersion-New_LL(J_out(0));
      
      // Block 4: UB3B  
      
      New_LL_log_disp=lm_log1+lm_log2*log(dispersion);
      
      UB3B=(max_New_LL_UB-max_LL_log_disp+New_LL_log_disp)-(lmc1+lmc2*dispersion);
      
      //test= LL_Test[0]-UB1;  // Should be all negative 
      
      //        Rcpp::Rcout << "test1 " << std::flush << test << std::endl;
      
      //test= LL_Test[0]-(UB1+UB2);  // Should be all negative 
      
      //      Rcpp::Rcout << "test2 " << std::flush << test << std::endl;
      
      //test= LL_Test[0]-(UB1+UB2+UB3A);  // Should be all negative 
      
      // Rcpp::Rcout << "test3 " << std::flush << test << std::endl;
      
      test= LL_Test[0]-(UB1+UB2+UB3A+UB3B);  // Should be all negative 
      
      //  Rcpp::Rcout << "test4 " << std::flush << test << std::endl;
      
      //      P4.print("P4 after step 1");  
      //      epsilon.print("epsilon after step 1");  
      
      // Block 1: UB1 
      //   Same form as in fixed dispersion case but thetabar is a function of the dispersion
      //   So all components that include thetabar must now be bounded as well
      
      
      test=test-log_U2;
      
      
      disp_out[i]=dispersion;  
      beta_out(i,_)=out(0,_);
      
      if(test>=0) a1=1;
      else{iters_out[i]=iters_out[i]+1;}        
      
    }  
    
    
  }
  
  // Temporarily just return non-sense constants equal to all 1
  
  return Rcpp::List::create(Rcpp::Named("beta_out")=beta_out,Rcpp::Named("disp_out")=disp_out,
                            Rcpp::Named("iters_out")=iters_out,Rcpp::Named("weight_out")=weight_out);  
  
  
  
}



// [[Rcpp::export(".rindep_norm_gamma_reg_std_V5_cpp")]]

Rcpp::List  rindep_norm_gamma_reg_std_v5_cpp(int n,NumericVector y,NumericMatrix x,
                                             NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,
                                             Function f2,Rcpp::List  Envelope,
                                             Rcpp::List  gamma_list,
                                             Rcpp::List  UB_list,
                                             Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, int progbar=1)
{
  
  int l1 = mu.nrow();
  int l2 = x.nrow();
  
  
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
  double lmc1 =UB_list["lmc1"];
  double lmc2 =UB_list["lmc2"];
  
  NumericVector shape3_vector=gamma_list["shape3_vector"];
  NumericVector lm_log1_vector=UB_list["lm_log1_vector"];
  NumericVector lm_log2_vector=UB_list["lm_log2_vector"];
  NumericVector lm_c1_vector=UB_list["lm_c1_vector"];
  NumericVector lm_c2_vector=UB_list["lm_c2_vector"];
  
  
  NumericVector lg_prob_factor =UB_list["lg_prob_factor"];
  NumericMatrix cbars=Envelope["cbars"];
  
  
  NumericVector iters_out(n);
  NumericVector disp_out(n);
  NumericVector weight_out(n);
  NumericMatrix beta_out(n,l1);
  double dispersion;
  NumericVector wt2(l1);
  
  
  arma::vec wt1b(wt.begin(), x.nrow());
  
  
  NumericMatrix cbarst(cbars.ncol(),cbars.nrow());
  NumericMatrix thetabars(cbars.nrow(),cbars.ncol());
  NumericMatrix thetabars_new(1,cbars.ncol());
  
  NumericVector New_LL(cbars.nrow());
  
  
  
  arma::mat cbarsb(cbars.begin(), cbars.nrow(), cbars.ncol(), false);
  arma::mat cbarstb(cbarst.begin(), cbarst.nrow(), cbarst.ncol(), false);
  
  arma::mat thetabarsb(thetabars.begin(), thetabars.nrow(), thetabars.ncol(), false);
  arma::mat thetabarsb_new(thetabars_new.begin(), thetabars_new.nrow(), thetabars_new.ncol(), false);
  cbarstb=trans(cbarsb);
  
  arma::vec y2(y.begin(),l2);
  arma::vec alpha2(alpha.begin(),l2);
  arma::mat x2(x.begin(),l2,l1);
  arma::mat P2(P.begin(),l1,l1);
  
  double UB1;
  double UB2;
  double UB3A;
  double UB3B;
  double New_LL_log_disp;
  
  int a1=0;
  double test=0;
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(1,l1);
  double a2=0;
  double U=0;
  double U2=0;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
  
  if(progbar==1){ Rcpp::Rcout << "Starting Simulation:" << std::endl;  };
  
  
  for(int i=0;i<n;i++){
    
    Rcpp::checkUserInterrupt();
    if(progbar==1){
      progress_bar3(i, n-1);
      if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
    }
    
    
    
    a1=0;
    iters_out[i]=1;  
    while(a1==0){
      
      
      
      // Simulate from discrete distribution
      
      U=R::runif(0.0, 1.0);
      a2=0;
      J(0)=0;    
      while(a2==0){
        if(U<=PLSD(J(0))) a2=1;
        if(U>PLSD(J(0))){ 
          U=U-PLSD(J(0));
          J(0)=J(0)+1;
          
        }
      }
      
      
      // Simulate for beta
      
      for(int j=0;j<l1;j++){  out(0,j)=ctrnorm_cpp(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);          }
      
      
      // Update this to make distribution contingent on component of the grid
      
//      dispersion=r_invgamma(shape3,rate2,disp_upper,disp_lower);
      dispersion=r_invgamma(shape3_vector(J(0)),rate2,disp_upper,disp_lower);
      
      wt2=wt/dispersion;
      NumericMatrix cbars_small = cbars( Range(J(0),J(0)) , Range(0,cbars.ncol()-1) );
      arma::mat theta2 =Inv_f3_gaussian(transpose(cbars_small), y,x, mu, P, alpha, wt2);  
      
      thetabarsb_new=theta2;
      NumericVector LL_New2=-f2_gaussian(transpose(thetabars_new),  y, x, mu, P, alpha, wt2);  
      
      
      
      
      
      U2=R::runif(0.0, 1.0);
      
      double log_U2=log(U2);
      NumericVector J_out=J;
      NumericVector b_out=out(0,_);
      arma::rowvec b_out2(b_out.begin(),l1,false);
      NumericVector thetabars_temp=thetabars_new(0,_); // Changed
      
      arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
      NumericVector cbars_temp=cbars(J_out(0),_);
      arma::vec  cbars_temp2(cbars_temp.begin(), l1);
      
      
      
      NumericVector LL_Test=-f2_gaussian(transpose(out),  y, x, mu, P, alpha, wt2);
      
      // Block 1: UB1 
      //   Same form as in fixed dispersion case but thetabar is a function of the dispersion
      //   So all components that include thetabar must now be bounded as well
      
      arma::colvec betadiff=trans(b_out2)-thetabars_temp2;
      UB1=LL_New2(0) -arma::as_scalar(trans(cbars_temp2)*betadiff);
      
      //Block 2: UB2 [RSS Term bounded by shifting it to the gamma candidate]
      
      
      arma::colvec yxbeta=(y2-alpha2-x2*thetabars_temp2)%sqrt(wt1b); 
      UB2=0.5*(1/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_ML);
      
      
      // Block 3: UB3A (adjusts because probabilities of components in grid are different from original grid)
      // Investigate whether changing probabilities of grid componets for proposal
      // allows us to do away with this term and to thereby improve the acceptance rate
      
      // This is likely time consuming part
      
      for(int j=J_out(0);j<(J_out(0)+1);j++){
        thetabars_temp=thetabars_new(0,_); // Changed
        
        
        cbars_temp=cbars(j,_);
        arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
        arma::vec  cbars_temp2(cbars_temp.begin(), l1);
        
        New_LL(j)=arma::as_scalar(-0.5*trans(thetabars_temp2)*P2*thetabars_temp2
                                    +trans(cbars_temp2)*thetabars_temp2);
        
      }
      
      
      // Modified UB3A 
      
      UB3A= lg_prob_factor(J_out(0))+lmc1+lmc2*dispersion-New_LL(J_out(0));
      
      // Block 4: UB3B  
      
      New_LL_log_disp=lm_log1+lm_log2*log(dispersion);
      
      UB3B=(max_New_LL_UB-max_LL_log_disp+New_LL_log_disp)-(lmc1+lmc2*dispersion);
      
      //test= LL_Test[0]-UB1;  // Should be all negative 
      
      //        Rcpp::Rcout << "test1 " << std::flush << test << std::endl;
      
      //test= LL_Test[0]-(UB1+UB2);  // Should be all negative 
      
      //      Rcpp::Rcout << "test2 " << std::flush << test << std::endl;
      
      //test= LL_Test[0]-(UB1+UB2+UB3A);  // Should be all negative 
      
      // Rcpp::Rcout << "test3 " << std::flush << test << std::endl;
      
      test= LL_Test[0]-(UB1+UB2+UB3A+UB3B);  // Should be all negative 
      
      //  Rcpp::Rcout << "test4 " << std::flush << test << std::endl;
      
      //      P4.print("P4 after step 1");  
      //      epsilon.print("epsilon after step 1");  
      
      // Block 1: UB1 
      //   Same form as in fixed dispersion case but thetabar is a function of the dispersion
      //   So all components that include thetabar must now be bounded as well
      
      
      test=test-log_U2;
      
      
      disp_out[i]=dispersion;  
      beta_out(i,_)=out(0,_);
      
      if(test>=0) a1=1;
      else{iters_out[i]=iters_out[i]+1;}        
      
    }  
    
    
  }
  
  // Temporarily just return non-sense constants equal to all 1
  
  return Rcpp::List::create(Rcpp::Named("beta_out")=beta_out,Rcpp::Named("disp_out")=disp_out,
                            Rcpp::Named("iters_out")=iters_out,Rcpp::Named("weight_out")=weight_out);  
  
  
  
}
