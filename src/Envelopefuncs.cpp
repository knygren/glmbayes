// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us

#include "RcppArmadillo.h"

using namespace Rcpp;


// [[Rcpp::export("Set_Grid")]]
Rcpp::List Set_Grid(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint) {

  // Get dimensions
  
  int l1=GIndex.ncol();
  int l2=GIndex.nrow();
 
  // Initialize Matrices 
  
  Rcpp::NumericMatrix Down(l2,l1);
  Rcpp::NumericMatrix Up(l2,l1);
  Rcpp::NumericMatrix lglt(l2,l1);
  Rcpp::NumericMatrix lgrt(l2,l1);
  Rcpp::NumericMatrix lgct(l2,l1);
  Rcpp::NumericMatrix logU(l2,l1);
  Rcpp::NumericMatrix logP(l2,2);

  // Main Function Code

  for(int j=0;j<l2;j++)
{
  logP(j,0)=0;
  
  for(int i=0;i<l1;i++)
{
  if(GIndex(j,i)==1){
  Down(j,i)=-INFINITY;
  Up(j,i)=Lint(0,i)+cbars(j,i);
  }
  if(GIndex(j,i)==2){
  Down(j,i)=Lint(0,i)+cbars(j,i);
  Up(j,i)=Lint(1,i)+cbars(j,i); 
  
  }
  if(GIndex(j,i)==3){
  Down(j,i)=Lint(1,i)+cbars(j,i); 
  Up(j,i)=INFINITY;
  }
  if(GIndex(j,i)==4){
    Down(j,i)=-INFINITY;
   Up(j,i)=INFINITY;
  }
}
}
  
  
  for(int i=0;i<l1;i++){
  lglt(_,i) = pnorm(Up(_,i),0.0,1.0,true,true);
  lgrt(_,i) = pnorm(Down(_,i),0.0,1.0,false,true);
  lgct(_,i) = pnorm(Up(_,i),0.0,1.0)-pnorm(Down(_,i),0.0,1.0); 
  }
  
  
  for(int j=0;j<l2;j++)
{
  for(int i=0;i<l1;i++)
{
  if(GIndex(j,i)==1){
  logU(j,i)=lglt(j,i);
  logP(j,0)=logP(j,0)+logU(j,i);
  }
  if(GIndex(j,i)==2){
  lgct(j,i)=log(lgct(j,i));
  logU(j,i)=lgct(j,i);
  logP(j,0)=logP(j,0)+logU(j,i);
  }
  if(GIndex(j,i)==3){
  logU(j,i)=lgrt(j,i);
  logP(j,0)=logP(j,0)+logU(j,i);
  }
  if(GIndex(j,i)==4){
  logU(j,i)=0;
  logP(j,0)=logP(j,0)+logU(j,i);
  }
}}
  
  // Return List

  return Rcpp::List::create(Rcpp::Named("Down")=Down
  ,Rcpp::Named("Up")=Up
  ,Rcpp::Named("lglt")=lglt
  ,Rcpp::Named("lgrt")=lgrt
  ,Rcpp::Named("lgct")=lgct
  ,Rcpp::Named("logU")=logU
  ,Rcpp::Named("logP")=logP);

}


// [[Rcpp::export]]
Rcpp::List   setlogP(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3) {

    int n = logP.nrow(), k = logP.ncol();
    int l1 =cbars.ncol();
//    int l2=cbars.nrow();
    
    arma::mat logP2(logP.begin(), n, k, false); 
    NumericVector cbartemp=cbars(0,_);  
    NumericVector G3temp=G3(0,_);  
    Rcpp::NumericMatrix LLconst(n,1);
    
    arma::colvec cbarrow(cbartemp.begin(),l1,false);
    arma::colvec G3row(G3temp.begin(),l1,false);
    
//    double v = arma::as_scalar(cbarrow.t() * cbarrow);
//    LLconst[j]<--t(as.matrix(cbars[j,1:l1]))%*%t(as.matrix(G3[j,1:l1]))+NegLL[j]    
        
    for(int i=0;i<n;i++){
    cbartemp=cbars(i,_);  
    G3temp=G3(i,_);  
      logP(i,1)=logP(i,0)-NegLL(i)+0.5*arma::as_scalar(cbarrow.t() * cbarrow)+arma::as_scalar(G3row.t() * cbarrow);
      LLconst(i,0)=NegLL(i)-arma::as_scalar(G3row.t() * cbarrow);
    }
    
    
//    return logP;
    return Rcpp::List::create(Rcpp::Named("logP")=logP,Rcpp::Named("LLconst")=LLconst);
    
}


//////////////////////////////////////////////////////////////////////////////

Rcpp::List Set_Grid_C(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, 
Rcpp::NumericMatrix Lint,
Rcpp::NumericMatrix Down,
Rcpp::NumericMatrix Up,
Rcpp::NumericMatrix lglt,
Rcpp::NumericMatrix lgrt,
Rcpp::NumericMatrix lgct,
Rcpp::NumericMatrix logU,
Rcpp::NumericMatrix logP) {



  // Get dimensions
  
  int l1=GIndex.ncol();
  int l2=GIndex.nrow();
   

  // Main Function Code

  for(int j=0;j<l2;j++)
{
  logP(j,0)=0;
  
  for(int i=0;i<l1;i++)
{
  if(GIndex(j,i)==1){
  Down(j,i)=-INFINITY;
  Up(j,i)=Lint(0,i)+cbars(j,i);
  }
  if(GIndex(j,i)==2){
  Down(j,i)=Lint(0,i)+cbars(j,i);
  Up(j,i)=Lint(1,i)+cbars(j,i); 
  
  }
  if(GIndex(j,i)==3){
  Down(j,i)=Lint(1,i)+cbars(j,i); 
  Up(j,i)=INFINITY;
  }
  if(GIndex(j,i)==4){
    Down(j,i)=-INFINITY;
   Up(j,i)=INFINITY;
  }
}
}
  
  
  for(int i=0;i<l1;i++){
  lglt(_,i) = pnorm(Up(_,i),0.0,1.0,true,true);
  lgrt(_,i) = pnorm(Down(_,i),0.0,1.0,false,true);
  lgct(_,i) = pnorm(Up(_,i),0.0,1.0)-pnorm(Down(_,i),0.0,1.0); 
  }
  for(int j=0;j<l2;j++)
{
  for(int i=0;i<l1;i++)
{
  if(GIndex(j,i)==1){
  logU(j,i)=lglt(j,i);
  logP(j,0)=logP(j,0)+logU(j,i);
  }
  if(GIndex(j,i)==2){
  lgct(j,i)=log(lgct(j,i));
  logU(j,i)=lgct(j,i);
  logP(j,0)=logP(j,0)+logU(j,i);
  }
  if(GIndex(j,i)==3){
  logU(j,i)=lgrt(j,i);
  logP(j,0)=logP(j,0)+logU(j,i);
  }
  if(GIndex(j,i)==4){
  logU(j,i)=0;
  logP(j,0)=logP(j,0)+logU(j,i);
  }
}}
  
  // Return List

  return Rcpp::List::create(Rcpp::Named("Down")=Down
  ,Rcpp::Named("Up")=Up
  ,Rcpp::Named("lglt")=lglt
  ,Rcpp::Named("lgrt")=lgrt
  ,Rcpp::Named("lgct")=lgct
  ,Rcpp::Named("logU")=logU
  ,Rcpp::Named("logP")=logP);

}



//////////////////////////////////////////////////////////////////////////////

Rcpp::List   setlogP_C(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst) {

    int n = logP.nrow(), k = logP.ncol();
    int l1 =cbars.ncol();
    
      arma::mat logP2(logP.begin(), n, k, false); 
      NumericVector cbartemp=cbars(0,_);  
      NumericVector G3temp=G3(0,_);  
    
      arma::colvec cbarrow(cbartemp.begin(),l1,false);
      arma::colvec G3row(G3temp.begin(),l1,false);
    
        
      for(int i=0;i<n;i++){
        cbartemp=cbars(i,_);  
        G3temp=G3(i,_);  

        logP(i,1)=logP(i,0)-NegLL(i)+0.5*arma::as_scalar(cbarrow.t() * cbarrow)+arma::as_scalar(G3row.t() * cbarrow);

      LLconst(i,0)=NegLL(i)-arma::as_scalar(G3row.t() * cbarrow);
      }
    
    
    return Rcpp::List::create(Rcpp::Named("logP")=logP,Rcpp::Named("LLconst")=LLconst);
    
}



void Set_Grid_C2(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, 
Rcpp::NumericMatrix Lint,
Rcpp::NumericMatrix Down,
Rcpp::NumericMatrix Up,
Rcpp::NumericMatrix lglt,
Rcpp::NumericMatrix lgrt,
Rcpp::NumericMatrix lgct,
Rcpp::NumericMatrix logU,
Rcpp::NumericMatrix logP){
 
 
  // Get dimensions
  
  int l1=GIndex.ncol();
  int l2=GIndex.nrow();
  double g1;
  double g2;
    
  Rcpp::NumericMatrix lgct1(l2,l1);
  Rcpp::NumericMatrix lgct2(l2,l1);
  Rcpp::NumericMatrix lgct3(l2,l1);
  Rcpp::NumericMatrix lgct4(l2,l1);
  
  
  // This might be the problem. - cbars here appear to be all zeros..
  // Corrected cbars issue 06/03/20
  
//  Rcpp::Rcout << "Lint inside Set_Grid_C2 :" << std::endl << Lint << std::endl;
//  Rcpp::Rcout << "cbars inside Set_Grid_C2 :" << std::endl << cbars << std::endl;
   

  // Main Function Code

  for(int j=0;j<l2;j++)
{
  logP(j,0)=0;
  
  for(int i=0;i<l1;i++)
{
  if(GIndex(j,i)==1){
  Down(j,i)=-INFINITY;
  Up(j,i)=Lint(0,i)+cbars(j,i);
  }
  if(GIndex(j,i)==2){
  Down(j,i)=Lint(0,i)+cbars(j,i);
  Up(j,i)=Lint(1,i)+cbars(j,i); 
  
  }
  if(GIndex(j,i)==3){
  Down(j,i)=Lint(1,i)+cbars(j,i); 
  Up(j,i)=INFINITY;
  }
  if(GIndex(j,i)==4){
    Down(j,i)=-INFINITY;
   Up(j,i)=INFINITY;
  }
}
}
  
  // Should log_ctnorm be used here for improvement
  
  for(int i=0;i<l1;i++){
  lglt(_,i) = pnorm(Up(_,i),0.0,1.0,true,true);
  lgrt(_,i) = pnorm(Down(_,i),0.0,1.0,false,true);
  lgct(_,i) = pnorm(Up(_,i),0.0,1.0)-pnorm(Down(_,i),0.0,1.0); 
  
  // Use loop and if statements here
  // These calculations seem to yield same values as lgct above
  // when lgct is correct - should be safer version when differences 
  // between Up and Down are small
    
   lgct1(_,i)=pnorm(Up(_,i),0.0,1.0)*
     (1-exp(pnorm(Down(_,i),0.0,1.0,true,true)-pnorm(Up(_,i),0.0,1.0,true,true)));
   
   lgct2(_,i)=pnorm(Down(_,i),0.0,1.0,false)*
     (1-exp(pnorm(Up(_,i),0.0,1.0,false,true)-pnorm(Down(_,i),0.0,1.0,false,true)));
  
  // Replace old lgct calculation with results from lgct3 
  // Keep redundant part for now to allow for testing of impact/further changes
  
  for(int j=0;j<l2;j++) {
    g1=-Down(j,i);
    g2=Up(j,i);
    if(g1>=g2) lgct3(j,i)=lgct1(j,i);
    if(g2>g1) lgct3(j,i)=lgct2(j,i);
    lgct4(j,i)=lgct3(j,i)-lgct(j,i);
    lgct(j,i)=lgct3(j,i);
    
  }
  
//     -pnorm(Down(_,i),0.0,1.0); 
//     exp(pnorm(q=a,mean=mu,sd=sigma,log.p=TRUE)-pnorm(q=b2,mean=mu,sd=sigma,log.p=TRUE)))
  
  }

//  Rcpp::Rcout << "Up:" << std::endl << Up << std::endl;
//  Rcpp::Rcout << "Down:" << std::endl << Down << std::endl;
//  Rcpp::Rcout << "lglt:" << std::endl << lglt << std::endl;
//  Rcpp::Rcout << "lgrt:" << std::endl << lgrt << std::endl;
//  Rcpp::Rcout << "lgct:" << std::endl << lgct << std::endl;
//  Rcpp::Rcout << "lgct3:" << std::endl << lgct3 << std::endl;
//  Rcpp::Rcout << "lgct3-lgct:" << std::endl << lgct4 << std::endl;
  
    
  
// Loop through i and j

  for(int j=0;j<l2;j++)
{
  for(int i=0;i<l1;i++)
{
  if(GIndex(j,i)==1){
  logU(j,i)=lglt(j,i);
  logP(j,0)=logP(j,0)+logU(j,i);
  }
  if(GIndex(j,i)==2){
  lgct(j,i)=log(lgct(j,i));   // log of lgct here 
  logU(j,i)=lgct(j,i);
  logP(j,0)=logP(j,0)+logU(j,i);
  }
  if(GIndex(j,i)==3){
  logU(j,i)=lgrt(j,i);
  logP(j,0)=logP(j,0)+logU(j,i);
  }
  if(GIndex(j,i)==4){
  logU(j,i)=0;
  logP(j,0)=logP(j,0)+logU(j,i);
  }
}}

//  Rcpp::Rcout << "lct:" << std::endl << lgct << std::endl;
//  Rcpp::Rcout << "logU:" << std::endl << logU << std::endl;
  
   
 
}

void setlogP_C2(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst){
 
    int n = logP.nrow(), k = logP.ncol();
    int l1 =cbars.ncol();
    
      arma::mat logP2(logP.begin(), n, k, false); 
      NumericVector cbartemp=cbars(0,_);  
      NumericVector G3temp=G3(0,_);  
    
      arma::colvec cbarrow(cbartemp.begin(),l1,false);
      arma::colvec G3row(G3temp.begin(),l1,false);
    
        
      for(int i=0;i<n;i++){
        cbartemp=cbars(i,_);  
        G3temp=G3(i,_);  

        logP(i,1)=logP(i,0)-NegLL(i)+0.5*arma::as_scalar(cbarrow.t() * cbarrow)+arma::as_scalar(G3row.t() * cbarrow);

      LLconst(i,0)=NegLL(i)-arma::as_scalar(G3row.t() * cbarrow);
      }

 
}