// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"

using namespace Rcpp;

/// New Pointwise implementation (suggested by COPILOT)

void Set_Grid_C2_pointwise(
    NumericMatrix GIndex,
    NumericMatrix cbars,
    NumericMatrix Lint,
    NumericMatrix Down,
    NumericMatrix Up,
    NumericMatrix lglt,
    NumericMatrix lgrt,
    NumericMatrix lgct,
    NumericMatrix logU,
    NumericMatrix logP
) {
  int nGrid  = GIndex.nrow();
  int kParam = GIndex.ncol();
  
  // zero‐initialize logP
  for (int j = 0; j < nGrid; ++j) {
    logP(j, 0) = 0.0;
  }
  
  for (int j = 0; j < nGrid; ++j) {
    double logPj = 0.0;
    
    for (int i = 0; i < kParam; ++i) {
      int code = static_cast<int>(GIndex(j, i));
      double cb  = cbars(j, i);
      
      // 1. Determine [down, up]
      double down, up;
      if (code == 1) {
        down = R_NegInf;
        up   = Lint(0, i) + cb;
      }
      else if (code == 2) {
        down = Lint(0, i) + cb;
        up   = Lint(1, i) + cb;
      }
      else if (code == 3) {
        down = Lint(1, i) + cb;
        up   = R_PosInf;
      }
      else {  // code == 4
        down = R_NegInf;
        up   = R_PosInf;
      }
      Down(j, i) = down;
      Up(j,   i) = up;
      
      // 2. Compute CDFs and log‐CDFs
      // P(Z <= up), log P(Z <= up)
      double p_up_nl    = R::pnorm(up,   0.0, 1.0, /*lower_tail=*/ true,  /*log_p=*/ false);
      double logp_up    = R::pnorm(up,   0.0, 1.0, /*lower_tail=*/ true,  /*log_p=*/ true);
      
      // P(Z > down), log P(Z > down)
      double p_down_gt  = R::pnorm(down, 0.0, 1.0, /*lower_tail=*/ false, /*log_p=*/ false);
      double logp_down_gt = R::pnorm(down, 0.0, 1.0, /*lower_tail=*/ false, /*log_p=*/ true);
      
      // log‐tails for pure lower/upper
      lglt(j, i) = logp_up;         // log P(Z <= up)
      lgrt(j, i) = logp_down_gt;   // log P(Z >  down)
      
      // 3. Stable mid‐interval probability
      // p_mid = P(down < Z <= up)
      double p_mid;
      double g1 = -down;
      double g2 =  up;
      
      if (g1 >= g2) {
        // use P(Z <= up) * (1 - exp(log P(Z<=down) - log P(Z<=up)))
        double logp_down_le = R::pnorm(down, 0.0,1.0, /*lower_tail=*/ true, /*log_p=*/ true);
        p_mid = p_up_nl * (1.0 - exp(logp_down_le - logp_up));
      } else {
        // use P(Z > down) * (1 - exp(log P(Z>up)  - log P(Z>down)))
        double logp_up_gt = R::pnorm(up, 0.0,1.0, /*lower_tail=*/ false, /*log_p=*/ true);
        p_mid = p_down_gt * (1.0 - exp(logp_up_gt - logp_down_gt));
      }
      
      // log of the mid‐interval
      double logp_mid = (p_mid > 0.0) ? std::log(p_mid) : R_NaReal;
      lgct(j, i) = logp_mid;
      
      // 4. Choose the correct log‐density for this cell
      double uji = 0.0;
      if      (code == 1) uji = logp_up;
      else if (code == 2) uji = logp_mid;
      else if (code == 3) uji = logp_down_gt;
      // code == 4 ⇒ uji = 0
      
      logU(j, i)  = uji;
      logPj      += uji;
    }
    
    logP(j, 0) = logPj;
  }
}


// Set_Grid creates new objects and returns them is a list
// Set_Grid_C is likely an intermediate implementation that takes objects as inputs 
// and the returns them
// Set_Grid_C2 take existing allocated objects and populates them
// In cases where the function is called repeatedly the last approach is preferred

void Set_Grid_C2(Rcpp::NumericMatrix GIndex,  
                 Rcpp::NumericMatrix cbars, 
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

/////////////////////////////////////////////////////////////////////////////////////////////////
//    Set Grid Does not update existing objects - It creates new ones    ////////////////////////
//   This could be bad when function is called many times                 //////////////  


// [[Rcpp::export(".Set_Grid_cpp")]]


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

/////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/// This function does not return all components used to build the grid   /////  

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


