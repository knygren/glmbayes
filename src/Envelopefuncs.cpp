// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us

#include "RcppArmadillo.h"
#include "famfuncs.h"

using namespace Rcpp;


// [[Rcpp::export]]


List glmbenvelope_c(NumericVector bStar,NumericMatrix A,
                    NumericVector y, 
                    NumericMatrix x,
                    NumericMatrix mu,
                    NumericMatrix P,
                    NumericVector alpha,
                    NumericVector wt,
                    std::string family="binomial",
                    std::string link="logit",
                    int Gridtype=2, 
                    int n=1,
                    bool sortgrid=false
){
  
  int l1 = A.nrow(), k = A.ncol();
  arma::mat A2(A.begin(), l1, k, false);
  arma::colvec bStar_2(bStar.begin(), bStar.size(), false);
  
  
  NumericVector a_1(l1);
  arma::vec a_2(a_1.begin(), a_1.size(), false);
  
  NumericVector xx_1(3, 1.0);
  NumericVector xx_2=NumericVector::create(-1.0,0.0,1.0);
  NumericVector yy_1(2, 1.0);
  NumericVector yy_2=NumericVector::create(-0.5,0.5);
  NumericMatrix G1(3,l1);
  NumericMatrix Lint1(2,l1);
  arma::mat G1b(G1.begin(), 3, l1, false);
  arma::mat Lint(Lint1.begin(), 2, l1, false);
  
  arma::colvec xx_1b(xx_1.begin(), xx_1.size(), false);
  arma::colvec xx_2b(xx_2.begin(), xx_2.size(), false);
  arma::colvec yy_1b(yy_1.begin(), yy_1.size(), false);
  arma::colvec yy_2b(yy_2.begin(), yy_2.size(), false);
  List G2(a_1.size());
  List GIndex1(a_1.size());
  Rcpp::Function opGrid("optgrid");
  Rcpp::Function expGrid("expand.grid");
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function EnvSort("EnvelopeSort");
  
  
  int i;  
  
  //  bStar_2.print("bstar part of Grid calculation");
  
  a_2=arma::diagvec(A2);
  arma::vec omega=(sqrt(2)-arma::exp(-1.20491-0.7321*sqrt(0.5+a_2)))/arma::sqrt(1+a_2);
  G1b=xx_1b*arma::trans(bStar_2)+xx_2b*arma::trans(omega);
  Lint=yy_1b*arma::trans(bStar_2)+yy_2b*arma::trans(omega);
  
  // Second row in G1b here is the posterior mode
  
  //  G1b.print("G1b part of Grid calculation"); 
  //  Lint.print("Lint b part of Grid calculation");
  
  
  
  NumericVector gridindex(l1);
  
  if(Gridtype==2){
    gridindex=opGrid(a_2,n);
  }
  
  NumericVector Temp1=G1( _, 0);
  double Temp2;
  
  for(i=0;i<l1;i++){
    
    
    
    if(Gridtype==1){
      if((1+a_2[i])<=(2/sqrt(M_PI))){ 
        Temp2=G1(1,i);
        G2[i]=NumericVector::create(Temp2);
        GIndex1[i]=NumericVector::create(4.0);
      }
      if((1+a_2[i])>(2/sqrt(M_PI))){
        Temp1=G1(_,i);
        G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
        GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
      }    
    }  
    if(Gridtype==2){
      if(gridindex[i]==1){
        Temp2=G1(1,i);
        G2[i]=NumericVector::create(Temp2);
        GIndex1[i]=NumericVector::create(4.0);
      }
      if(gridindex[i]==3){
        Temp1=G1(_,i);
        G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
        GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
      }
    }
    
    if(Gridtype==3){
      Temp1=G1(_,i);
      G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
      GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
    }
    
    if(Gridtype==4){
      Temp2=G1(1,i);
      G2[i]=NumericVector::create(Temp2);
      GIndex1[i]=NumericVector::create(4.0);
    }
    
    
    
  }
  
  NumericMatrix G3=asMat(expGrid(G2));
  NumericMatrix GIndex=asMat(expGrid(GIndex1));
  NumericMatrix G4(G3.ncol(),G3.nrow());
  int l2=GIndex.nrow();
  
  arma::mat G3b(G3.begin(), G3.nrow(), G3.ncol(), false);
  arma::mat G4b(G4.begin(), G4.nrow(), G4.ncol(), false);
  
  //    G3b.print("expanded Grid - Should have info from G1b");
  
  G4b=trans(G3b);
  
  NumericMatrix cbars(l2,l1);
  NumericMatrix Up(l2,l1);
  NumericMatrix Down(l2,l1);
  NumericMatrix logP(l2,2);
  NumericMatrix logU(l2,l1);
  NumericMatrix loglt(l2,l1);
  NumericMatrix logrt(l2,l1);
  NumericMatrix logct(l2,l1);
  
  NumericMatrix LLconst(l2,1);
  NumericVector NegLL(l2);    
  arma::mat cbars2(cbars.begin(), l2, l1, false); 
  arma::mat cbars3(cbars.begin(), l2, l1, false); 
  
  // Note: NegLL_2 only added to allow for QC printing of results 
  
  arma::colvec NegLL_2(NegLL.begin(), NegLL.size(), false);
  
  //    G4b.print("tangent points");
  
  Rcpp::Rcout << "Gridtype is :"  << Gridtype << std::endl;
  Rcpp::Rcout << "Number of Variables in model are :"  << l1 << std::endl;
  Rcpp::Rcout << "Number of points in Grid are :"  << l2 << std::endl;
  
  //    Rcpp::Rcout << "mu passed to LL and Gradient Functions:" << std::flush << mu << std::endl;
  
  if( family=="binomial" && link=="logit"){
    //      Rcpp::Rcout << "Finding Values of Log-posteriors and Gradients - New function:" << std::endl;
    
    //    f4_binomial_logit(G4,y, x,mu,P,alpha,wt,NegLL,cbars,1);
    //    cbars2.print("Value of cbars2 after f4_binomial_logit");  
    
    //    NegLL_2.print("Value of NegLL_2 after f4_binomial_logit");  
    
    //    Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_logit(G4,y, x, mu, P, alpha, wt,1);  
    
    //    NegLL_2.print("Value of NegLL_2 after f2_binomial_logit");  
    
    
    //        Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    
    // This might point cbars2 to a different part of memory so cbars does not get updated
    cbars2=f3_binomial_logit(G4,y, x,mu,P,alpha,wt,1);
  }
  if(family=="binomial"  && link=="probit"){
    Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_probit(G4,y, x, mu, P, alpha, wt,1);  
    Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_binomial_probit(G4,y, x,mu,P,alpha,wt,1);
  }
  if(family=="binomial"   && link=="cloglog"){
    Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_cloglog(G4,y, x, mu, P, alpha, wt,1);  
    Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    
    cbars2=f3_binomial_cloglog(G4,y, x,mu,P,alpha,wt,1);
  }
  
  if(family=="quasibinomial"  && link=="logit"){
    Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_logit(G4,y, x, mu, P, alpha, wt,1);  
    Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_binomial_logit(G4,y, x,mu,P,alpha,wt,1);
  }
  if(family=="quasibinomial" && link=="probit"){
    Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_probit(G4,y, x, mu, P, alpha, wt,1);  
    Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_binomial_probit(G4,y, x,mu,P,alpha,wt,1);
  }
  if(family=="quasibinomial" && link=="cloglog"){
    Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_cloglog(G4,y, x, mu, P, alpha, wt,1);  
    Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_binomial_cloglog(G4,y, x,mu,P,alpha,wt,1);
  }
  
  if(family=="poisson" ){
    Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_poisson(G4,y, x, mu, P, alpha, wt,1);  
    Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_poisson(G4,y, x,mu,P,alpha,wt,1);
  }
  
  if(family=="quasipoisson" ){
    Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_poisson(G4,y, x, mu, P, alpha, wt,1);  
    Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_poisson(G4,y, x,mu,P,alpha,wt,1);
  }
  
  if(family=="Gamma" ){
    Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_gamma(G4,y, x, mu, P, alpha, wt,1);  
    Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_gamma(G4,y, x,mu,P,alpha,wt,1);
  }
  
  // 
  
  //    cbars2.print("cbars2 out of Gradient Valuations");
  //    Rcpp::Rcout << "Negative Log-likelihood at tangents:" << std::endl << NegLL << std::endl;
  
  Rcpp::Rcout << "Finished Log-posterior evaluations:" << std::endl;
  
  // Do a temporary correction here cbars3 should point to correct memory
  // See if this sets cbars
  
  
  cbars3=cbars2;
  
  // why aren't cbars being passed correctly?
  
  //    Rcpp::Rcout << "cbars being passed to Set_Grid_C2 :" << std::endl << cbars << std::endl;
  
  
  Set_Grid_C2(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);
  
  // Earlier version of this did not pass LLconst --> Did not carry it along
  // 
  
  setlogP_C2(logP,NegLL,cbars,G3,LLconst);
  
  NumericMatrix::Column logP2 = logP( _, 1);
  
  //    logP.print("logP matrix - before and after setlogP_C2");
  
  //    Rcpp::Rcout << "logP - before and after setlogP_C2 :" << std::endl << logP << std::endl;
  //    Rcpp::Rcout << "LLconst - after setlogP_C2 :" << std::endl << LLconst << std::endl;
  
  double  maxlogP=max(logP2);
  
  NumericVector PLSD=exp(logP2-maxlogP);
  
  double sumP=sum(PLSD);
  
  PLSD=PLSD/sumP;
  
  //    Rcpp::Rcout << "PLSD Vector - probabilities:" << std::endl << PLSD << std::endl;
  
  if(sortgrid==true){
    
    Rcpp::List outlist=EnvSort(l1,l2,GIndex,G3,cbars,logU,logrt,loglt,logP,LLconst,PLSD,a_1);
    
    return(outlist);
    
  }
  
  
  return Rcpp::List::create(Rcpp::Named("GridIndex")=GIndex,
                            Rcpp::Named("thetabars")=G3,
                            Rcpp::Named("cbars")=cbars,
                            Rcpp::Named("logU")=logU,
                            Rcpp::Named("logrt")=logrt,
                            Rcpp::Named("loglt")=loglt,
                            Rcpp::Named("LLconst")=LLconst,
                            Rcpp::Named("logP")=logP(_,0),
                            Rcpp::Named("PLSD")=PLSD,
                            Rcpp::Named("a1")=a_1
  );
  
  
}




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