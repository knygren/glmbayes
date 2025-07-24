// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Set_Grid.h"

using namespace Rcpp;

// [[Rcpp::export(".EnvelopeBuild_cpp")]]

List EnvelopeBuild_c(NumericVector bStar,NumericMatrix A,
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
  

  Rcpp::Rcout << "Entering EnvelopeBuild_c: "
              << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
              << "\n";
  
  int progbar=0;
  
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
  Rcpp::Function EnvelopeOpt("EnvelopeOpt");
  Rcpp::Function expGrid("expand.grid");
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function EnvSort("EnvelopeSort");

  int i;  
  
  a_2=arma::diagvec(A2);
  arma::vec omega=(sqrt(2)-arma::exp(-1.20491-0.7321*sqrt(0.5+a_2)))/arma::sqrt(1+a_2);
  G1b=xx_1b*arma::trans(bStar_2)+xx_2b*arma::trans(omega);
  Lint=yy_1b*arma::trans(bStar_2)+yy_2b*arma::trans(omega);
  
  // Second row in G1b here is the posterior mode
  
  NumericVector gridindex(l1);
  
  if(Gridtype==2){
    gridindex=EnvelopeOpt(a_2,n);
  }
  
  NumericVector Temp1=G1( _, 0);
  double Temp2;

  
  Rcpp::Rcout << "Entering Envelope Loop: "
              << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
              << "\n";
  
      
  // Should write a small note with logic behind types 1 and 2
  
  for(i=0;i<l1;i++){
  
    if(Gridtype==1){
      
      // For Gridtype==1, small 1+a[i]<=(2/sqrt(M_PI) yields grid over full line
      // Can check speed for simulation when Gridtype=1 vs. Gridtyp=2 or 3     
      
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

  
//  Rcpp::Rcout << "Exiting Envelope Loop: "
//              << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
//              << "\n";
  
    
  NumericMatrix G3=asMat(expGrid(G2));
  NumericMatrix GIndex=asMat(expGrid(GIndex1));
  NumericMatrix G4(G3.ncol(),G3.nrow());
  int l2=GIndex.nrow();
  
  arma::mat G3b(G3.begin(), G3.nrow(), G3.ncol(), false);
  arma::mat G4b(G4.begin(), G4.nrow(), G4.ncol(), false);
  
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
  
//  Rcpp::Rcout << "Gridtype is :"  << Gridtype << std::endl;
//  Rcpp::Rcout << "Number of Variables in model are :"  << l1 << std::endl;
//  Rcpp::Rcout << "Number of points in Grid are :"  << l2 << std::endl;
  
  if( family=="binomial" && link=="logit"){
    
    Rcpp::Rcout << "Initiating NegLL Calculations: "
                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                << "\n";
    
        
    //Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_logit(G4,y, x, mu, P, alpha, wt,progbar);  
    //Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    //Rcout << "Initiating Gradient Calculations: " << Rcpp::as<std::string>(Rcpp::Function("Sys.time")()) << "\n";

    Rcpp::Rcout << "Initiating Gradient Calculations: "
                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                << "\n";
    
        
        cbars2=f3_binomial_logit(G4,y, x,mu,P,alpha,wt,progbar);
  }
  if(family=="binomial"  && link=="probit"){
  //  Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_probit(G4,y, x, mu, P, alpha, wt,progbar);  
  //  Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_binomial_probit(G4,y, x,mu,P,alpha,wt,progbar);
  }
  if(family=="binomial"   && link=="cloglog"){
  //  Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_cloglog(G4,y, x, mu, P, alpha, wt,progbar);  
  //  Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_binomial_cloglog(G4,y, x,mu,P,alpha,wt,progbar);
  }
  
  if(family=="quasibinomial"  && link=="logit"){
    //Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_logit(G4,y, x, mu, P, alpha, wt,progbar);  
    //Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_binomial_logit(G4,y, x,mu,P,alpha,wt,progbar);
  }
  if(family=="quasibinomial" && link=="probit"){
    //Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_probit(G4,y, x, mu, P, alpha, wt,progbar);  
    //Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_binomial_probit(G4,y, x,mu,P,alpha,wt,progbar);
  }
  if(family=="quasibinomial" && link=="cloglog"){
    //Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_binomial_cloglog(G4,y, x, mu, P, alpha, wt,progbar);  
    //Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_binomial_cloglog(G4,y, x,mu,P,alpha,wt,progbar);
  }
  
  if(family=="poisson" ){
    //Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_poisson(G4,y, x, mu, P, alpha, wt,progbar);  
    //Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_poisson(G4,y, x,mu,P,alpha,wt,progbar);
  }
  
  if(family=="quasipoisson" ){
    //Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_poisson(G4,y, x, mu, P, alpha, wt,progbar);  
    //Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_poisson(G4,y, x,mu,P,alpha,wt,progbar);
  }
  
  if(family=="Gamma" ){
    //Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_gamma(G4,y, x, mu, P, alpha, wt,progbar);  
    //Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_gamma(G4,y, x,mu,P,alpha,wt,progbar);
  }
  
  if(family=="gaussian" ){
    //Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    NegLL=f2_gaussian(G4,y, x, mu, P, alpha, wt);  
    //Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_gaussian(G4,y, x,mu,P,alpha,wt);
  }
  

//  Rcpp::Rcout << "Finished cbars Calculations: "
//              << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
//              << "\n";
  
    
  //  Rcpp::Rcout << "Finished Log-posterior evaluations:" << std::endl;
  
  // Do a temporary correction here cbars3 should point to correct memory
  // See if this sets cbars

  cbars3=cbars2;
  
  // July 2025 - Parallelization Implementation in steps
  
  // 1) Set_Grid_C2_pointwise changes loop to enable parallel processing (suggested by Copilot)


//  Rcpp::Rcout << "Entering Set grid C2 pointwise: "
//              << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
//              << "\n";
  
  
//  Set_Grid_C2(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);
  Set_Grid_C2_pointwise(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);
  
    

//    Rcpp::Rcout << "Entering setlogP_C2: "
//                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
//                << "\n";
    
        
      setlogP_C2(logP,NegLL,cbars,G3,LLconst);

      
//      Rcpp::Rcout << "Exiting setlogP_C2: "
//                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
//                  << "\n";
      

  
  NumericMatrix::Column logP2 = logP( _, 1);
  

  double  maxlogP=max(logP2);
  
  NumericVector PLSD=exp(logP2-maxlogP);
  
  double sumP=sum(PLSD);
  
  PLSD=PLSD/sumP;

//  Rcout << "Entering Enveloped sort: " << Rcpp::as<std::string>(Rcpp::Function("Sys.time")()) << "\n";
  
    
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





// [[Rcpp::export(".EnvelopeBuild_Ind_Normal_Gamma")]]

List EnvelopeBuild_Ind_Normal_Gamma(NumericVector bStar,NumericMatrix A,
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
  
  
//  int progbar=0;
  
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
  Rcpp::Function EnvelopeOpt("EnvelopeOpt");
  Rcpp::Function expGrid("expand.grid");
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function EnvSort("EnvelopeSort");
  
  int i;  
  
  a_2=arma::diagvec(A2);
  arma::vec omega=(sqrt(2)-arma::exp(-1.20491-0.7321*sqrt(0.5+a_2)))/arma::sqrt(1+a_2);
  G1b=xx_1b*arma::trans(bStar_2)+xx_2b*arma::trans(omega);
  Lint=yy_1b*arma::trans(bStar_2)+yy_2b*arma::trans(omega);
  
  // Second row in G1b here is the posterior mode
  
  NumericVector gridindex(l1);
  
  if(Gridtype==2){
    gridindex=EnvelopeOpt(a_2,n);
  }
  
  NumericVector Temp1=G1( _, 0);
  double Temp2;
  
  // Should write a small note with logic behind types 1 and 2
  
  for(i=0;i<l1;i++){
    
    if(Gridtype==1){
      
      // For Gridtype==1, small 1+a[i]<=(2/sqrt(M_PI) yields grid over full line
      // Can check speed for simulation when Gridtype=1 vs. Gridtyp=2 or 3     
      
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
  
  G4b=trans(G3b);
  
  NumericMatrix cbars(l2,l1);
  NumericMatrix cbars_slope(l2,l1);
  NumericMatrix Up(l2,l1);
  NumericMatrix Down(l2,l1);
  NumericMatrix logP(l2,2);
  NumericMatrix logU(l2,l1);
  NumericMatrix loglt(l2,l1);
  NumericMatrix logrt(l2,l1);
  NumericMatrix logct(l2,l1);
  
  NumericMatrix LLconst(l2,1);
  NumericVector NegLL(l2);    
  NumericVector NegLL_slope(l2);    
  NumericVector RSS_Out(l2);
  arma::mat cbars2(cbars.begin(), l2, l1, false); 
  arma::mat cbars3(cbars.begin(), l2, l1, false); 
  
  arma::mat cbars_slope2(cbars_slope.begin(), l2, l1, false); 
  arma::mat cbars_slope3(cbars_slope.begin(), l2, l1, false); 
  
  
    // Note: NegLL_2 only added to allow for QC printing of results 
  
  arma::colvec NegLL_2(NegLL.begin(), NegLL.size(), false);
  
  //    G4b.print("tangent points");
  
//  Rcpp::Rcout << "Gridtype is :"  << Gridtype << std::endl;
//  Rcpp::Rcout << "Number of Variables in model are :"  << l1 << std::endl;
//  Rcpp::Rcout << "Number of points in Grid are :"  << l2 << std::endl;
  

  if(family=="gaussian" ){
    //Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;

    // Adjust the slope calculations to split into several terms:
    // (i) Terms from shifted "prior" that does not depend on the dispersion
    // (ii) Constant terms from the actual LL that do not depend on dispersion or beta
    // (iii) Term from the LL that depends on the dispersion but not beta
    // (iv) Term from the LL that depends on beta and the dispersion (scaled RSS)
    
    NegLL=f2_gaussian(G4,y, x, mu, P, alpha, wt);  
    NegLL_slope=f2_gaussian(G4,y, x, mu, 0*P, alpha, wt);  
    //Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_gaussian(G4,y, x,mu,P,alpha,wt);
    cbars_slope2=f3_gaussian(G4,y, x,mu,0*P,alpha,wt);
    RSS_Out=RSS(y, x,G4,alpha,wt); // Note currenly includes the dispersion in the weight
    
  }
  
  
  //  Rcpp::Rcout << "Finished Log-posterior evaluations:" << std::endl;
  
  // Do a temporary correction here cbars3 should point to correct memory
  // See if this sets cbars
  
  cbars3=cbars2;
  cbars_slope3=cbars_slope2;
  
  Set_Grid_C2(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);
  
  setlogP_C2(logP,NegLL,cbars,G3,LLconst);
  
  NumericMatrix::Column logP2 = logP( _, 1);
  
  
  double  maxlogP=max(logP2);
  
  NumericVector PLSD=exp(logP2-maxlogP);
  
  double sumP=sum(PLSD);
  
  PLSD=PLSD/sumP;
  
  // Add sorting step back later after modifying EnvSort function
  // Should accomodate ready List
  
//  if(sortgrid==true){
//    Rcpp::List outlist=EnvSort(l1,l2,GIndex,G3,cbars,logU,logrt,loglt,logP,LLconst,PLSD,a_1);
//    return(outlist);
//  }
  
  
  return Rcpp::List::create(Rcpp::Named("GridIndex")=GIndex,
                            Rcpp::Named("thetabars")=G3,
                            Rcpp::Named("cbars")=cbars,
                            Rcpp::Named("cbars_slope")=cbars_slope,
                            Rcpp::Named("NegLL")=NegLL,
                            Rcpp::Named("NegLL_slope")=NegLL_slope,
                            Rcpp::Named("Lint1")=Lint1,
                            Rcpp::Named("RSS_Out")=RSS_Out,
                            Rcpp::Named("logU")=logU,
                            Rcpp::Named("logrt")=logrt,
                            Rcpp::Named("loglt")=loglt,
                            Rcpp::Named("LLconst")=LLconst,
                            Rcpp::Named("logP")=logP(_,0),
                            Rcpp::Named("PLSD")=PLSD,
                            Rcpp::Named("a1")=a_1
  );
  
  
}


