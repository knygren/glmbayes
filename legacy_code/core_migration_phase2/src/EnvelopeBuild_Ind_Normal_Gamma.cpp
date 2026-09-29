// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// #include "famfuncs.h"
#include "Envelopefuncs.h"
#include <RcppParallel.h>
#include "openclPort.h"
#include "progress_utils.h"
#include "R_interface.h"

using namespace Rcpp;
using namespace openclPort;
// using namespace glmbayes::famfuncs;
using namespace glmbayes::progress;



namespace glmbayes {

namespace env {

List EnvelopeBuild_Ind_Normal_Gamma(NumericVector bStar,NumericMatrix A,
                                    NumericVector y, 
                                    NumericMatrix x,
                                    NumericMatrix mu,
                                    NumericMatrix P,
                                    NumericVector alpha,
                                    NumericVector wt,
                                    std::string family,
                                    std::string link,
                                    int Gridtype, 
                                    int n,
                                    int n_envopt,
                                    bool sortgrid,
                                    bool use_opencl    ,
                                    bool verbose       
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
  Rcpp::Function EnvelopeOpt = glmbayes_R::r_envelope_opt();
  Rcpp::Function expGrid("expand.grid");
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function EnvSort = glmbayes_R::r_envelope_sort();
  
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
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Starting EnvelopeEval (NegLL, cbars) at " << now_hms() << " <<<\n";
    }
    Timer t_eval1; if (verbose) t_eval1.begin();
    
    Rcpp::List eval_info = EnvelopeEval(G4, y, x, mu, P, alpha, wt, family, link, use_opencl, verbose);
    NegLL = eval_info["NegLL"];
    cbars2 = Rcpp::as<arma::mat>(eval_info["cbars"]);
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting EnvelopeEval (NegLL, cbars) at " << now_hms() << " <<<\n";
      print_completed("[EnvelopeBuild] EnvelopeEval (NegLL, cbars)", t_eval1);
    }
    
    
    //    Rcpp::List eval_info = EnvelopeEval(G4, y, x, mu, P, alpha, wt,
    //                                        family, link, use_opencl, verbose);
    
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Starting EnvelopeEval (slope variants) at " << now_hms() << " <<<\n";
    }
    Timer t_eval2; if (verbose) t_eval2.begin();
    
    Rcpp::List eval_info2 = EnvelopeEval(G4, y, x, mu, 0*P, alpha, wt, family, link, use_opencl, verbose);
    NegLL_slope  = eval_info2["NegLL"];
    cbars_slope2 = Rcpp::as<arma::mat>(eval_info2["cbars"]);
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting EnvelopeEval (slope variants) at " << now_hms() << " <<<\n";
      print_completed("[EnvelopeBuild] EnvelopeEval (slope variants)", t_eval2);
      Rcpp::Rcout << "[EnvelopeBuild] Finished assigning NegLL_slope and cbars_slope2\n";
    }
    
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Starting RSS evaluation at " << now_hms() << " <<<\n";
    }
    Timer t_rss; if (verbose) t_rss.begin();
    
    RSS_Out = RSS(y, x, G4, alpha, wt); // includes dispersion in weight
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting RSS evaluation at " << now_hms() << " <<<\n";
      print_completed("[EnvelopeBuild] RSS evaluation", t_rss);
    } 
  }
  
  
  //  Rcpp::Rcout << "Finished Log-posterior evaluations:" << std::endl;
  
  // Do a temporary correction here cbars3 should point to correct memory
  // See if this sets cbars
  
  cbars3=cbars2;
  cbars_slope3=cbars_slope2;
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Entering EnvelopeSet_Grid_C2 at " << now_hms() << " <<<\n";
  }
  Timer t_setgrid; if (verbose) t_setgrid.begin();
  
  EnvelopeSet_Grid_C2(GIndex, cbars, Lint1, Down, Up, loglt, logrt, logct, logU, logP);
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting EnvelopeSet_Grid_C2 at " << now_hms() << " <<<\n";
    print_completed("[EnvelopeBuild] EnvelopeSet_Grid_C2", t_setgrid);
  }
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Entering Set_logP_C2 at " << now_hms() << " <<<\n";
  }
  Timer t_setlogp; if (verbose) t_setlogp.begin();
  
  setlogP_C2(logP, NegLL, cbars, G3, LLconst);
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting Set_logP_C2 at " << now_hms() << " <<<\n";
    print_completed("[EnvelopeBuild] Set_logP_C2", t_setlogp);
  }  
  
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Starting PLSD computation at " << now_hms() << " <<<\n";
  }
  Timer t_plsd; if (verbose) t_plsd.begin();
  
  NumericMatrix::Column logP2 = logP(_, 1);
  double maxlogP = max(logP2);
  NumericVector PLSD = exp(logP2 - maxlogP);
  double sumP = sum(PLSD);
  PLSD = PLSD / sumP;
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting PLSD computation at " << now_hms() << " <<<\n";
    print_completed("[EnvelopeBuild] PLSD computation", t_plsd);
  }
  
  
  
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


} // env
} // glmbayes
