// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;

// Function declarations

NumericVector dbinom_glmb( NumericVector x, NumericVector N, NumericVector means, int lg);
NumericVector  f2_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
arma::mat  f3_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
arma::mat  f3_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
arma::mat  f3_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);


NumericVector dpois_glmb( NumericVector x, NumericVector means, int lg);
NumericVector  f2_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
arma::mat  f3_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);

NumericVector dgamma_glmb( NumericVector x, NumericVector shape, NumericVector scale, int lg);
NumericVector  f2_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
arma::mat  f3_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);

void Set_Grid_C2(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, 
Rcpp::NumericMatrix Lint,
Rcpp::NumericMatrix Down,
Rcpp::NumericMatrix Up,
Rcpp::NumericMatrix lglt,
Rcpp::NumericMatrix lgrt,
Rcpp::NumericMatrix lgct,
Rcpp::NumericMatrix logU,
Rcpp::NumericMatrix logP);

void setlogP_C2(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst);




double ctrnorm_cpp(double lgrt,double lglt,double mu,double sigma){
  RNGScope scope;
    
    double U=0;
    double out=0;
    double lgU2=0;
  
  if(lgrt>=lglt){
    U=R::runif(0.0, 1.0);
    double  u1=1-exp(lgrt);
		double lgu1=log(u1);
    lgU2=log(U)+lglt+log(1-exp(lgu1-lglt));
  	double lgU3=lgU2+log(1+exp(lgu1-lgU2));
  	out=R::qnorm(lgU3,mu,sigma,TRUE,TRUE);

  
  }  
  
  if(lgrt<lglt){
    U=R::runif(0.0, 1.0);
    double e1mu2=1-exp(lglt);
		double lg1mu2=log(e1mu2);
		double lgU2=log(U)+lgrt+log(1-exp(lg1mu2-lgrt));
		double lgU3=lgU2+log(1+exp(lg1mu2-lgU2));
     out=R::qnorm(lgU3,mu,sigma,FALSE,TRUE);

  }

 
  return out;

}

// [[Rcpp::export]]
Rcpp::List  glmbsim_cpp(int n,NumericVector y,NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,Function f2,Rcpp::List  Envelope,Rcpp::CharacterVector   family,Rcpp::CharacterVector   link)
{
  RNGScope scope;
  int l1 = mu.nrow();
//  int l2=pow(3,l1);
  
  std::string family2 = Rcpp::as<std::string>(family);
  std::string link2 = Rcpp::as<std::string>(link);  
  
  
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(n,l1);
  double a1=0;
  double a2=0;
  double U=0;
  double test=0;
  double U2=0;
  
  //out(0,0)=1;

  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
  NumericMatrix cbars=Envelope["cbars"];
  NumericVector LLconst=Envelope["LLconst"]; 

  NumericVector outtemp=out(0,_);
  arma::rowvec outtemp2(outtemp.begin(),l1,false);
  NumericVector cbartemp=cbars(0,_);
  arma::rowvec cbartemp2(cbartemp.begin(),l1,false);
  NumericMatrix testtemp(1,1);
  arma::mat testtemp2(testtemp.begin(),1,1,false);
  NumericMatrix btemp(l1,1);
  arma::mat btemp2(btemp.begin(),l1,1,false); 
  NumericVector testll(1);

    for(int i=0;i<n;i++){
      
    Rcpp::checkUserInterrupt();

    a1=0;
    draws(i)=1;
    while(a1==0){
      
      U=R::runif(0.0, 1.0);
      a2=0;
      J(i)=0;    
      while(a2==0){
        if(U<=PLSD(J(i))) a2=1;
          		if(U>PLSD(J(i))){ 
                U=U-PLSD(J(i));
							J(i)=J(i)+1;
             
							}
        //a2=1; 
          }
       for(int j=0;j<l1;j++){  
          
          out(i,j)=ctrnorm_cpp(logrt(J(i),j),loglt(J(i),j),-cbars(J(i),j),1.0);    

          
        }
     outtemp=out(i,_);
     cbartemp=cbars(J(i),_);
     testtemp2=outtemp2 * trans(cbartemp2);
     U2=R::runif(0.0, 1.0);
     btemp2=trans(outtemp2);    

    // Need to modify to call correct f2 function based on family and link

    
      if(family2=="binomial"){
      if(link2=="logit"){  
          testll=f2_binomial_logit(btemp,y, x,mu,P,alpha,wt);
            }
      if(link2=="probit"){  
          testll=f2_binomial_probit(btemp,y, x,mu,P,alpha,wt);
            }
      if(link2=="cloglog"){  
          testll=f2_binomial_cloglog(btemp,y, x,mu,P,alpha,wt);
            }
      }
      
      if(family2=="quasibinomial"){
      if(link2=="logit"){  
          testll=f2_binomial_logit(btemp,y, x,mu,P,alpha,wt);
            }
      if(link2=="probit"){  
          testll=f2_binomial_probit(btemp,y, x,mu,P,alpha,wt);
            }
      if(link2=="cloglog"){  
          testll=f2_binomial_cloglog(btemp,y, x,mu,P,alpha,wt);
            }
      }


      if(family2=="poisson"){  
      testll=f2_poisson(btemp,y, x,mu,P,alpha,wt);
      }
       if(family2=="quasipoisson"){  
      testll=f2_poisson(btemp,y, x,mu,P,alpha,wt);
      }

      if(family2=="Gamma"){  
      testll=f2_gamma(btemp,y, x,mu,P,alpha,wt);
      }


      test=LLconst(J(i))+testtemp(0,0)-log(U2)-testll(0);

      if(test>=0) a1=1;
	    if(test<0) draws(i)=draws(i)+1;

      
      }
    
    
    }

//return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws,Rcpp::Named("J")=J,Rcpp::Named("PLSD")=PLSD,Rcpp::Named("famout")=family);
return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("draws")=draws);

}

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
  
  a_2=arma::diagvec(A2);
  arma::vec omega=(sqrt(2)-arma::exp(-1.20491-0.7321*sqrt(0.5+a_2)))/arma::sqrt(1+a_2);
  G1b=xx_1b*arma::trans(bStar_2)+xx_2b*arma::trans(omega);
  Lint=yy_1b*arma::trans(bStar_2)+yy_2b*arma::trans(omega);
  
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



    if( family=="binomial" && link=="logit"){
    NegLL=f2_binomial_logit(G4,y, x, mu, P, alpha, wt);  
    cbars2=f3_binomial_logit(G4,y, x,mu,P,alpha,wt);
    }
    if(family=="binomial"  && link=="probit"){
    NegLL=f2_binomial_probit(G4,y, x, mu, P, alpha, wt);  
    cbars2=f3_binomial_probit(G4,y, x,mu,P,alpha,wt);
    }
    if(family=="binomial"   && link=="cloglog"){
    NegLL=f2_binomial_cloglog(G4,y, x, mu, P, alpha, wt);  
    cbars2=f3_binomial_cloglog(G4,y, x,mu,P,alpha,wt);
    }

    if(family=="quasibinomial"  && link=="logit"){
    NegLL=f2_binomial_logit(G4,y, x, mu, P, alpha, wt);  
    cbars2=f3_binomial_logit(G4,y, x,mu,P,alpha,wt);
    }
    if(family=="quasibinomial" && link=="probit"){
    NegLL=f2_binomial_probit(G4,y, x, mu, P, alpha, wt);  
    cbars2=f3_binomial_probit(G4,y, x,mu,P,alpha,wt);
    }
    if(family=="quasibinomial" && link=="cloglog"){
    NegLL=f2_binomial_cloglog(G4,y, x, mu, P, alpha, wt);  
    cbars2=f3_binomial_cloglog(G4,y, x,mu,P,alpha,wt);
    }

    if(family=="poisson" ){
    NegLL=f2_poisson(G4,y, x, mu, P, alpha, wt);  
    cbars2=f3_poisson(G4,y, x,mu,P,alpha,wt);
    }
    
    if(family=="quasipoisson" ){
    NegLL=f2_poisson(G4,y, x, mu, P, alpha, wt);  
    cbars2=f3_poisson(G4,y, x,mu,P,alpha,wt);
    }

    if(family=="Gamma" ){
    NegLL=f2_gamma(G4,y, x, mu, P, alpha, wt);  
    cbars2=f3_gamma(G4,y, x,mu,P,alpha,wt);
    }


    Set_Grid_C2(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);

    setlogP_C2(logP,NegLL,cbars,G3,LLconst);


    NumericMatrix::Column logP2 = logP( _, 1);

    double  maxlogP=max(logP2);
  
    NumericVector PLSD=exp(logP2-maxlogP);

    double sumP=sum(PLSD);

    PLSD=PLSD/sumP;

    
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


        p1=pnorm(xbtemp,0.0,1.0);
        p2=pnorm(-xbtemp,0.0,1.0);
        d1=dnorm(xbtemp,0.0,1.0);

        xbtemp2=alpha2+ x2 * b2;

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

// [[Rcpp::export]]

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





// [[Rcpp::export]]

Rcpp::List glmbsim_NGauss_cpp(int n,NumericVector y,NumericMatrix x, 
NumericVector mu,NumericMatrix P,NumericVector offset2,NumericVector wt,double dispersion,Rcpp::List
famfunc, Function f1,Function f2,Function f3,NumericVector start,
      std::string family="binomial",
      std::string link="logit",
      int Gridtype=2      
) {

    Rcpp::Function asMat("as.matrix");
    Rcpp::Function asVec("as.vector");
    int l1=x.ncol();
    int l2=x.nrow();

    double dispersion2;
    NumericVector alpha(l2);
    NumericMatrix mu2a=asMat(mu);

    arma::mat x2(x.begin(), l2, l1, false);
    arma::vec alpha2(alpha.begin(),l2,false);  
    arma::vec offset2b(offset2.begin(),l2,false);  
    arma::mat mu2(mu2a.begin(), mu2a.nrow(), mu2a.ncol(), false);

    NumericMatrix x2b(clone(x));
    arma::mat P2(P.begin(), P.nrow(), P.ncol(), false);

    if(family=="poisson"||family=="binomial")dispersion2=1;
	  else dispersion2=dispersion;
  
    NumericVector  wt2=wt/dispersion2;
    arma::vec wt3(wt2.begin(), x.nrow());
    
    alpha2=x2*mu2+offset2b;

    int i;

    
    NumericVector parin=start-mu;
    NumericVector mu1=mu-mu;
    Rcpp::Function optfun("optim");
  
//  Rcpp::List funclist=Rcpp::List::create(Rcpp::Named("f2")=f2,
//  Rcpp::Named("f3")=f3);  
    
//  return(funclist);

    List opt=optfun(_["par"]=parin,_["fn"]=f2, _["gr"]=f3,_["y"]=y,
    _["x"]=x,
    _["mu"]=mu1,_["P"]=P,_["alpha"]=alpha,_["wt"]=wt2,_["method"]="BFGS",_["hessian"]=true);

 

    NumericMatrix b2a=asMat(opt[0]);
    arma::mat b2(b2a.begin(), b2a.nrow(), b2a.ncol(), false);

//    b2.print("b2 inside rglmb");

    NumericVector min1=opt[1];
    int conver1=opt[3];
    NumericMatrix A1=opt[5];

    
    if(conver1>0){Rcpp::stop("Posterior Optimization failed");}
    
    arma::mat A1_b(A1.begin(), l1, l1, false); 
    arma::vec mu_0(mu.begin(), l1, false);

//    A1_b.print("A1_b inside rglmb");

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
//		eStemp <- eigen(P4_temp, symmetric = TRUE)
//    		evtemp <- min(eStemp$values)
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

    if(n==1){
    Envelope=glmbenvelope_c(b5, A4_1,y, x4_1,mu4_1,
    P5_1,alpha,wt2,family,link,Gridtype, n,false);
    }
    if(n>1){
    Envelope=glmbenvelope_c(b5, A4_1,y, x4_1,mu4_1,P5_1,alpha,wt2,family,link,Gridtype, n,true);
    }
    

    Rcpp::List sim=glmbsim_cpp(n,y,x4_1,mu4_1,P5_1,alpha,wt2,f2,Envelope,family,link);

    NumericMatrix sim2=sim[0];
    arma::mat sim2b(sim2.begin(), sim2.nrow(), sim2.ncol(), false);
    NumericMatrix out(l1,n);
    arma::mat out2(out.begin(), out.nrow(), out.ncol(), false);
    
    out2=L2Inv*L3Inv*trans(sim2b);
    NumericVector LL(n);
    
    for(i=0;i<n;i++){
    out(_,i)=out(_,i)+mu;
      LL[i]=as<double>(f1(_["b"]=out(_,i),_["y"]=y,_["x"]=x,offset2,wt2));
    }
    

  Rcpp::List Prior=Rcpp::List::create(Rcpp::Named("mean")=mu,Rcpp::Named("Precision")=P);  
    
  Rcpp::List outlist=Rcpp::List::create(Rcpp::Named("coefficients")=trans(out2),
          Rcpp::Named("PostMode")=b2a+mu,
          Rcpp::Named("Prior")=Prior,
          Rcpp::Named("iters")=sim[1],
          Rcpp::Named("famfunc")=famfunc,
          Rcpp::Named("Envelope")=Envelope,
          Rcpp::Named("dispersion")=dispersion2,
          Rcpp::Named("loglike")=LL
          );  



    return(outlist);


}


// [[Rcpp::export]]

Rcpp::List glmbsim_Gauss_cpp(int n,NumericVector y,NumericMatrix x, 
NumericVector mu,NumericMatrix P,NumericVector offset2,NumericVector wt,double dispersion,Rcpp::List
famfunc, Function f1,Function f2,Function f3,NumericVector start,
      std::string family="binomial",
      std::string link="logit",
      int Gridtype=2      
) {

    // Need to check combination of weighting and offset working properly

    Rcpp::Function asMat("as.matrix");
    Rcpp::Function asVec("as.vector");
    int l1=x.ncol();
    int l2=x.nrow();

    double dispersion2=dispersion;
//    NumericVector alpha(l2);
    NumericMatrix mu2a=asMat(mu);

    arma::mat x2(x.begin(), l2, l1, false);
//    arma::vec alpha2(alpha.begin(),l2,false);  
    arma::vec offset2b(offset2.begin(),l2,false);  
    arma::mat mu2(mu2a.begin(), mu2a.nrow(), mu2a.ncol(), false);

    NumericMatrix x2b(clone(x));
    arma::mat x2bb(x2b.begin(), l2, l1, false);
    arma::mat P2(P.begin(), P.nrow(), P.ncol(), false);

  
    NumericVector  wt2=wt/dispersion2;
    arma::vec wt3(wt2.begin(), x.nrow());
    
    // Is this needed (Check Offset code works properly-Both new and old code)
    

    // Move inside If Statement once confirmed works
    
        int i;
        
        // Should this subtract alpha2 [subtraction of offset likely makes this ok]?
        
        NumericVector  y1=y-offset2;
        arma::vec y2b(y1.begin(),l2,false);
        NumericMatrix W1(l2+l1,l1);
        arma::mat W(W1.begin(), W1.nrow(), W1.ncol(), false);
        NumericVector z1(l2+l1);
        arma::vec z(z1.begin(),l2+l1,false);

    
        for(i=0;i<l2;i++){
          x2b(i,_) =x2b(i,_)*sqrt(wt2[i]);
          y1(i)=y1(i)*sqrt(wt2[i]);
          }        
          
        arma::mat RA=arma::chol(P2);

        W.rows(0,l2-1)=x2bb;
        W.rows(l2,l2+l1-1)=RA;
        z.rows(0,l2-1)=y2b;
        z.rows(l2,l1+l2-1)=RA*mu2;
        arma::mat IR=arma::inv(trimatu(chol(trans(W)*W)));
        arma::mat b2=(IR*trans(IR))*(trans(W)*z);
    
        NumericMatrix out(n,l1);
        arma::mat out2(out.begin(), out.nrow(), out.ncol(), false);
        NumericVector draws(n,1.0);
        NumericVector LL(n);
        
        NumericMatrix U1(l1,n);
        arma::mat U2(U1.begin(), U1.nrow(), U1.ncol(), false);     
        
        for(i=0;i<n;i++){
           U1( _, i)=rnorm(l1);
          out2.row(i)=trans(b2+IR*U2.col(i));
//  Log-likelihood not set correctly if wrong famfuncs passed to this function
// Poisson throws warning....

//        LL[i]=as<double>(f1(_["b"]=out(i,_),_["y"]=y,_["x"]=x,offset2,wt2));

          }
      Rcpp::List Prior=Rcpp::List::create(Rcpp::Named("mean")=mu,Rcpp::Named("Precision")=P);  

      Rcpp::List outlist=Rcpp::List::create(Rcpp::Named("coefficients")=out,
              Rcpp::Named("PostMode")=b2,
              Rcpp::Named("Prior")=Prior,
              Rcpp::Named("iters")=draws,
              Rcpp::Named("famfunc")=famfunc,
              Rcpp::Named("dispersion")=dispersion2,
//              Rcpp::Named("U")=U1,
              Rcpp::Named("loglike")=LL
          );  

    return(outlist);

}



// [[Rcpp::export]]

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

    
    
    if(n==1){
    Envelope=glmbenvelope_c(b5, A4_1,y, x4_1,mu4_1,
    P5_1,alpha,wt2,family,link,Gridtype, n,false);
    }
    if(n>1){
    Envelope=glmbenvelope_c(b5, A4_1,y, x4_1,mu4_1,P5_1,alpha,wt2,family,link,Gridtype, n,true);
    }


    Rcpp::List sim=glmbsim_cpp(n,y,x4_1,mu4_1,P5_1,alpha,wt2,f2,Envelope,family,link);


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
                  double mu_const){
  
  
  double golden_ratio = 2/(sqrt(5) + 1);
  
  // Use the golden ratio to set the initial test points
  
  double  n1 = upper_bound - golden_ratio*(upper_bound - lower_bound);
  double  n2 = lower_bound + golden_ratio*(upper_bound - lower_bound);
  
  
  // ### Evaluate the function at the test points
  
  double epsilon_temp2=get_epsilon1( rstar2,epsilon, U_out,alpha_out,n1,  gammastar,t_star, mu_const);
  double  f1 = (epsilon_temp2-0.01)*(epsilon_temp2-0.01);
  
  
  
  epsilon_temp2=get_epsilon1( rstar2,epsilon, U_out,alpha_out,n2,  gammastar,t_star, mu_const);
  double  f2 = (epsilon_temp2-0.01)*(epsilon_temp2-0.01);
  
  
  
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
      f1 = (epsilon_temp2-0.01)*(epsilon_temp2-0.01);
      
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
      f2 = (epsilon_temp2-0.01)*(epsilon_temp2-0.01);
      
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
double beta_const){

    
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
    
    while(epsilon_temp2<0.01){
      
      nstar_temp=nstar_temp/2;
      
      epsilon_temp2=get_epsilon1( rstar2,epsilon, U_out,alpha_out,nstar_temp,  gammastar,t_star, mu_const);
      
    }
    
    
    double upper_bound=nstar;
    double lower_bound=nstar_temp;
    
    
    
    nstar=find_nstar(upper_bound,lower_bound,rstar2,
                     epsilon,U_out, alpha_out,gammastar,t_star,mu_const);

//    Rcpp::Rcout << "nstar:                                " << std::flush << nstar << std::endl;
    
        
    rstar1=rstar2;

    
    
//    Rcpp::Rcout << "rstar1:                                " << std::flush << rstar1 << std::endl;
    
    i=i+1;

    
    
    }
    
    
    
//    double rstar3=golden_r(rstar1,0,epsilon, U_out,alpha_out,nstar, gammastar,t_star, mu_const);
    
    
        
//    Rcpp::Rcout << "rstar_third:                                " << std::flush << rstar3 << std::endl;
    
        
//    Rcpp::Rcout << "nstar_Out:                                " << std::flush << nstar << std::endl;
//    Rcpp::Rcout << "rstar_Out:                                " << std::flush << rstar1 << std::endl;

            
    return nstar;

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


 
 Rcpp::Rcout << "Gamma Opt: Initial Gammastar"  << std::endl << gammastar  << std::endl ;
 Rcpp::Rcout << "Gamma Opt: Initial n"  << std::endl << val  << std::endl ;
 
    
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
    
    Rcpp::Rcout << "gammastar:Proposed"  << std::endl << gammastar  << std::endl ;
    
// Temporarily edit this out        
    val=get_n(gammastar,trace_const, lambda_star, epsilon1,  epsilon_converge,gamma_star_lower,mu_const,beta_const);

//    Rcpp::Rcout << "min-value"  << std::endl << min  << std::endl ;
      Rcpp::Rcout << "proposed val at upper"  << std::endl << val  << std::endl ;


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
    
    
    Rcpp::Rcout << "lower_bound  "  <<  std::endl << lower_bound <<std::endl;
    Rcpp::Rcout << "upper_bound  "  <<  std::endl << upper_bound <<std::endl;
    Rcpp::Rcout << "Value at current min  "  <<  std::endl << min <<std::endl;
    Rcpp::Rcout << "Value at upper_bound  "  <<  std::endl << val <<std::endl;
    
    
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

//          Rcpp::Rcout << "lower_bound  "  << std::endl <<  lower_bound <<std::endl;
//          Rcpp::Rcout << "upper_bound  "  << std::endl <<  upper_bound <<std::endl;
      
            
      }

    
    gammastar=  (g1+g2)/2;
    
    min = get_n(gammastar,trace_const, lambda_star, epsilon1,  epsilon_converge,gamma_star_lower,mu_const,beta_const);
    
    Rcpp::Rcout << "min: Out of Golden Section search"  << std::endl << min  << std::endl ;
    
    

    // Return Final Estimate of gammastar

    Rcpp::List outlist=Rcpp::List::create(
    Rcpp::Named("gammastar")=gammastar,
            Rcpp::Named("nstar")=min);  


    return(outlist);

  
}





Rcpp::List set_nstar(NumericMatrix x, NumericMatrix P, NumericMatrix P_0,
arma::vec mu_0, arma::vec mu_star,arma::vec beta_star,arma::vec beta_star2,double epsilon_converge){


    int l1=x.ncol();

    arma::mat P2(P.begin(),P.nrow(),P.ncol(),false);
    arma::mat x2(x.begin(), x.nrow(), x.ncol(), false);
    arma::mat P_0b(P_0.begin(), P_0.nrow(), P_0.ncol(), false);
    arma::mat PX=P2*x2;
    arma::mat XTPX=x2.t()*PX;
    arma::mat P_Inner=P_0b+XTPX;
    arma::vec mu_diff=mu_0-mu_star;
    arma::vec mu_const=0.5*mu_diff.t()*XTPX*mu_diff;
    arma::vec beta_diff=beta_star2-beta_star;
    
    arma::vec beta_const=0.5*beta_diff.t()*PX*inv_sympd(XTPX)*PX.t()*beta_diff;
      
    arma::vec eigval;
    arma::mat eigvec;
    
    eig_sym(eigval, eigvec, XTPX);


    arma::mat eigvec2=eigvec.t();

      for(int i=0;i<l1;i++)    eigvec2.row(i)=eigvec2.row(i)/sqrt(eigval(i));


    arma::mat InvXTPX_1_2=eigvec*eigvec2;

    arma::mat P_AA=InvXTPX_1_2*P_Inner*InvXTPX_1_2;
    arma::mat P_AB=-InvXTPX_1_2*PX.t();
    arma::mat P_BA=P_AB.t();
    arma::mat P_BB=P2;
    arma::mat P2_AB=inv_sympd(P_AA)*P_AB*inv_sympd(P_BB)*P_BA;
    

    arma::vec  eigen_out=eig_sym(P2_AB.t()*P2_AB) ;
    
    double lambda_star=sqrt(eigen_out(l1-1));

    arma::mat P_Initial=P_0b+XTPX;    
    
        
        
    arma::mat P_Upper=P_0b+XTPX;    
    arma::mat P_Lower=P_Initial-PX.t()*inv_sympd(P2+PX*inv_sympd(P_Initial)*PX.t())*PX;
  
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
    
        
    
    
    eigvec2=eigvec.t();
    
    for(int i=0;i<l1;i++)    eigvec2.row(i)=eigvec2.row(i)*sqrt(eigval(i));

    arma::mat W_1_2=eigvec*eigvec2;
    
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
    

    double t_star=sqrt(beta_const(0)/gammastar);
    double gammastar2=(1+t_star)*gammastar;
    double trace_const2=(1+t_star)*trace_const;
    double mu_const2=(1+t_star)*mu_const(0);
    double beta_const2=0;
    if(t_star>0) beta_const2=beta_const(0)/t_star;
    double gamma_star_lower2=(1+t_star)*gamma_star_lower;
      
    
//    double alpha_out=(1+gammastar)/(1+trace_const+lambda_star*gammastar);
//    double U_out=(1+trace_const+2*lambda_star*gammastar);
//    double epsilon=exp(log(epsilon1)-gammastar);
//    double A1_out=(1-exp(log(epsilon1)-gammastar));  
//    double rstar1=log(alpha_out)/(log(U_out)+log(alpha_out)-log(A1_out));
//    double A3=pow(A1_out,rstar1);
//    double log_A3=log(A3);
//    double log_A3_2=-rstar1*exp(log(epsilon1)-gammastar);

    double alpha_out=(1+gammastar2)/(1+trace_const2+lambda_star*gammastar2);
    double U_out=(1+trace_const2+2*lambda_star*gammastar2);
    double epsilon=exp(log(epsilon1)-beta_const2-gammastar2);
//    double A1_out=(1-exp(log(epsilon1)-beta_const2-gammastar2));
    double lg_A1_out=log(1-exp(log(epsilon1)-beta_const2-gammastar2));
    
    double x1=exp(log(epsilon1)-beta_const2-gammastar2);
    double qc_lg_A1_out=-x1-0.5*x1*x1;
  
    if(lg_A1_out>-2.47036e-012)
    {
      lg_A1_out=-2.47036e-012;
      if(qc_lg_A1_out>-2.47036e-012)  lg_A1_out=qc_lg_A1_out;
      
    }
    
    double lg_alpha_out=log(1+(1/gammastar2))-log(lambda_star+( (1+trace_const2)/gammastar2));
    double x3=((lambda_star-1)+ (trace_const2)/gammastar2);
    
    double qc_lg_alpha_out2= -x3+0.5*x3*x3;
    
    
          Rcpp::Rcout << "U_out"  << std::endl << U_out  << std::endl ;
          Rcpp::Rcout << "alpha_out"  << std::endl << alpha_out  << std::endl ;
          
          
    //      Rcpp::Rcout << "x1"  << std::endl << x1  << std::endl ;
    //    Rcpp::Rcout << "x2"  << std::endl << x2  << std::endl ;
    //    Rcpp::Rcout << "lg_alpha_out"  << std::endl << lg_alpha_out  << std::endl ;
    
    
    if (lg_alpha_out<5.00028e-014){
      
      lg_alpha_out=5.00028e-014;
      if(qc_lg_alpha_out2<5.00028e-014) lg_alpha_out=qc_lg_alpha_out2;
      
    } 
    
  
    
//    double rstar1=log(alpha_out)/(log(U_out)+log(alpha_out)-log(A1_out));
//    double A3=pow(A1_out,rstar1);
//    double log_A3=log(A3);
//    double log_A3_2=-rstar1*exp(log(epsilon1)-beta_const2-gammastar2);


    double rstar1=lg_alpha_out/(log(U_out)+lg_alpha_out-lg_A1_out);
//    double A3=pow(A1_out,rstar1);
    double log_A3=rstar1*lg_A1_out;
    double log_A3_2=-rstar1*exp(log(epsilon1)-beta_const2-gammastar2);
    
    
    

// Adjustment in case log_A3=0 
// May want to use even if log_A3 if only close to zero

    if(log_A3==0){
          log_A3=log_A3_2;
    }

    if (log_A3==0){
      Rcpp::stop("Function Returned positive infinity");
      
    } 
    
    
    
    
//      W_1_2

//    Rcpp::Rcout << "gamma_star_lower" << std::endl << gamma_star_lower << std::endl;


//    double nstar=((log(epsilon_converge))-log(2+gammastar))/log_A3;
//    double nstar=((log(epsilon_converge))-log(2+gamma_star_lower+mu_const(0)))/log_A3;


// Initialize nstar 

    double nstar=((log(epsilon_converge))-log(2+gamma_star_lower2+mu_const2))/log_A3;
    
    double rstar2=golden_r(rstar1,0,epsilon, U_out,alpha_out,nstar, gammastar,t_star, mu_const(0));
    

    double nstar_temp=nstar;
    
//    get_epsilon1
    double epsilon_temp=get_epsilon1( rstar2,epsilon, U_out,alpha_out,nstar,  gammastar,t_star, mu_const(0));

    double epsilon_temp2=epsilon_temp;

    while(epsilon_temp2<0.01){
    
    nstar_temp=nstar_temp/2;
        
    epsilon_temp2=get_epsilon1( rstar2,epsilon, U_out,alpha_out,nstar_temp,  gammastar,t_star, mu_const(0));
    
    }
    

    double upper_bound=nstar;
    double lower_bound=nstar_temp;
  
  Rcpp::Rcout << "rstar_first:                                " << std::flush << rstar1 << std::endl;
  Rcpp::Rcout << "nstar_first:                                " << std::flush << nstar << std::endl;
  
  
  nstar=find_nstar(upper_bound,lower_bound,rstar2,
                    epsilon,U_out, alpha_out,gammastar,t_star,mu_const(0));

  rstar1=rstar2;

  
    Rcpp::Rcout << "rstar_Second:                                " << std::flush << rstar1 << std::endl;
    Rcpp::Rcout << "nstar_Second:                                " << std::flush << nstar << std::endl;
    
  
    


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


void progress_bar2(double x, double N)
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

  for(int i=0;i<1000;i++){

    Rcpp::checkUserInterrupt();


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


// [[Rcpp::export]]

Rcpp::List rglmb_rand_cpp(int n,NumericVector y,NumericMatrix x, 
NumericVector mu,NumericMatrix P_0,NumericMatrix P,
NumericVector offset2,NumericVector wt,double dispersion,Rcpp::List
famfunc, Function f1,Function f2,Function f3,NumericVector start,
      std::string family="binomial",
      std::string link="logit",
      int Gridtype=2,      
      double epsilon_converge=0.01
                            ) {
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function asVec("as.vector");
  Rcpp::Function asDob("as.double");
//  Rcpp::Function set_lambda("set_nstar");

  unsigned int l1=x.ncol();
  unsigned int l2=x.nrow();
  unsigned int i;
  unsigned int j;
  int k;
  double dispersion2;
  NumericVector alpha(l2);

  NumericMatrix mutemp(l2,1);   
  NumericMatrix xtemp(1,1); 
  std::fill(xtemp.begin(), xtemp.end(), 1);
  
  arma::mat x2(x.begin(), x.nrow(), x.ncol(), false);
  arma::vec alpha2(alpha.begin(),l2,false);  
  arma::vec offset2b(offset2.begin(),l2,false);  
  arma::vec start2(start.begin(),l1,false);
  arma::mat mutemp2(mutemp.begin(), mutemp.nrow(), mutemp.ncol(), false);     
  
  mutemp2=x2*start2;
  
  NumericMatrix betatemp(l2,1);
  NumericVector alphatemp(l1);
  NumericVector offset3(l2);
  NumericVector wt3(l2);
  NumericVector temp(1);
  NumericVector temp2(1);
  NumericMatrix Ptilde=clone(P);
  Ptilde(0,0)=2*P(0,0);

  Rcpp::List out1;
  Rcpp::List out2;
  std::fill(wt3.begin(), wt3.end(), 1);

  if(family=="poisson"||family=="binomial") dispersion2=1;
  else dispersion2=dispersion;
  
    NumericVector  wt2=wt/dispersion2;

    NumericVector parin(1);
    NumericVector mu1=mu-mu;
    Rcpp::Function optfun("optim");
    arma::vec y2(y.begin(),l2, false);
    
    
//    arma::vec ystar=log(y2/(1-y2));
    NumericVector qc1(l2);
    arma::vec qc(qc1.begin(),l2,false);
    
//    ystar.print("ystar:");

    NumericVector low(l2);
    NumericVector high(l2);
    arma::vec low2(low.begin(),l2, false);
    arma::vec high2(high.begin(),l2, false);

    List opt;
    NumericVector b2a(1);
    NumericVector b2b(l2);

    arma::vec b2(b2b.begin(), l2, false);
    arma::mat P_0b(P_0.begin(), P_0.nrow(), P_0.ncol(), false);
    arma::mat P2(P.begin(), P.nrow(), P.ncol(), false);
    arma::vec mu2(mu.begin(), l1);


//    NumericMatrix LL(nstar4,0);

//    NumericMatrix betaout(n,l2);
//    NumericMatrix alphaout(n,l1);


// Find posterior mode 

    arma::mat P_Post=P2(0,0)*x2.t()*x2+P_0b;
    arma::mat Var_Post=inv_sympd(P_Post);    
    arma::vec mu_star2=Var_Post*(P2(0,0)*x2.t()*b2+P_0b*mu2);    
      
  for(k=0;k<10;k++){

  for(j=0;j<l2;j++){

    opt=optfun(_["par"]=parin,_["fn"]=f2, _["gr"]=f3,_["y"]=asVec(y[j]),
    _["x"]=xtemp,_["mu"]=asMat(mutemp[j]),_["P"]=P,_["alpha"]=asVec(alpha[j]),
    _["wt"]=asVec(wt2[j]),_["method"]="BFGS",
        _["hessian"]=true);

    b2a=asVec(opt[0]);
    b2(j)=b2a(0);  

  }

  mu_star2=Var_Post*(P2(0,0)*x2.t()*b2+P_0b*mu2);
  mutemp2=x2*mu_star2;

  }

//      b2.print("b2_Opt:");
  
    //  Calculate Initial Constants for Convergence rate calculations

    NumericMatrix P_rand(l2,l2);
    P_rand.fill_diag(P(0,0));

    Rcpp::Rcout << " " << std::endl;
    Rcpp::Rcout << "############  Deriving Initial Constants for Convergence Bounds (assuming symmetry) ############" << std::endl;
    Rcpp::Rcout << " " << std::endl;

List beta_stars=set_beta_const(x2,offset2,mu_star2,
1,y,xtemp,mutemp,P,offset2,wt,
dispersion,famfunc,f1,f2,f3,betatemp,x,mu,P_0,offset3,wt3,
family=family,link=link,
Gridtype=Gridtype);

arma::vec beta_star=beta_stars(0);
arma::vec beta_star2=beta_stars(1);

    Rcpp::Rcout << " " << std::endl;
    Rcpp::Rcout << "############  Initial Constants for Convergence Bounds (assuming symmetry) ############" << std::endl;

    Rcpp::Rcout.precision(5);

    
    List nstar_lambda_star=set_nstar(x, P_rand,P_0,mu_star2,mu_star2,beta_star,beta_star2,epsilon_converge);


    NumericVector temp3=nstar_lambda_star(0);
    double nstar=temp3(0);

    
    
    if (nstar-1>LLONG_MAX){
      Rcpp::stop("Number of Required Burn-in Iterations Exceeds Integer Long Long Limitations");
      
    } 
    
    
    temp3=nstar_lambda_star(1);
    
  
    double tau=temp3(0);
  
    Rcpp::List simconstants=nstar_lambda_star(2);
  
    
    double uppern;
    uppern=nstar+n*tau+1;
    

    unsigned long long int nstar2 = nstar+0.5;    
    unsigned long long int nstar3=n*tau+0.5;    
    unsigned long long int nstar4=nstar2+nstar3;

//    long long int nstar4=nstar+n*tau+1;
    
    
//    Rcpp::Rcout << "LLONG_MAX:         " << std::flush << LLONG_MAX << std::endl;
//    Rcpp::Rcout << "ULONG_MAX:         " << std::flush << ULONG_MAX << std::endl;
    
    
    if (uppern>LLONG_MAX){
      Rcpp::stop("Number of Total Iterations Exceeds Integer Long Long Limitations");
      
    } 
    
    unsigned long long int Ulong_Max2;
    Ulong_Max2=2*ULONG_MAX;
    
    if (nstar4>Ulong_Max2){
      Rcpp::stop("Number of Total Iterations Exceeds Integer Limitations");
      
    } 
    
    
    
//    Rcpp::Rcout << "nstar:         " << std::flush << nstar << std::endl;
    
    Rcpp::Rcout.precision(5);

    Rcpp::Rcout << "############ Initial Estimate of Number of Required iterations  #######################" << std::endl;
    Rcpp::Rcout << "" << std::endl;
    Rcpp::Rcout << "Burn-in iterations:         " << std::flush << nstar2 << std::endl;
    Rcpp::Rcout << "Post Burn-in iterations:    " << std::flush << nstar3 << std::endl;
    Rcpp::Rcout << "Total Iterations:    " << std::flush << nstar4 << std::endl;
    

    NumericMatrix betaout(nstar4,l2);
    NumericMatrix alphaout(nstar4,l1);
    NumericMatrix betastarout(1000,l2);
    NumericVector LL(nstar4);
    NumericVector LL2(nstar3);

  
    
  for(i=0;i<nstar4;i++){

    Rcpp::checkUserInterrupt();


  if(i==0) {
    Rcpp::Rcout << "Running burnin simulation:" << std::endl;}
  
 

  if(i<nstar2) {  
    progress_bar2(i, nstar2-1);}
  
  
  if(i==nstar2) {Rcpp::Rcout << "" << std::endl
  << "Running main simulation:" << std::endl;}
  if(i>(nstar2-1)) {  progress_bar2(i-nstar2, nstar4-nstar2);
  }
  if(i==nstar4-1) {Rcpp::Rcout << "" << std::endl;}



    alpha2=mutemp2+offset2b;
    

  for(j=0;j<l2;j++){

      if(i==0){  parin=  asVec(0);}
      if(i>0){parin=asVec(betatemp(j-1,0)-mutemp(j,0));}

    if(i==0){   
        out1=glmbsim_NGauss2_cpp(1,asVec(y[j]),xtemp,
                                 asMat(mutemp[j]),P,
                                 asVec(offset2[j]),
                                 asVec(wt[j]),
                                 dispersion,
                                 famfunc,f1,f2,f3,
                                 asMat(mutemp[j]),
                                 family=family,
                                 link=link,
                                 Gridtype=Gridtype);
                                 }
    if(i>0){   
        out1=glmbsim_NGauss2_cpp(1,asVec(y[j]),xtemp,
                                 asMat(mutemp[j]),P,
                                 asVec(offset2[j]),
                                 asVec(wt[j]),
                                 dispersion,
                                 famfunc,f1,f2,f3,asMat(betatemp(j-1,0)),
                                 family=family,
                                 link=link,
                                 Gridtype=Gridtype);
                                 }



            temp=out1(0);
            betatemp(j,0)=temp(0);
           betaout(i,j)=temp(0);
           
        //   Note: Index to which temp 2 should be set depends on what 
        //  index LL has
           
        temp2=out1(1);
//        LL(i)=LL(i)+temp2(0);
        if(i>(nstar2-1)) LL2(i-nstar2)=LL2(i-nstar2)+temp2(0); 
//    Rcpp::Rcout << "Setting LL" << std::endl << LL(i,0) << std::endl;
        
      }

    // Should edit famfunc passed here....

    out2=glmbsim_Gauss_cpp(1,betatemp,x,
                        mu,P_0,offset3
                        ,wt3,
                        1/P(0,0),
                        famfunc,f1,f2,f3,
                        mu);   
    alphatemp=out2(0);    
    arma::vec alphatemp2(alphatemp.begin(), l1, false);
    
    for(j=0;j<l1;j++) {alphaout(i,j)=alphatemp(j);}

      mutemp2=x2*alphatemp2;
    
    
    }

NumericMatrix alphaout2 = alphaout( Range(nstar2,nstar4-1),Range(0,l1-1));
NumericMatrix betaout2= betaout(Range(nstar2,nstar4-1),Range(0,l2-1));
  
//  NumericMatrix betaout(nstar4,l2);
  
  
    
Rcpp::Rcout << " " << std::endl;
Rcpp::Rcout << "############  Confirming convergence... ############" << std::endl;
Rcpp::Rcout << " " << std::endl;

// Calculate mustar

arma::vec mu_0=mu_star2;
NumericVector b3(l2);
arma::vec mu_star3(b3.begin(), l1);


for(i=0;i<l1;i++){mu_star2(i)=mean(alphaout2(_,i));}

beta_stars=set_beta_const(x2,offset2,mu_star2,
1,y,xtemp,mutemp,P,offset2,wt,
dispersion,famfunc,f1,f2,f3,betatemp,x,mu,P_0,offset3,wt3,
family=family,link=link,
Gridtype=Gridtype);

arma::vec b_star=beta_stars(0);
arma::vec b_star2=beta_stars(1);


Rcpp::Rcout.precision(5);

Rcpp::Rcout << " " << std::endl;
Rcpp::Rcout << "############  Revised Constants for Convergence Bounds (Controlling for asymmetry) ############" << std::endl;


nstar_lambda_star=set_nstar(x, P_rand,P_0,mu_0,mu_star2,b_star,b_star2,epsilon_converge);

    temp3=nstar_lambda_star(0);
double nstarb=temp3(0);

    Rcpp::Rcout << "Old nstar:         " << std::flush << nstar << std::endl;
    Rcpp::Rcout << "New nstar:         " << std::flush << nstarb << std::endl;

if(nstar>=nstarb){
Rcpp::Rcout << "No additional iterations required:" <<  std::endl;

}

if(nstar<nstarb){
Rcpp::Rcout << "An additional " << std::flush << nstarb-nstar << std::flush  <<" iterations required:" <<  std::endl;
}


Rcpp::List Prior=Rcpp::List::create(Rcpp::Named("mean")=mu,
Rcpp::Named("Precision")=P_0);


Rcpp::List Out=Rcpp::List::create(Rcpp::Named("coefficients")=alphaout2,
Rcpp::Named("randcoefficients")=betaout2,
Rcpp::Named("PostMode")=mu_0,
Rcpp::Named("Prior")=Prior,
Rcpp::Named("iters")=1,
Rcpp::Named("famfunc")=famfunc,
Rcpp::Named("dispersion")=dispersion,
Rcpp::Named("loglike")=LL2,
Rcpp::Named("simconstants")=simconstants
);  
  
  
return(Out  ) ; 
  
}




