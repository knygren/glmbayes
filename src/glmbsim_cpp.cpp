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
      p1(j)=1-exp(-exb(j));
      p2(j)=exp(-exb(j));
      atemp(j)=exp(xb(j)-exb(j));

      xrow2=x2.row(j);
        Ptemp2=Ptemp2
        +wt(j)*atemp(j)*(
          (  ( y(j)*(1-exb(j))/p1(j) ) ) - ( y(j)*p2(j)*exb(j)/(p1(j)*p1(j)))
  +( (1-y(j))*(1-exb(j)) /p2(j) ) - ( (1-y(j))*exb(j)/p2(j)) )*trans(xrow2)*xrow2;

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
      xb2=exp(-alpha2- x2 * b2);
      for(j=0;j<l1;j++){
        xrow2=x2.row(j);
        Ptemp2=Ptemp2+wt(j)*xb(j)*trans(xrow2)*xrow2;
        }
  
    if(arma::is_finite(bstar2)){
      b2=inv_sympd(Ptemp2)*((Ptemp2-P2)*bstar2+P2*mu2); 
      xb2=exp(-alpha2- x2 * b2);
     
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
      xb2=exp(-alpha2- x2 * b2);
      for(j=0;j<l1;j++){
        xrow2=x2.row(j);
        Ptemp2=Ptemp2+wt(j)*xb(j)*trans(xrow2)*xrow2;
        }
  
    if(arma::is_finite(bstar2)){
      b2=inv_sympd(Ptemp2)*((Ptemp2-P2)*bstar2+P2*mu2); 
      xb2=exp(-alpha2- x2 * b2);
     
    bmu2=b2-mu2;
    double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    yy=-dpois_glmb(y,xb,true);    
        
    for(int j=0;j<l1;j++){
    yy[j]=yy[j]*wt[j];  
    }

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
      atemp(j)=exp(xb(j)-exb(j));

      xrow2=x2.row(j);
        Pout2=Pout2
        +wt(j)*atemp(j)*(
          (  ( y(j)*(1-exb(j))/p1(j) ) ) - ( y(j)*p2(j)*exb(j)/(p1(j)*p1(j)))
  +( (1-y(j))*(1-exb(j)) /p2(j) ) - ( (1-y(j))*exb(j)/p2(j)) )*trans(xrow2)*xrow2;

      }

        Varout=inv_sympd(Pout2);

        double res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
        yy=-dbinom_glmb(y,wt,xb,true);    
        res2=std::accumulate(yy.begin(), yy.end(), res1);

        for(int j=0;j<l1;j++){
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

      xb2=exp(-alpha2- x2 * b2);
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
        

        grad2=(P2 * bmu2+x2.t() * xb2);

    }

  /////////////////// quasipoisson /////////////////////////////
  
  if(family=="quasipoisson" )
  {

      xb2=exp(-alpha2- x2 * b2);
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
        

        grad2=(P2 * bmu2+x2.t() * xb2);

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
        btemp2=b2-stepsize*Varout*(P2 * bmu2+x2.t() * xb2);    
        bmutemp2=btemp2-mu2;

        xbtemp2=exp(-alpha2- x2 * btemp2);


        double res1=0.5*arma::as_scalar(bmutemp2.t() * P2 *  bmutemp2);
        yy=-dpois_glmb(y,xbtemp,true); 
        for(int j=0;j<l1;j++){
        yy[j]=yy[j]*wt[j];  }
        res3=std::accumulate(yy.begin(), yy.end(), res1);
    }

  /////////////////// quasipoisson /////////////////////////////
  
  if(family=="quasipoisson" )
  {
        btemp2=b2-stepsize*Varout*(P2 * bmu2+x2.t() * xb2);    
        bmutemp2=btemp2-mu2;

        xbtemp2=exp(-alpha2- x2 * btemp2);


        double res1=0.5*arma::as_scalar(bmutemp2.t() * P2 *  bmutemp2);
        yy=-dpois_glmb(y,xbtemp,true); 
        for(int j=0;j<l1;j++){
        yy[j]=yy[j]*wt[j];  }
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
      p1(j)=1-exp(-exb(j));
      p2(j)=exp(-exb(j));
      atemp(j)=exp(xbtemp(j)-exb(j));

      xrow2=x2.row(j);
        Pout2=Pout2
        +wt(j)*atemp(j)*(
          (  ( y(j)*(1-exb(j))/p1(j) ) ) - ( y(j)*p2(j)*exb(j)/(p1(j)*p1(j)))
  +( (1-y(j))*(1-exb(j)) /p2(j) ) - ( (1-y(j))*exb(j)/p2(j)) )*trans(xrow2)*xrow2;

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

    double res2=Initialize_bstar(y, x2, mu2,P2, alpha2,wt,b2, xb,xb2,
    Ptemp2,bmu2,bstar2,yy,family,link);


    ///////////////////////////////////////////////////////
    
       
    // Initialize while loop
    
    int i=0;
    int reset=0;
    int check=0;
    int check2=0;
    
    while(i<30 && check==0){

    /////////////////////////////////////////////////////////////

    // Calculate Function Value and gradient At Latest Iteration 
      
    
      Pout2=P2;

      res2=Find_Value(y,x2, mu2, P2,alpha2,  wt,  b2, xb,yy,grad2,Pout2,Varout,
      xb2,bmu2,xbtemp2,family,link);
      res_final=res2;
      
      reset=0;
 
    if(arma::any(grad2)==false){
    check=1;
    }
   
    gradb2=inv_sympd(P2)*grad2;
    gradb2=abs(2*gradb2/(2*b2+gradb2));
    maxgrad=max(gradb2);
    if(maxgrad<0.0001){
      check=1;
    }
    
    
    
    
    // Update b2 (Newton-Rhapson update)
    // Reduce stepsize if function value increases
    
    stepsize=1;
    
    check2=0;
    k=0;
    
    while(check2==0&& k<10 && check==0){

    //////////////   Set candidate point and check function value 

    res3=set_candidate( b2,  stepsize, Pout2, Varout,P2, bmu2, alpha2, 
    x2,xb2, mu2, btemp2, bmutemp2, xbtemp2, y, wt, xbtemp, yy,res2,family,link);


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
      
      if(reset==1){set_Pout(b2,y,alpha2,l1,P2,x2,wt,xbtemp,xbtemp2,xrow2,Pout2,family,link);}

    
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
    
    List opt=optfun(_["par"]=parin,_["fn"]=f2, _["gr"]=f3,_["y"]=y,
    _["x"]=x,
    _["mu"]=mu1,_["P"]=P,_["alpha"]=alpha,_["wt"]=wt2,_["method"]="BFGS",_["hessian"]=true);

 

    NumericMatrix b2a=asMat(opt[0]);
    arma::mat b2(b2a.begin(), b2a.nrow(), b2a.ncol(), false);

    NumericVector min1=opt[1];
    int conver1=opt[3];
    NumericMatrix A1=opt[5];

    
    if(conver1>0){Rcpp::stop("Posterior Optimization failed");}
    
    arma::mat A1_b(A1.begin(), l1, l1, false); 
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


//   return(Rcpp::List::create(Rcpp::Named("mean")=mu,Rcpp::Named("Precision")=P));

//    return(Envelope);

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
        LL[i]=as<double>(f1(_["b"]=out(i,_),_["y"]=y,_["x"]=x,offset2,wt2));

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

      NumericMatrix b2a(l1);
      NumericVector parin2(clone(parin));
      List opt1=optPostMode(y,x,mu1, P, alpha,wt2,
      parin2,log(y/(1-y)),family,link);

      b2a=asMat(opt1(0));
      NumericMatrix A1=opt1(1);
    NumericVector min1=asVec(opt1[2]);
    int conver1=0;
    
    arma::mat b2(b2a.begin(), b2a.nrow(), b2a.ncol(), false);
//    arma::mat Atemp_b(Atemp.begin(), l1, l1, false); 
    arma::mat A1_b(A1.begin(), l1, l1, false); 
//      b2.print("b2 - New Optimization:");
//    A1_b.print("A1 -  New Optimization");


//    arma::vec parin2b(parin2.begin(),l1);
//    parin2b.print("New Optimization:");


//    List opt=optfun(_["par"]=parin,_["fn"]=f2, _["gr"]=f3,_["y"]=y,_["x"]=x,
//    _["mu"]=mu1,_["P"]=P,_["alpha"]=alpha,_["wt"]=wt2,_["method"]="BFGS",_["hessian"]=true);


//    b2a=asMat(opt[0]);
//    arma::mat b2(b2a.begin(), b2a.nrow(), b2a.ncol(), false);
//    b2.print("Old Optimization:");
    
//    NumericVector min1=opt[1];
//    min1=opt[1];
//    int conver1=opt[3];
//    A1=asMat(opt[5]);

//    NumericMatrix  A1=opt[5];
//    A1=asMat(opt[5]);

    if(conver1>0){Rcpp::stop("Posterior Optimization failed");}




//    arma::mat A1_b(A1.begin(), l1, l1, false); 
//    b2.print("b2 - Old Optimization");
//    A1_b.print("A1 -  Old Optimization");
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
    
    for(i=0;i<n;i++){
    out(_,i)=out(_,i)+mu;
    }
      
  Rcpp::List outlist=Rcpp::List::create(Rcpp::Named("coefficients")=trans(out2),
  Rcpp::Named("Envelope")=Envelope);  

    return(outlist);
}






// [[Rcpp::export]]

Rcpp::List rglmb_rand_cpp(int n,NumericVector y,NumericMatrix x, 
NumericVector mu,NumericMatrix P_0,NumericMatrix P,
NumericVector offset2,NumericVector wt,double dispersion,Rcpp::List
famfunc, Function f1,Function f2,Function f3,NumericVector start,
      std::string family="binomial",
      std::string link="logit",
      int Gridtype=2      
) {
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function asVec("as.vector");


  int l1=x.ncol();
  int l2=x.nrow();
  int i;
  int j;
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
  
  
  NumericMatrix betaout(n,l2);
  NumericMatrix alphaout(n,l1);
  NumericMatrix betatemp(l2,1);
  NumericVector alphatemp(l1);
  NumericVector offset3(l2);
  NumericVector wt3(l2);
  NumericVector temp(1);
  
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
    arma::vec ystar=log(y2/(1-y2));
    NumericVector qc1(l2);
    arma::vec qc(qc1.begin(),l2,false);
    
//    ystar.print("ystar:");

    NumericVector low(l2);
    NumericVector high(l2);
    arma::vec low2(low.begin(),l2, false);
    arma::vec high2(high.begin(),l2, false);



    List opt;
  
  for(i=0;i<n;i++){

    alpha2=mutemp2+offset2b;
      qc=ystar-alpha2;
//    qc.print("qc:");  
    low=pmin(qc1,0);
    high=pmax(qc1,0);


//    low2.print("low2:");
//    high2.print("high2:");
    

  for(j=0;j<l2;j++){

  if(low(j)==high(j)){
    low(j)=-0.001;
    high(j)=0.001;
    
  }
  

if(i==0){  parin=  asVec(0);}

if(i>0){parin=asVec(betatemp(j-1,0)-mutemp(j,0));}

//if(arma::is_finite(low2(j)))
//{  opt=optfun(_["par"]=parin,_["fn"]=f2, _["gr"]=f3,_["y"]=asVec(y[j]),
//  _["x"]=xtemp,_["mu"]=asMat(mutemp[j]),_["P"]=P,_["alpha"]=asVec(alpha[j]),
//  _["wt"]=asVec(wt2[j]),_["method"]="Brent",_["lower"]=low[j],
//  _["upper"]=high[j],
//  _["hessian"]=true);
//}
//else{
  // Is this needed?
  
//opt=optfun(_["par"]=parin,_["fn"]=f2, _["gr"]=f3,_["y"]=asVec(y[j]),
//  _["x"]=xtemp,_["mu"]=asMat(mutemp[j]),_["P"]=P,_["alpha"]=asVec(alpha[j]),
//  _["wt"]=asVec(wt2[j]),_["method"]="BFGS",
//  _["hessian"]=true);



//}



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
            
      }
    
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
  
Rcpp::List Out=Rcpp::List::create(Rcpp::Named("betaout")=betaout,Rcpp::Named("alphaout")=alphaout,
Rcpp::Named("Envelope")=out1[1]);  
  
  
return(Out  ) ; 
  
}




