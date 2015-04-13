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

    // Move inside If Statement once confirmed works
    
        int i;
        NumericVector  y1=y-offset2;
        NumericVector  y2=y1;
        arma::vec y2b(y2.begin(),l2,false);

    
        for(i=0;i<l2;i++){
          x2b(i,_) =x2b(i,_)*sqrt(wt2[i]);
          y2(i)=y2(i)*sqrt(wt2[i]);
          }        
        arma::mat RA=arma::chol(P2);

        NumericMatrix W1(l2+l1,l1);
        arma::mat W(W1.begin(), W1.nrow(), W1.ncol(), false);
        NumericVector z1(l2+l1);
        arma::vec z(z1.begin(),l2+l1,false);

        W.rows(0,l2-1)=x2;
        W.rows(l2,l2+l1-1)=RA;
        z.rows(0,l2-1)=y2b;
        z.rows(l2,l1+l2-1)=RA*mu2;
        arma::mat IR=arma::trans(arma::inv(trimatu(chol(trans(W)*W))));
      
//      arma::mat W;
//      .rows(first_row, last_row)
    
    NumericVector parin=start-mu;
    NumericVector mu1=mu-mu;
    Rcpp::Function optfun("optim");
    


    List opt=optfun(_["par"]=parin,_["fn"]=f2, _["gr"]=f3,_["y"]=y,_["x"]=x,
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

//     return(Rcpp::List::create(
//    Rcpp::Named("b5")=b5,
//    Rcpp::Named("A4_1")=A4_1,
//    Rcpp::Named("y")=y,
//    Rcpp::Named("x4_1")=x4_1,
//    Rcpp::Named("mu4_1")=mu4_1,
//    Rcpp::Named("P5_1")=P5_1,
//    Rcpp::Named("alpha")=alpha,
//    Rcpp::Named("wt2")=wt2,
//    Rcpp::Named("family")=family,
//    Rcpp::Named("link")=link,
//    Rcpp::Named("Gridtype")=Gridtype,
//    Rcpp::Named("n")=n
//));


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
    
    //alpha2=x2*mu2+offset2b;

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




