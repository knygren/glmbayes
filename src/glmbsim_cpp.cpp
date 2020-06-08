// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "famfuncs.h"
#include "Envelopefuncs.h"

using namespace Rcpp;

// Putside Function declarations (move or delete)

void  f4_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, NumericVector NegLL, NumericMatrix cbars, int progbar=0);


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

// [[Rcpp::export]]

Rcpp::List glmb_Standardize_Model(
    NumericVector y, 
    NumericMatrix x,   // Original design matrix (to be adjusted)
    NumericMatrix P,   // Prior Precision Matrix (to be adjusted)
    NumericMatrix bstar, // Posterior Mode from optimization (to be adjusted)
    NumericMatrix A1  // Precision for Log-Posterior at posterior mode (to be adjusted)
  ) {
  
  // Get dimensions for x matrix
  
  int l1=x.ncol();
  int l2=x.nrow();

  // Arma matrix versions of most inputs
    
  arma::mat x2(x.begin(), l2, l1, false);
  arma::mat P2(P.begin(), P.nrow(), P.ncol(), false);
  arma::mat b2(bstar.begin(), bstar.nrow(), bstar.ncol(), false);
  arma::mat A1_b(A1.begin(), l1, l1, false); 
  
  // For now keep these name (legacy from older code inside glmbsim_NGauss
  
  NumericMatrix A4_1(l1, l1);  
  NumericMatrix b4_1(l1,1);  // Maybe this can be set as vector to avoid the below conversi
  NumericMatrix x4_1(l2, l1);
  NumericMatrix mu4_1(l1,1);  
  NumericMatrix P5_1(l1, l1);
  NumericMatrix P6Temp_1(l1, l1);
  NumericMatrix L2Inv_1(l1, l1); 
  NumericMatrix L3Inv_1(l1, l1); 

  double scale=1;
  int check=0;
  double eigval_temp;
  
  Rcpp::Function asVec("as.vector");

  arma::mat L2Inv(L2Inv_1.begin(), L2Inv_1.nrow(), L2Inv_1.ncol(), false);
  arma::mat L3Inv(L3Inv_1.begin(), L3Inv_1.nrow(), L3Inv_1.ncol(), false);
  arma::mat b4(b4_1.begin(), b4_1.nrow(), b4_1.ncol(), false);
//  arma::mat mu4(mu4_1.begin(), mu4_1.nrow(), mu4_1.ncol(), false);
  arma::mat x4(x4_1.begin(), x4_1.nrow(), x4_1.ncol(), false);
  arma::mat A4(A4_1.begin(), A4_1.nrow(), A4_1.ncol(), false);
  arma::mat P5(P5_1.begin(), P5_1.nrow(), P5_1.ncol(), false);
  arma::mat P6Temp(P6Temp_1.begin(), P6Temp_1.nrow(), P6Temp_1.ncol(), false);
  arma::vec eigval_1;
  arma::mat eigvec_1;
  arma::vec eigval_2;
  arma::mat eigvec_2;
  arma::mat ident=arma::mat (l1,l1,arma::fill::eye);
  
// Standardize Model to Have Diagonal Variance-Covariance Matrix at Posterior Mode

      eig_sym(eigval_1, eigvec_1, A1_b);
      arma::mat D1=arma::diagmat(eigval_1);
      arma::mat L2= arma::sqrt(D1)*trans(eigvec_1);
      L2Inv=eigvec_1*sqrt(inv_sympd(D1));  // Also used to undo normalization later

      // output variables used in latter step

      arma::mat b3=L2*b2;   
//      arma::mat mu3=L2*mu2; // These are needed but will not be used to pass 
      //  arma::mat mu3=L2*mu1b; // Corrected - mu1b is zero vector

      arma::mat x3=x2*L2Inv;
      arma::mat P3=trans(L2Inv)*P2*L2Inv;

      // Seup for loop that sets epsilon
    
      arma::mat P3Diag=arma::diagmat(arma::diagvec(P3));// diagonal part of P3
      arma::mat epsilon=P3Diag;
      arma::mat P4=P3Diag;   

    // Find scaled matrix epsilon

    while(check==0){
      epsilon=scale*P3Diag;  // scaled version of diagonal matrix
  
      // Checks if difference between Prior precision and diagonal matrix
      // is positive definite
      // is positive definite - to be added to likelihood 
  
      P4=P3-epsilon;				
      eig_sym(eigval_2, eigvec_2, P4);
      eigval_temp=arma::min(eigval_2);
      if(eigval_temp>0){check=1;
    
      //      Rcpp::Rcout << "scale after step1 " << std::flush << scale << std::endl;
      //      P4.print("P4 after step 1");  
      //      epsilon.print("epsilon after step 1");  
      }
      else{scale=scale/2;}
    }

      // Setup prior to Eigenvalue decomposition

      arma::mat A3=ident-epsilon;	// This should be a diagonal matrix and represents "data" precision in transformed model

    //   Put into Standard form where prior is identity matrix

    // Perform second eigenvalue decomposition

    eig_sym(eigval_2, eigvec_2, epsilon);
    arma::mat D2=arma::diagmat(eigval_2);
    
    // Apply second transformation

    arma::mat L3= arma::sqrt(D2)*trans(eigvec_2);
    L3Inv=eigvec_2*sqrt(inv_sympd(D2));
    b4=L3*b3; 
//    mu4=L3*mu3; 
    x4=x3*L3Inv;
    A4=trans(L3Inv)*A3*L3Inv;  // Should be transformed data precision matrix
    P5=trans(L3Inv)*P4*L3Inv;  // Should be precision matrix without epsilon
    P6Temp=P5+ident;           // Should be precision matrix for posterior (may not be used)

    // "Oddball extra steps due to legacy code 
    
    NumericVector b5=asVec(b4_1); // Maybe this causes error?
    
    NumericMatrix mu5_1=0*mu4_1; // Does this modify mu4_1? Set mu5_1 to 0 more efficiently
    
    return Rcpp::List::create(
      Rcpp::Named("bstar2")=b5,       // Transformed posterior mode (untransposed also used)
      Rcpp::Named("A")=A4_1,                 // Precision for Standardized data precision
      Rcpp::Named("x2")=x4_1,                // Transformed Design matrix
      Rcpp::Named("mu2")=mu5_1,               // Transformed prior mean (should really always be 0)
      Rcpp::Named("P2")=P5_1,               // Precision matrix for Normal component shifted to Log-Likelihood
      Rcpp::Named("L2Inv")=L2Inv,               // Precision matrix for Normal component shifted to Log-Likelihood
      Rcpp::Named("L3Inv")=L3Inv               // Precision matrix for Normal component shifted to Log-Likelihood
    );
  
}

// [[Rcpp::export]]

Rcpp::List  glmbsim_cpp(int n,NumericVector y,NumericMatrix x,
NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,
Function f2,Rcpp::List  Envelope,Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, int progbar=1)
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

  if(progbar==1){
  Rcpp::Rcout << "Starting Simulation:" << std::endl;  };
    for(int i=0;i<n;i++){
      
      Rcpp::checkUserInterrupt();
     // if(progbar==1){
        progress_bar2(i, n-1);
      if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
      //}
      
      
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
    
    int l1b=mu.length();
    int l1c=P.ncol();
    int l1d=P.nrow();
    
    if(l1b!=l1) Rcpp::stop("Number of rows in mu not consistent with number of columns in matrix x");;
    if(l1c!=l1) Rcpp::stop("Number of columns in matrix P not consistent with number of columns in matrix x");;
    if(l1d!=l1) Rcpp::stop("Number of rows in matrix P not consistent with number of columns in matrix x");;
    
    int l2b=y.length();
    int l2c=offset2.length();
    int l2d=wt.length();
    
    if(l2b!=l2) Rcpp::stop("Number of rows in y not consistent with number of rows in matrix x");;
    if(l2c!=l2) Rcpp::stop("Number of rows in offset2 vector not consistent with number of rows in matrix x");;
    if(l2d!=l2) Rcpp::stop("Number of rows in wt vector not consistent with number of rows in matrix x");;
    
    
    double dispersion2;
    NumericVector alpha(l2);
    NumericMatrix mu2a=asMat(mu);
    NumericMatrix out(l1,n);   
    NumericVector LL(n);   
    
    
    arma::mat x2(x.begin(), l2, l1, false);
    arma::vec alpha2(alpha.begin(),l2,false);  
    arma::vec offset2b(offset2.begin(),l2,false);  
    arma::mat mu2(mu2a.begin(), mu2a.nrow(), mu2a.ncol(), false);

    NumericMatrix x2b(clone(x));
    arma::mat P2(P.begin(), P.nrow(), P.ncol(), false);

    if(family=="poisson"||family=="binomial")dispersion2=1;
	  else dispersion2=dispersion;

	  int i;  // This can be likely be shifted towards top of function
	 
    NumericVector  wt2=wt/dispersion2; // Adjusts weight for dispersion
    arma::vec wt3(wt2.begin(), x.nrow());

    // Step 1: Shifts mean and offset to alpha --> Modified offset - Set inputs for optimization
    // Transforme Model now has prior mean 0
    
    alpha2=x2*mu2+offset2b; 
    NumericVector parin=start-mu;  // Starting value for optimization is now start - mu
    NumericVector mu1=mu-mu;       // new prior means are zero
    Rcpp::Function optfun("optim");
  
    arma::vec mu1b(mu1.begin(),l2,false);
  
    // Step 2: Run posterior optimization with log-posterior function and gradient functions
    // glmsim_NGauss2 seems to use a function for optimization (optPostMode)
      
    List opt=optfun(_["par"]=parin,_["fn"]=f2, _["gr"]=f3,_["y"]=y,
    _["x"]=x,
    _["mu"]=mu1,_["P"]=P,_["alpha"]=alpha,_["wt"]=wt2,_["method"]="BFGS",_["hessian"]=true);

  
    NumericMatrix b2a=asMat(opt[0]);  // optimized value
    NumericVector min1=opt[1]; // Not clear this is used - should be minimum
    int conver1=opt[3]; // check on convergence
    
    // Approximate hessian - Should consider replacing with Hessian based on 
    // known Hessian formula (when available)
    // This could be a source of error 
    
    NumericMatrix A1=opt[5]; 

    // Return Error if Optimizaton failed
    
    if(conver1>0){Rcpp::stop("Posterior Optimization failed");}
  
   // Standardize the model 

   Rcpp::Rcout << "Standardizing the model:" << std::endl;
   
   Rcpp::List Standard_Mod=glmb_Standardize_Model(
       y, 
       x,   // Original design matrix (to be adjusted)
       P,   // Prior Precision Matrix (to be adjusted)
       b2a, // Posterior Mode from optimization (to be adjusted)
       A1  // Precision for Log-Posterior at posterior mode (to be adjusted)
   );

    // Get output from call to glmb_Standardize_Model (not sure if they really need to be allocated)     
   // Advantage of allocating may be due to clarity of code in below

     NumericVector bstar2_temp=Standard_Mod[0];
     NumericMatrix A_temp=Standard_Mod[1];
     NumericMatrix x2_temp=Standard_Mod[2];
     NumericMatrix mu2_temp=Standard_Mod[3];
     NumericMatrix P2_temp=Standard_Mod[4];
     arma::mat L2Inv_temp=Standard_Mod[5];
     arma::mat L3Inv_temp=Standard_Mod[6];

    // Calls seem to be different because of setting for gridsort component!
    // When doing samples of just 1, sorting grid is slow when running many samples,
    // using a sorted grid is much faster
    // Most applications where n=1 are likely to be Gibbs samplers
    // There are likely to be few instances where someone runs a small 
    // number of samples greater than 1

    Rcpp::Rcout << "Starting Envelope Creation:" << std::endl;

    Rcpp::List Envelope; // Can move this towards top of the function
    
        
    if(n==1){
      Envelope=glmbenvelope_c(bstar2_temp, A_temp,y, x2_temp,mu2_temp,
                              P2_temp,alpha,wt2,family,link,Gridtype, n,false);
          }
    
    if(n>1){
      Envelope=glmbenvelope_c(bstar2_temp, A_temp,y, x2_temp,mu2_temp,
                              P2_temp,alpha,wt2,family,link,Gridtype, n,true);
          }
    
    Rcpp::Rcout << "Finished Envelope Creation:" << std::endl;

    // Run simulation 
    
    Rcpp::List sim=glmbsim_cpp(n,y,x2_temp,mu2_temp,P2_temp,alpha,wt2,
                               f2,Envelope,family,link);
    
        
    //  Post processing
    
    //  1) Undo-Standardization of Posterior Precision
    //  2) Undo shifting of prior mean to offset
    //  3) Calculate Log_likelihood (used in model diagnostics)

    // These two can likely be shifted towards top of function to make code simpler
    
    NumericMatrix sim2=sim[0];

    // Arma matrices pointing to matrices above 
    
    arma::mat sim2b(sim2.begin(), sim2.nrow(), sim2.ncol(), false);
    arma::mat out2(out.begin(), out.nrow(), out.ncol(), false);
    
    out2=L2Inv_temp*L3Inv_temp*trans(sim2b); // reverse transformation

    // Add mean back in and compute LL (for post processing)
    // Can add option to not compute LL
  
    for(i=0;i<n;i++){
    out(_,i)=out(_,i)+mu;  // Add mean vector back 
    LL[i]=as<double>(f1(_["b"]=out(_,i),_["y"]=y,_["x"]=x,offset2,wt2)); // Calculate log_likelihood
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

Rcpp::List rnorm_reg_cpp(int n,NumericVector y,NumericMatrix x, 
NumericVector mu,NumericMatrix P,NumericVector offset2,NumericVector wt,double dispersion,Rcpp::List
famfunc, Function f1,Function f2,Function f3,NumericVector start,
      std::string family="binomial",
      std::string link="logit",
      int Gridtype=2      
) {

    // add checks to make sure that dimensions are consistent
    // (i) number of rows in y, offset2, and wt2 should equal rows of x
    // (ii) number of rows in mu should equal number of columns in x
    // (iii) number of rows and columns in P should equal number of columns in x
    
      
    // Need to check combination of weighting and offset working properly

    Rcpp::Function asMat("as.matrix");
    Rcpp::Function asVec("as.vector");
    int l1=x.ncol();
    int l2=x.nrow();

    int l1b=mu.length();
    int l1c=P.ncol();
    int l1d=P.nrow();
    
    if(l1b!=l1) Rcpp::stop("Number of rows in mu not consistent with number of columns in matrix x");;
    if(l1c!=l1) Rcpp::stop("Number of columns in matrix P not consistent with number of columns in matrix x");;
    if(l1d!=l1) Rcpp::stop("Number of rows in matrix P not consistent with number of columns in matrix x");;
    
    int l2b=y.length();
    int l2c=offset2.length();
    int l2d=wt.length();
    
    if(l2b!=l2) Rcpp::stop("Number of rows in y not consistent with number of rows in matrix x");;
    if(l2c!=l2) Rcpp::stop("Number of rows in offset2 vector not consistent with number of rows in matrix x");;
    if(l2d!=l2) Rcpp::stop("Number of rows in wt vector not consistent with number of rows in matrix x");;
    
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
        
        // WTW apparently take on very large values
        
        arma::mat WTW=trans(W)*W;

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


















