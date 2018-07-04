// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;


void progress_bar2(double x, double N);


NumericVector dbinom_glmb( NumericVector x, NumericVector N, NumericVector means, int lg){
    int n = x.size() ;
    NumericVector res(n) ;
    for( int i=0; i<n; i++) res[i] = R::dbinom( round(x[i]*N[i]),round(N[i]), means[i], lg ) ;
    return res ;
}


///////////////////////// Logit Functions ///////////////////////////////////////

// [[Rcpp::export]]
NumericVector  f1_binomial_logit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);


    for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
  
    
    xb2=exp(-alpha2- x2 * b2);

    for(int j=0;j<l1;j++){
    xb(j)=1/(1+xb(j));
    }

    yy=-dbinom_glmb(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);

    }
    
    return res;      
}




// [[Rcpp::export]]
NumericVector  f2_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
      
    NumericVector invwt=1/sqrt(wt);

    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);

    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 

    double res1=0;


    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 

    bmu2=b2-mu2;
        
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
    xb2=exp(-alpha2- x2 * b2);

    for(int j=0;j<l1;j++){
    xb(j)=1/(1+xb(j));
    }


    yy=-dbinom_glmb(y,wt,xb,true);


    res(i) =std::accumulate(yy.begin(), yy.end(), res1);

    }
    
    return res;      
}




// [[Rcpp::export]]
arma::mat  f3_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat y2(y.begin(), l1, 1, false);
    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
       
   
       
    NumericMatrix Ptemp(l1,l1);  
      
    for(int i=0;i<l1;i++){
     Ptemp(i,i)=wt(i); 
    }  
    
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);
    NumericMatrix out(l2,m1);


    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 
    arma::mat Ptemp2(Ptemp.begin(), l1, l1, false);
    arma::mat out2(out.begin(), l2, m1, false);
    
    NumericMatrix::Column outtemp=out(_,0);
    //NumericMatrix::Row outtempb=out(0,_);

    arma::mat outtemp2(outtemp.begin(),1,l2,false);
    //arma::mat outtempb2(outtempb.begin(),1,l2,false);
    
    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
    if(progbar==1){ 
    progress_bar2(i, m1-1);
    if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
    };  
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);

    bmu2=b2-mu2;

//  		p<-1/(1+t(exp(-alpha-x%*%b)))
//		t(x)%*%((t(p)-y)*wt)+P%*%(b-mu)


    xb2=exp(-alpha2- x2 * b2);
    
    

    for(int j=0;j<l1;j++){
    xb(j)=1/(1+xb(j));  
    xb(j)=(xb(j)-y(j))*wt(j);
    }



    outtemp2= P2 * bmu2+x2.t() * xb2;
    }
    
   // return  b;
  
    return trans(out2);      
}


// [[Rcpp::export]]
arma::mat    f4_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, NumericVector NegLL, int progbar=0)
{
  
  // Get dimensions of x - Note: should match dimensions of
  //  y, b, alpha, and wt (may add error checking)
  
  // May want to add method for dealing with alpha and wt when 
  // constants instead of vectors
  
  int l1 = x.nrow(), l2 = x.ncol();
  int m1 = b.ncol();
  
  //    int lalpha=alpha.nrow();
  //    int lwt=wt.nrow();
  
  Rcpp::NumericMatrix b2temp(l2,1);
  
  arma::mat y2(y.begin(), l1, 1, false);
  arma::mat x2(x.begin(), l1, l2, false); 
  arma::mat alpha2(alpha.begin(), l1, 1, false); 
  
  Rcpp::NumericVector xb(l1);
  arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below

  
  NumericMatrix Ptemp(l1,l1);  
  
  for(int i=0;i<l1;i++){
    Ptemp(i,i)=wt(i); 
  }  
  
  // Moving Loop inside the function is key for speed
  
  NumericVector yy(l1);
  //  NumericVector res(m1);
  NumericMatrix bmu(l2,1);
  NumericMatrix out(l2,m1);
  
  
  arma::mat mu2(mu.begin(), l2, 1, false); 
  arma::mat bmu2(bmu.begin(), l2, 1, false); 
  arma::mat P2(P.begin(), l2, l2, false); 
  arma::mat Ptemp2(Ptemp.begin(), l1, l1, false);
  arma::mat out2(out.begin(), l2, m1, false);
  
  NumericMatrix::Column outtemp=out(_,0);
  //NumericMatrix::Row outtempb=out(0,_);
  
  arma::mat outtemp2(outtemp.begin(),1,l2,false);
  //arma::mat outtempb2(outtempb.begin(),1,l2,false);
  
  double res1;
  
  for(int i=0;i<m1;i++){
    Rcpp::checkUserInterrupt();
    
    if(progbar==1){ 
      progress_bar2(i, m1-1);
      if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
    };  
    
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(),l2,l2,false);
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);
    
    bmu2=b2-mu2;
    
    // First part is for log-prior
    
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
    
    xb2=exp(-alpha2- x2 * b2);
    
    for(int j=0;j<l1;j++){
      xb(j)=1/(1+xb(j));  
    }
    
    yy=-dbinom_glmb(y,wt,xb,true);
    NegLL(i) =std::accumulate(yy.begin(), yy.end(), res1);
    
    for(int j=0;j<l1;j++){
      xb(j)=(xb(j)-y(j))*wt(j);
    }
    
    outtemp2= P2 * bmu2+x2.t() * xb2;
  }
  
  // return  b;
  
  
    return trans(out2);      
}



///////////////////////// Probit Functions ///////////////////////////////////////

// [[Rcpp::export]]
NumericVector  f1_binomial_probit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);


    for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
 
    xb2=alpha2+ x2 * b2;   
    xb=pnorm(xb,0.0,1.0);
	

    yy=-dbinom_glmb(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);

    }
    
    return res;      
}




// [[Rcpp::export]]
NumericVector  f2_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
      
    NumericVector invwt=1/sqrt(wt);

    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);

    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 

    double res1=0;


    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 

    bmu2=b2-mu2;
        
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
  
    xb2=alpha2+ x2 * b2;   
    xb=pnorm(xb,0.0,1.0);

    yy=-dbinom_glmb(y,wt,xb,true);


    res(i) =std::accumulate(yy.begin(), yy.end(), res1);

    }
    
    return res;      
}




// [[Rcpp::export]]
arma::mat  f3_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat y2(y.begin(), l1, 1, false);
    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
       
   
       
    NumericMatrix Ptemp(l1,l1);  
      
    for(int i=0;i<l1;i++){
     Ptemp(i,i)=wt(i); 
    }  
    
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);
    NumericMatrix out(l2,m1);
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector d1(l1);


    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 
    arma::mat Ptemp2(Ptemp.begin(), l1, l1, false);
    arma::mat out2(out.begin(), l2, m1, false);
    
    NumericMatrix::Column outtemp=out(_,0);
    //NumericMatrix::Row outtempb=out(0,_);

    arma::mat outtemp2(outtemp.begin(),1,l2,false);
    //arma::mat outtempb2(outtempb.begin(),1,l2,false);
    
    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);

    bmu2=b2-mu2;
    
    xb2=alpha2+ x2 * b2;
    p1=pnorm(xb,0.0,1.0);
    p2=pnorm(-xb,0.0,1.0);
    d1=dnorm(xb,0.0,1.0);
    
//    -t(x)%*%as.matrix(((y*dnorm(alpha+x%*%b)/p1)-(1-y)*dnorm(alpha+x%*%b)/p2)*wt)+P%*%(b-mu)


    for(int j=0;j<l1;j++){
    xb(j)=(y(j)*d1(j)/p1(j)-(1-y(j))*d1(j)/p2(j))*wt(j);    
    }


    outtemp2= P2 * bmu2-x2.t() * xb2;
    }
    
   // return  b;
  
    return trans(out2);      
}

///////////////////////// cLOGLOG FUNCTION ///////////////////////////////////////

// [[Rcpp::export]]
NumericVector  f1_binomial_cloglog(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);


    for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
 
    xb2=alpha2+ x2 * b2;   
    xb=exp(-exp(xb));
  
    for(int j=0;j<l1;i++){
    xb(j)=1-xb(j);
    }
  

    yy=-dbinom_glmb(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);

    }
    
    return res;      
}




// [[Rcpp::export]]
NumericVector  f2_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
      
    NumericVector invwt=1/sqrt(wt);

    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);

    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 

    double res1=0;


    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 

    bmu2=b2-mu2;
        
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
    xb2=alpha2+ x2 * b2;   
    
    for(int j=0;j<l1;j++){
    xb(j)=1-exp(-exp(xb(j)));
    }

    yy=-dbinom_glmb(y,wt,xb,true);


    res(i) =std::accumulate(yy.begin(), yy.end(), res1);

    }
    
    return res;      
}




// [[Rcpp::export]]
arma::mat  f3_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat y2(y.begin(), l1, 1, false);
    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
       
   
       
    NumericMatrix Ptemp(l1,l1);  
      
    for(int i=0;i<l1;i++){
     Ptemp(i,i)=wt(i); 
    }  
    
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);
    NumericMatrix out(l2,m1);
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector atemp(l1);


    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 
    arma::mat Ptemp2(Ptemp.begin(), l1, l1, false);
    arma::mat out2(out.begin(), l2, m1, false);
    
    NumericMatrix::Column outtemp=out(_,0);
    //NumericMatrix::Row outtempb=out(0,_);

    arma::mat outtemp2(outtemp.begin(),1,l2,false);
    //arma::mat outtempb2(outtempb.begin(),1,l2,false);
    
    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);

    bmu2=b2-mu2;
    
    xb2=alpha2+ x2 * b2;

    for(int j=0;j<l1;j++){
    p1(j)=1-exp(-exp(xb(j)));
    p2(j)=exp(-exp(xb(j)));
    atemp(j)=exp(xb(j)-exp(xb(j)));
    xb(j)=((y(j)*atemp(j)/p1(j))-((1-y(j))*atemp(j)/p2(j)))*wt(j);    
    }


    outtemp2= P2 * bmu2-x2.t() * xb2;
    }
    
  
    return trans(out2);      
}


