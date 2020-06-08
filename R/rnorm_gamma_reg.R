#' The Bayesian Normal-Gamma Regression Distribution
#'
#' \code{rnorm_gamma_reg} is used to generate iid samples from Bayesian linear models with a normal-gamma prior. 
#' The model is specified by providing a data vector, a design matrix, and 4 prior constants.
#' @param n An integer
#' @param y A vector
#' @param x A matrix
#' @param mu A vector
#' @param P A matrix
#' @param nu A numeric constant>0
#' @param V A numeric constant>0
#' @param offset2 A vector
#' @param wt A vector or a numeric constant
#' @return The sum of \code{x} and \code{y}
#' @examples
#' 1+1
#' 10+1


rnorm_gamma_reg<-function(n,y,x,mu,P,nu,V,offset2=NULL,wt=1){

if(is.numeric(n)==FALSE||is.numeric(y)==FALSE||is.numeric(x)==FALSE||
is.numeric(mu)==FALSE||is.numeric(P)==FALSE) stop("non-numeric argument to numeric function")
  
  x <- as.matrix(x)
  mu<-as.matrix(as.vector(mu))
  P<-as.matrix(P)    
  
## Allow function to be called without offset2
  
if(length(n)>1) n<-length(n)	   
nobs <- NROW(y)
nvars <- ncol(x)

if(is.null(offset2)) offset2=rep(0,nobs)
nvars2<-length(mu)	
if(!nvars==nvars2) stop("incompatible dimensions")
if (!all(dim(P) == c(nvars2, nvars2))) 
  stop("incompatible dimensions")
if(!isSymmetric(P))stop("matrix P must be symmetric")
if(length(wt)==1) wt=rep(wt,nobs)
nobs2=NROW(wt)
nobs3=NROW(x)
nobs4=NROW(offset2)
if(nobs2!=nobs) stop("weighting vector must have same number of elements as y")
if(nobs3!=nobs) stop("matrix X must have same number of rows as y")
if(nobs4!=nobs) stop("offset vector must have same number of rows as y")

tol<- 1e-06 # Link this to Magnitude of P	
eS <- eigen(P, symmetric = TRUE,only.values = FALSE)
ev <- eS$values
if (!all(ev >= -tol * abs(ev[1L]))) 
  stop("'P' is not positive definite")

## Should add dimension checks here  
## Should move core part of rmultireg inside this code to avoid call

famfunc=glmbfamfunc(gaussian())  
f1=famfunc$f1
    
sim<-rmultireg(n=n,Y=matrix((y-offset2)*sqrt(wt),nrow=length(y)),X=x*sqrt(wt),Bbar=mu,A=P,nu=nu,V=V)


draws<-matrix(1,n)
LL<-matrix(1,n)

for(i in 1:n){
  
    ## This function should return the negative log-likelihood 
  
    LL[i]=f1(b=sim$B[i,],y=y,x=x,alpha=offset2,wt=wt/sim$Sigma[i])	
}

outlist<-list(coefficients=sim$B,PostMode=sim$BStar,
              Prior=list(mean=as.numeric(mu),Precision=P),
              iters=draws,famfunc=famfunc,Envelope=NULL,
              dispersion=sim$Sigma,loglike=LL)

colnames(outlist$coefficients)<-colnames(x)

outlist$call<-match.call()

class(outlist)<-c(outlist$class,"rglmb")

return(outlist)

}