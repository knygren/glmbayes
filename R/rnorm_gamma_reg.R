#' The Bayesian Normal-Gamma Regression Distribution
#'
#' \code{rnorm_gamma_reg} is used to generate iid samples from Bayesian linear models with a normal-gamma prior. 
#' The model is specified by providing a data vector, a design matrix, and 4 prior constants.
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x a design matrix of dimension \code{m * p}.
#' @param mu a vector of length \code{p} giving the prior means of the variables in the design matrix.
#' @param P a positive-definite symmetric matrix of dimension \code{p * p} specifying the prior precision matrix of the variable.
#' @param nu Prior shape parameter for the dispersion parameter (gaussian model only).
#' @param V Prior rate parameter for the dispersion parameter (gaussian model only).
#' @param offset2 this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.
#' @param wt an optional vector of \sQuote{prior weights} to be used in the fitting process. Should be NULL or a numeric vector.
#' @return \code{rnorm_gamma_reg} returns a object of class \code{"rglmb"}.  The function \code{summary} 
#' (i.e., \code{\link{summary.rglmb}}) can be used to obtain or print a summary of the results.
#' The generic accessor functions \code{\link{coefficients}}, \code{\link{fitted.values}},
#' \code{\link{residuals}}, and \code{\link{extractAIC}} can be used to extract
#' various useful features of the value returned by \code{\link{rnorm_gamma_reg}}.
#' An object of class \code{"rglmb"} is a list containing at least the following components:
#' \item{coefficients}{a \code{n} by \code{length(mu)} matrix with one sample in each row}
#' \item{PostMode}{a vector of \code{length(mu)} with the estimated posterior mode coefficients}
#' \item{Prior}{A list with two components. The first being the prior mean vector and the second the prior precision matrix}
#' \item{iters}{an \code{n} by \code{1} matrix giving the number of candidates generated before acceptance for each sample.}
#' \item{famfunc}{an object of class \code{"famfunc"}}
#' \item{Envelope}{an object of class \code{"envelope"}  }
#' \item{dispersion}{an \code{n} by \code{1} matrix with simulated values for the dispersion}
#' \item{loglike}{a \code{n} by \code{1} matrix containing the negative loglikelihood for each sample.}
#' @details The \code{rglmb} function produces iid samples for Bayesian generalized linear 
#' models. It has a more minimialistic interface than than the \code{\link{glmb}} 
#' function. Core required inputs for the function include the data vector, the design  
#' matrix and a prior specification. In addition, the dispersion parameter must 
#' currently be provided for the gaussian, Gamma, quasipoisson, and quasibinomial 
#' families (future implementations may incorporate a prior for these into the 
#' \code{rglmb} function).
#' 
#' Current implemented models include the gaussian family (identity link function), the
#' poisson and quasipoisson families (log link function), the gamma family (log link 
#' function), as well as the binomial and quasibinomial families (logit, probit, and 
#' cloglog link functions).  The function returns the simulated Bayesian coefficients 
#' and some associated outputs.
#' 
#' For the gaussian family, iid samples from the posterior density is genererated using 
#' standard simulation procedures for multivariate normal densities. For all other 
#' families, the samples are generated using accept-reject procedures by leveraging the 
#' likelihood subgradient approach of Nygren and Nygren (2006). This approach relies on
#' tight enveloping functions that bound the posterior density from above. The Gridtype 
#' parameter is used to select the method used for determining the size of a Grid used 
#' to build the enveloping function. See \code{\link{EnvelopeBuild_c}} for details. 
#' Depending on the selection, the time to build the envelope and the acceptance rate 
#' during the simulation process may vary. The returned value \code{iters} contains the 
#' number of candidates generated before acceptance for each draw.
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


#coefficients2=cbind(sim$B,sim$Sigma)

outlist<-list(
  coefficients=sim$B
#  coefficients=coefficients2
  ,PostMode=sim$BStar,
              Prior=list(mean=as.numeric(mu),Precision=P),
              iters=draws,famfunc=famfunc,Envelope=NULL,
              dispersion=sim$Sigma,loglike=LL)

colnames(outlist$coefficients)<-colnames(x)

outlist$call<-match.call()

class(outlist)<-c(outlist$class,"rglmb")

return(outlist)

}