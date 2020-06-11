#' The Bayesian Normal-Gamma Regression Distribution
#'
#' \code{rnorm_gamma_reg_temp} is used to generate iid samples from Bayesian linear models with a normal-gamma prior. 
#' The model is specified by providing a data vector, a design matrix, and 4 prior constants.
#'
#' 
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
#' @details The \code{rnorm_gamma} function produces iid samples for Bayesian generalized linear 
#' models from the gaussian family (identity link) with a conjugate multivariate normal-gamma prior
#' for the regression coefficients and the dispersion (variance).
#' See Raiffa and Schlaifer (1961) for details on conjugate priors. 
#' 
#' Core required inputs for the function include the data vector, the design  
#' matrix and a prior specification. The function returns the simulated Bayesian coefficients 
#' and some associated outputs. The iid samples from the posterior density is genererated using 
#' standard simulation procedures for multivariate normal and gamma distributions. 
#' 
#' @references 
#' Dobson, A. J. (1990)
#' \emph{An Introduction to Generalized Linear Models.}
#' London: Chapman and Hall.
#' 
#' Hastie, T. J. and Pregibon, D. (1992)
#' \emph{Generalized linear models.}
#' Chapter 6 of \emph{Statistical Models in S}
#' eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#' McCullagh P. and Nelder, J. A. (1989)
#' \emph{Generalized Linear Models.}
#' London: Chapman and Hall.
#' 
#' Nygren, K.N. and Nygren, L.M (2006)
#' Likelihood Subgradient Densities. \emph{Journal of the American Statistical Association}.
#' vol.101, no.475, pp 1144-1156.
#' 
#' Raiffa, Howard and Schlaifer, R (1961)
#' \emph{Applied Statistical Decision Theory.}
#' Boston: Clinton Press, Inc.
#' 
#' Venables, W. N. and Ripley, B. D. (2002)
#' \emph{Modern Applied Statistics with S.}
#' New York: Springer.
#' 
#' @example inst/examples/Ex_rnorm_gamma_reg.R


rnorm_gamma_reg_temp<-function(n,y,x,mu,P,nu,V,offset2=NULL,wt=1){

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

## Convert as in call to rmultireg    

Y=matrix((y-offset2)*sqrt(wt),nrow=length(y))
X=x*sqrt(wt)
Bbar=mu
A=P
nu=nu
V=V

#####    Beginning of rmultireg operations

l0=length(Y)
Ytemp=matrix(Y,ncol=1)

l1=nrow(Ytemp)
if(l0>l1) stop("Dimensions of y not correct")

m=1

k=ncol(X)
l2=nrow(X)

if(l2!=l1) stop("Dimensions of X and Y are inconsistent")

k1=dim(A)[1]
k2=dim(A)[2]
k3=length(Bbar)

if(k1!=k) stop("dimensions of X and A are inconsistent")
if(k2!=k) stop("dimensions of X and A are inconsistent")
if(k3!=k) stop("dimensions of X and Bbar are inconsistent")


#
# first draw Sigma
#
RA=chol(A)
W=rbind(X,RA)

Z=rbind(Ytemp,matrix(RA%*%Bbar,ncol=1)) ## Force dimension to 1

#   note:  Y,X,A,Bbar must be matrices!
IR=backsolve(chol(crossprod(W)),diag(k))
#                      W'W = R'R  &  (W'W)^-1 = IRIR'  -- this is the UL decomp!
Btilde=crossprod(t(IR))%*%crossprod(W,matrix(Z,ncol=1))   
#                      IRIR'(W'Z) = (X'X+A)^-1(X'Y + ABbar)
S=crossprod(Z-W%*%Btilde)
#                      E'E

out1<-matrix(0,nrow=n,ncol=k)

if(m==1){
  
  out2=matrix(0,nrow=n,ncol=1)  
}

if(m>1){
  
  out2<-vector("list", n)
  
}

for(i in 1:n){
  
  rwout=rwishart(nu+l1,chol2inv(chol(V+S)))
  #
  # now draw B given Sigma
  #   note beta ~ N(vec(Btilde),Sigma (x) Covxxa)
  #       Cov=(X'X + A)^-1  = IR t(IR)  
  #       Sigma=CICI'    
  #       therefore, cov(beta)= Omega = CICI' (x) IR IR' = (CI (x) IR) (CI (x) IR)'
  #  so to draw beta we do beta= vec(Btilde) +(CI (x) IR)vec(Z_mk)  
  #			Z_mk is m x k matrix of N(0,1)
  #	since vec(ABC) = (C' (x) A)vec(B), we have 
  #		B = Btilde + IR Z_mk CI'
  #
  out1[i,1:k]<- t(Btilde + IR%*%matrix(rnorm(m*k),ncol=m)%*%t(rwout$CI))
  
  if(m==1){
    
    out2[i,1]<-rwout$IW
  }
  else{
    out2[i]<-rwout$IW
  }
  
}

sim=list(B=out1,Sigma=out2,BStar=Btilde)

####    End of rmultireg operations  


#sim<-rmultireg(n=n,Y=matrix((y-offset2)*sqrt(wt),nrow=length(y)),X=x*sqrt(wt),Bbar=mu,A=P,nu=nu,V=V)


draws<-matrix(1,n)
LL<-matrix(1,n)

for(i in 1:n){
  
    ## This function should return the negative log-likelihood 
  
    LL[i]=f1(b=sim$B[i,],y=y,x=x,alpha=offset2,wt=wt/sim$Sigma[i])	
}


outlist<-list(
  coefficients=sim$B
  ,PostMode=sim$BStar,
              Prior=list(mean=as.numeric(mu),Precision=P),
              iters=draws,famfunc=famfunc,Envelope=NULL,
              dispersion=sim$Sigma,loglike=LL)

colnames(outlist$coefficients)<-colnames(x)

outlist$call<-match.call()

class(outlist)<-c(outlist$class,"rglmb")

return(outlist)

}