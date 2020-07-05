#' The Bayesian Normal-Gamma Regression Distribution
#'
#' \code{rnorm_gamma_reg} is used to generate iid samples from Bayesian linear models with a normal-gamma prior. 
#' The model is specified by providing a data vector, a design matrix, and 4 prior constants.
#'
#' 
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x a design matrix of dimension \code{m * p}.
#' @param prior_list a list with the prior parameters (mu, Sigma, shape and rate) for the 
#' prior distribution.
#' @param offset this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.
#' @param weights an optional vector of \sQuote{prior weights} to be used in the fitting process. Should be NULL or a numeric vector.
#' @inheritParams glmb
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
#' @family simfuncs 
#' @seealso \code{\link{pfamily}} for the list of available pfamilies and how 
#' they are specified. Each pfamily has a specific assigned simulation function. This
#' particular simulation function is the simulation function for the \code{\link{dNormal_Gamma}} 
#' pfamily.
#' @seealso \code{\link{rglmb}} and \code{\link{glmb}} for functions that depend on 
#' simulation functions.
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



rnorm_gamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian()){

## Added for consistency with earlier verion of function
  
offset2=offset
wt=weights

## Below code used precision matrix (not Sigma)
## Code checks for the presence of P in the prior
## if not present, it imputes Precision by inverting the Sigma matrix

if(missing(prior_list)) stop("Prior Specification Missing")
if(!missing(prior_list)){
  if(!is.null(prior_list$mu)) mu=prior_list$mu
  if(!is.null(prior_list$Sigma)) Sigma=prior_list$Sigma
  if(!is.null(prior_list$P)) P=prior_list$P
  if(is.null(prior_list$P)) P=solve(prior_list$Sigma)
  if(!is.null(prior_list$dispersion)) dispersion=prior_list$dispersion
  else dispersion=NULL
  if(!is.null(prior_list$shape)) shape=prior_list$shape
  else shape=NULL
  if(!is.null(prior_list$rate)) rate=prior_list$rate
  else rate=NULL
}

  

  
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
#V=V

#####    Beginning of rmultireg operations

#
# revision history:
#    Modified 4/7/15 to produce produce vectorized output  
#    changed 1/11/05 by P. Rossi to fix sum of squares error
#
# purpose:
#    draw from posterior for Multivariate Regression Model with
#    natural conjugate prior
# arguments:
#    Y is n x m matrix
#    X is n x k
#    Bbar is the prior mean of regression coefficients  (k x m)
#    A is prior precision matrix
#    nu, V are parameters for prior on Sigma
# output:
#    list of B, Sigma draws of matrix of coefficients and Sigma matrix
# model:
#    Y=XB+U  cov(u_i) = Sigma
#    B is k x m matrix of coefficients
# priors:  beta|Sigma  ~ N(betabar,Sigma (x) A^-1)
#                   betabar=vec(Bbar)
#                   beta = vec(B) 
#          Sigma ~ IW(nu,V) or Sigma^-1 ~ W(nu, V^-1)

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



P_Post=P+t(X)%*%X
mu_Post=solve(P_Post)%*%(mu+t(X)%*%Ytemp)

a_prior=shape     ## Should be relationship to shape in Wishart  
b_prior=rate  ## Should be relationship to scale in Wishart (could also be V/2)

a_post=a_prior+(nobs/2) # Posterior Shape parameter
b_post=b_prior+0.5*S # Posterior rate  (S is scaled differently than V?)

out1<-matrix(0,nrow=n,ncol=k)
dispersion=1/rgamma(n=n,shape=a_post,rate=b_post)
out2<-matrix(dispersion,nrow=n,ncol=1)

for(i in 1:n){out1[i,1:k]<- t(Btilde + IR%*%matrix(rnorm(m*k),ncol=m)*sqrt(dispersion[i])) }

draws<-matrix(1,n)
LL<-matrix(1,n)

for(i in 1:n){
  
  ## This function should return the negative log-likelihood 
  
  LL[i]=f1(b=out1[i,],y=y,x=x,alpha=offset2,wt=wt/out2[i])	
}


outlist<-list(
  coefficients=out1
  ,coef.mode=Btilde,
  dispersion=dispersion,
  Prior=list(mean=as.numeric(mu),Precision=P),
  prior.weights=wt,
  y=y,
  x=x,
  famfunc=famfunc,
  iters=draws,
  Envelope=NULL,
  loglike=LL)

colnames(outlist$coefficients)<-colnames(x)

outlist$call<-match.call()

class(outlist)<-c(outlist$class,"rglmb")

return(outlist)

#return(dispersion)

}


