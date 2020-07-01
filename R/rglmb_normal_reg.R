#' The Bayesian Generalized Linear Model with Normal Prior Distribution
#'
#' \code{rglmb_norm_reg} is used to generate iid samples from Bayesian Generalized linear 
#' models with a normal prior. The model is specified by providing a data vector, a design matrix, 
#' and 2 prior constants (mu and Sigma) for the normal prior.
#'
#' 
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x a design matrix of dimension \code{m * p}.
#' @param prior_list a list with the two parameters (mu and Sigma) for the prior distribution.
#' @param offset this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.
#' @param weights an optional vector of \sQuote{prior weights} to be used in the fitting process. Should be NULL or a numeric vector.
#' @inheritParams glmb
#' @return \code{rglmb_norm_reg} returns a object of class \code{"rglmb"}.  The function \code{summary} 
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
#' @details The \code{rglmb_norm_reg} function produces iid samples for Bayesian generalized linear 
#' models. with multivariate-normal priors. Core required inputs for the function include the data vector, 
#' the design matrix and a prior specification (provided in a list when this function is called directly). Specifying
#' the \code{\link{pfamily}} as \code{\link{dNormal}} in the \code{\link{lmb}} or \code{\link{glmb}} is equivalent
#' to calling this function directly.
#' 
#' The function returns the simulated Bayesian coefficients and some associated outputs. While it is possible to 
#' call this function directly, it is generally recommended that the \code{\link{lmb}}, \code{\link{rlmb}}, \code{\link{glmb}} or \code{\link{rglmb}} functions be used 
#' instead as those functions have more documentation and the resulting objects come with more methods.
#' 
#' When the specified family is gaussian, the multivariate normal is the conjugate prior distribution for the 
#' likelihood function and the posterior distribution is therefore also a multivariate normal density. Fairly
#' standard procedures are applied in order to generate samples in that case. Currently, this includes using a Cholesky 
#' decomposition. Later implementations may switch to a QR decomposition to increase consistency with the 
#' \code{\link{lm}} function and to increase associated numerical accuracy.
#' 
#' When the specified \code{\link{family}} is any log-concave non-gaussian family, then the estimation uses the Likelihood
#' subgradient density approach of Nygren and Nygren (2006).  This approaches uses tangencies to the 
#' log-likelihood function in order to construct an enveloping function from which candidate draws are 
#' generated and then either accepted/rejected using accept/reject methods. The core C function performing
#' this simulation essentially goes through the following steps: 
#'  
#' 
#' 1) The model is standardized to have prior mean vector equal to 0 (i.e., offsets and any prior mean are combined
#' into a constant term).
#' 
#' 
#' 2) The posterior mode for this transformed model is found. Currently this uses the \code{\link{optim}} function. 
#' Later implementations may replace this with iteratively reweighted least squares (IWLS) to increase 
#' consistency with the \code{\link{glm}} function and to enhance numerical accuracy.
#' 
#' 
#' 3) The model is further standardized so that (a) the precision matrix at the posterior mode is diagonal and (b)
#' the prior variance-covariance matrix is the identity matrix (see \code{\link{glmb_Standardize_Model}} for details).
#' 
#' 
#' 4) An enveloping function is built for the the standardized model, containing constants needed during simulation (see
#' \code{\link{EnvelopeBuild}} for details).
#' 
#' 
#' 5) Samples for the standardized model are generated using accept-reject methods (see \code{\link{rnnorm_reg_std}} for
#' details).
#' 
#' 
#' 6) The output from the standardized model are transformed back to the original scale by reversing the two 
#' eigenvalue decompositions and by adding back the prior mean.  
#' 
#' @family simfuncs 

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

rglmb_norm_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian()){

## Added for consistency with earlier verion of function
  
offset2=offset
wt=weights

## Missing control variables (add option to pass these)
## Setting for default Gridtype might be important

Gridtype=2


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
  start <- mu
  
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

if (is.character(family)) 
  family <- get(family, mode = "function", envir = parent.frame())
if (is.function(family)) 
  family <- family()
if (is.null(family$family)) {
  print(family)
  stop("'family' not recognized")
}

okfamilies <- c("gaussian","poisson","binomial","quasipoisson","quasibinomial","Gamma")
if(family$family %in% okfamilies){
  if(family$family=="gaussian") oklinks<-c("identity")
  if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-c("log")		
  if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-c("logit","probit","cloglog")		
  if(family$family=="Gamma") oklinks<-c("log")		
  if(family$link %in% oklinks){
    
    famfunc<-glmbfamfunc(family)
    f1<-famfunc$f1
    f2<-famfunc$f2
    f3<-famfunc$f3
    #      f5<-famfunc$f5
    #      f6<-famfunc$f6
  }
  else{
    stop(gettextf("link \"%s\" not available for selected family; available links are %s", 
                  family$link , paste(sQuote(oklinks), collapse = ", ")), 
         domain = NA)
    
  }	
  
}		
else {
  stop(gettextf("family \"%s\" not available in glmb; available families are %s", 
                family$family , paste(sQuote(okfamilies), collapse = ", ")), 
       domain = NA)
}


if(family$family=="gaussian"){ 
  outlist<-.rnorm_reg_cpp(n=n,y=y,x=x,mu=mu,P=P,offset2=offset2,wt=wt,dispersion=dispersion,famfunc=famfunc,f1=f1,f2=f2,f3=f3,start=mu)
}
else{
  if(is.null(dispersion)){dispersion2=1}
  else{dispersion2=dispersion}
  outlist<-.rnnorm_reg_cpp(n=n,y=y,x=x,mu=mu,P=P,offset2=offset2,wt=wt,
dispersion=dispersion2,famfunc=famfunc,f1=f1,f2=f2,f3=f3,
start=start,family=family$family,link=family$link,Gridtype=Gridtype)
}


colnames(outlist$coefficients)<-colnames(x)

# include family in final list
rglmb_df=as.data.frame(cbind(y,x))
rglmb_f=DF2formula(rglmb_df)
rglmb_mf=model.frame(rglmb_f,rglmb_df)

outlist$family=family  
outlist$call<-match.call()
outlist$offset2<-offset2
outlist$formula<-rglmb_f
outlist$model<-rglmb_mf
outlist$data<-rglmb_df

class(outlist)<-c(outlist$class,c("rglmb","glmb","glm","lm"))
outlist

}


