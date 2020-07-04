#' Prior Family Objects for Bayesian Models
#'
#' Prior family objects provide a convenient way to specify the details of the priors 
#' used by functions such as \code{\link{glmb}}. See the documentation for \code{\link{glmb}}
#' for the details of how such model fitting takes place.
#' @param object the function \code{pfamily} accesses the \code{pfamily} objects which
#' are stored within objects created by modelling functions (e.g., \code{glmb}).
#' @param mu a prior mean vector for the model coefficients
#' @param Sigma a prior Variance-Covariance matrix for the model coefficients
#' @param dispersion the dispersion to be assumed when it is not given a prior. Should be provided
#' when the Normal prior is for the \code{gaussian()}, \code{Gamma()}, \code{quasibinomial},
#' or \code{quasipoisson} families. The \code{binomial()} and \code{poisson()} families
#' do not have dispersion coefficients. 
#' @param shape the prior shape parameter used by the gamma component of the prior. 
#' The gamma distribution is used as a prior for the inverse dispersion coefficients.
#' @param rate the rate parameter used by the gamma component of the prior.
#' @param beta the regression coefficients to be assumed when it is not given a prior. 
#' Needs to be provided when the Gamma prior is used for the dispersion. This
#' specification is typically only used as part of Gibbs sampling where the beta and 
#' dispersion parameters are updated separately. 
#' @param x an object, a family function that is to be printed
#' @param \ldots additional argument(s) for methods.
#' @details \code{pfamily} is a generic function with methods for classed \code{glmb} and 
#' \code{lmb}. Many \code{glmb} models currently only have implementations for the \code{dNormal()} 
#' prior family. The \code{Gamma()} family also works with the \code{dGamma()} prior 
#' family while the \code{gaussian()} family works with the \code{dGamma()}, 
#' \code{dNormal_Gamma()},and \code{dIndependent_Normal_Gamma()} families.  
#' @return An object of class \code{"pfamily"} (which has a concise print method). This is a
#' list with elements.
#' \item{pfamily}{character: the pfamily name}
#' \item{prior_list}{a list with the prior parameters associated with the prior specification}
#' \item{okfamilies}{currently implemented families for which the prior family can be used.}
#' \item{simfun}{function: the function used to generate samples from the posterior density. 
#' All currently implemented pfamiles have simulation functions that generate iid samples
#' for the associated posterior distribution.}
#' @example inst/examples/Ex_confint.glmb.R
#' @export 
#' @exportClass pfamily 
#' @rdname pfamily
#' @order 1

pfamily <- function(object, ...) UseMethod("pfamily")

#' @export 
#' @method pfamily default

pfamily.default <- function(object, ...){

  if(is.null(object$pfamily)) stop("no pfamily object found")
  if(!class(object$pfamily)=="pfamily") stop("Object named pfamily is not of class pfamily")
  return(object$pfamily)
}


#' @export 
#' @method print pfamily

#' @rdname pfamily
#' @order 6

print.pfamily <- function(x, ...)
{
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Prior Family:", x$pfamily, "\n\n")
  cat("Prior List:\n\n")
  print(x$prior_list)
  
  invisible(x)
}

#' @export 
#' @rdname pfamily
#' @order 2

dNormal<-function(mu,Sigma,dispersion=NULL){
  
  ## Check that the inputs are numeric
  
  if(is.numeric(mu)==FALSE||is.numeric(Sigma)==FALSE) stop("non-numeric argument to numeric function")

  mu=as.matrix(mu,ncol=1) ## Force mu to matrix
  Sigma=as.matrix(Sigma)  ## Force Sigma to matrix 
  
  nvar=length(mu)
  nvar1=nrow(Sigma)
  nvar2=ncol(Sigma)
  
  if(!nvar==nvar1||!nvar==nvar2) stop("dimensions of mu and Sigma are not consistent")

  ## Check for symmetry and positive definiteness
  if(!isSymmetric(Sigma))stop("matrix Sigma must be symmetric")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(Sigma, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  
  if(is.null(dispersion)) dispersion=1
  if(!is.null(dispersion)){
    if(!is.numeric(dispersion)) stop("non-numeric argument to numeric function")
    if(!length(dispersion)==1) stop("dispersion has length>1")
    if(!length(dispersion)>0) stop("dispersion must be >0")
  }
    
  okfamilies <- c("gaussian","poisson","binomial","quasipoisson","quasibinomial","Gamma")

  plinks<-function(family){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-c("log")		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-c("logit","probit","cloglog")		
    if(family$family=="Gamma") oklinks<-c("log")	
    return(oklinks)
  }
  
  prior_list=list(mu=mu,Sigma=Sigma,dispersion=dispersion)
  attr(prior_list,"Prior Type")="dNormal"  

  outlist=list(pfamily="dNormal",prior_list=prior_list,okfamilies=okfamilies,
  plinks=plinks,             
  simfun=rglmb_norm_reg)
  attr(outlist,"Prior Type")="dNormal"             
  class(outlist)="pfamily"
  outlist$call<-match.call()
  return(outlist)
  }

#' @export 
#' @rdname pfamily
#' @order 3

dGamma<-function(shape,rate,beta){

  if(is.numeric(shape)==FALSE||is.numeric(rate)==FALSE||is.numeric(beta)==FALSE) stop("non-numeric argument to numeric function")
  
  if(length(shape)>1) stop("shape is not of length 1")
  if(length(shape)>1) stop("rate is not of length 1")
  if(shape<=0) stop("shape must be>0")
  if(rate<=0) stop("rate must be>0")
  
  beta=as.matrix(beta,ncol=1)
  
  okfamilies <- c("gaussian","Gamma")
  
  plinks<-function(family){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-NULL		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-NULL		
    if(family$family=="Gamma") oklinks<-c("log")	
    return(oklinks)
  }
  
  prior_list=list(shape=shape,rate=rate,beta=beta)
  attr(prior_list,"Prior Type")="dGamma"  
  outlist=list(pfamily="dGamma",prior_list=prior_list,okfamilies=okfamilies,
               plinks=plinks,             
               simfun=rglmb_dispersion)
               
  attr(outlist,"Prior Type")="dGamma"
  class(outlist)="pfamily"
  outlist$call<-match.call()
  
  return(outlist)

}

#' @export 
#' @rdname pfamily
#' @order 4

dNormal_Gamma<-function(mu, Sigma,shape, rate){

  ############################################################  
  
  if(is.numeric(mu)==FALSE||is.numeric(Sigma)==FALSE) stop("non-numeric argument to numeric function")
  if(is.numeric(shape)==FALSE||is.numeric(rate)==FALSE) stop("non-numeric argument to numeric function")
  
  if(length(shape)>1) stop("shape is not of length 1")
  if(length(shape)>1) stop("rate is not of length 1")
  if(shape<=0) stop("shape must be>0")
  if(rate<=0) stop("rate must be>0")
  
  mu=as.matrix(mu,ncol=1) ## Force mu to matrix
  Sigma=as.matrix(Sigma)  ## Force Sigma to matrix 
    
  nvar=length(mu)
  nvar1=nrow(Sigma)
  nvar2=ncol(Sigma)
  
  if(!nvar==nvar1||!nvar==nvar2) stop("dimensions of mu and Sigma are not consistent")
  
  ## Check for symmetry and positive definiteness
  if(!isSymmetric(Sigma))stop("matrix Sigma must be symmetric")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(Sigma, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  
  
  ############################################################
  
  okfamilies <- c("gaussian") # Unclear if this could be used for Gamma  or quasi-families

  plinks<-function(family){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-NULL		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-NULL		
    if(family$family=="Gamma") oklinks<-NULL	
    return(oklinks)
  }
  
  prior_list=list(mu=mu,Sigma=Sigma,shape=shape,rate=rate)
  attr(prior_list,"Prior Type")="dNormal_Gamma"  
  outlist=list(pfamily="dNormal_Gamma",call=call,prior_list=prior_list,
    okfamilies=okfamilies,plinks=plinks,simfun=rnorm_gamma_reg)
  
  attr(outlist,"Prior Type")="dNormal_Gamma"             
  class(outlist)="pfamily"
  outlist$call<-match.call()
  
  return(outlist)
  }

