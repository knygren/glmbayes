#' The Bayesian Non-Gaussian Generalized Linear Model Distribution
#'
#' \code{rnnorm_reg} is used to generate iid samples from Non-Gaussian Bayesian Generalized Linear Models with a normal prior. 
#' The model is specified by providing a data vector, a design matrix, and 2 prior constants.
#' @aliases
#' rnnorm_reg
#' rnnorm_reg_cpp
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x a design matrix of dimension \code{m * p}.
#' @param mu a vector of length \code{p} giving the prior means of the variables in the design matrix.
#' @param P a positive-definite symmetric matrix of dimension \code{p * p} specifying the prior precision matrix of the variable.
#' @param wt an optional vector of \sQuote{prior weights} to be used in the fitting process. Should be NULL or a numeric vector.
#' @param dispersion the dispersion parameter. Either a single numerical value or NULL (the default). Must be provided here, use \code{\link{rnorm_gamma_reg}} to give the dispersion a prior.
#' @param nu Prior shape parameter for the dispersion parameter (gaussian model only).
#' @param V Prior rate parameter for the dispersion parameter (gaussian model only).
#' @param family a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)
#' @param offset2 this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.
#' @param start an optional argument providing starting values for the posterior mode optimization.
#' @param Gridtype an optional argument specifying the method used to determine number of likelihood subgradient densities used to construct the enveloping function.
#' @return \code{rnnorm_reg} returns a object of class \code{"rglmb"}.  The function \code{summary} 
#' (i.e., \code{\link{summary.rglmb}}) can be used to obtain or print a summary of the results.
#' The generic accessor functions \code{\link{coefficients}}, \code{\link{fitted.values}},
#' \code{\link{residuals}}, and \code{\link{extractAIC}} can be used to extract
#' various useful features of the value returned by \code{\link{rglmb}}.
#' An object of class \code{"rglmb"} is a list containing at least the following components:
#' \item{coefficients}{a \code{n} by \code{length(mu)} matrix with one sample in each row}
#' \item{PostMode}{a vector of \code{length(mu)} with the estimated posterior mode coefficients}
#' \item{Prior}{A list with two components. The first being the prior mean vector and the second the prior precision matrix}
#' \item{iters}{an \code{n} by \code{1} matrix giving the number of candidates generated before acceptance for each sample.}
#' \item{famfunc}{an object of class \code{"famfunc"}}
#' \item{Envelope}{an object of class \code{"envelope"}  }
#' \item{dispersion}{the dispersion parameter used in the model}
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
#' For the gaussian family, iid samples from the posterior density are genererated using 
#' standard simulation procedures for multivariate normal densities. For all other 
#' families, the samples are generated using accept-reject procedures by leveraging the 
#' likelihood subgradient approach of Nygren and Nygren (2006). This approach relies on
#' tight enveloping functions that bound the posterior density from above. The Gridtype 
#' parameter is used to select the method used for determining the size of a Grid used 
#' to build the enveloping function. See \code{\link{EnvelopeBuild_c}} for details. 
#' Depending on the selection, the time to build the envelope and the acceptance rate 
#' during the simulation process may vary. The returned value \code{iters} contains the 
#' number of candidates generated before acceptance for each draw.
#' @example inst/examples/Ex_rnnorm_reg.R


rnnorm_reg<-function(n=1,y,x,mu,P,wt=1,dispersion=NULL,nu=NULL,V=NULL,family=gaussian(),offset2=NULL,start=NULL,Gridtype=3)
  {
  
  if(is.numeric(n)==FALSE||is.numeric(y)==FALSE||is.numeric(x)==FALSE||
     is.numeric(mu)==FALSE||is.numeric(P)==FALSE) stop("non-numeric argument to numeric function")
  
  x <- as.matrix(x)
  mu<-as.matrix(as.vector(mu))
  P<-as.matrix(P)    
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) 
    rownames(y)
  else names(y)
  if(length(n)>1) n<-length(n)	   
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(offset2)) 
    offset2 <- rep(0, nobs)
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
  
  if (is.null(start)) 
    start <- mu
  if (is.null(offset2)) 
    offset2 <- rep.int(0, nobs)
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  okfamilies <- c("poisson","binomial","quasipoisson","quasibinomial","Gamma")
  if(family$family %in% okfamilies){
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
  
  
  if(family$family=="poisson"||family$family=="binomial")dispersion2<-1
  else dispersion2<-dispersion
  

    if(is.null(dispersion)){dispersion2=1}

      # f1, f2, and f3 passed here - Likely legacy of R code
    ## Can eliminate and replace with calling of corresponding c++ functions
    outlist<-rnnorm_reg_cpp(n=n,y=y,x=x,mu=mu,P=P,offset2=offset2,wt=wt,dispersion=dispersion2,
                                famfunc=famfunc,f1=f1,f2=f2,f3=f3,
                                start=start,family=family$family,link=family$link,Gridtype=Gridtype)
  
  colnames(outlist$coefficients)<-colnames(x)
  
  outlist$call<-match.call()
  
  class(outlist)<-c(outlist$class,"rglmb")
  outlist
  
}

