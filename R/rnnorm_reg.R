#' The Bayesian Non-Gaussian Generalized Linear Model Distribution
#'
#' \code{rnnorm_reg} is used to generate iid samples from Non-Gaussian Bayesian Generalized Linear Models with a normal prior. 
#' The model is specified by providing a data vector, a design matrix, and 2 prior constants.
#' @aliases
#' rnnorm_reg
#' rnnorm_reg_cpp
#' @param n An integer
#' @param y A vector
#' @param x A matrix
#' @param mu A vector
#' @param P A matrix
#' @param wt A vector or a numeric constant
#' @param dispersion A vector or a numeric constant
#' @param nu An optional numeric constant
#' @param V An optional numeric constant
#' @param family A family
#' @param offset2 A vector or a numeric constant
#' @param start A vector
#' @param Gridtype An integer
#' @return The sum of \code{x} and \code{y}
#' @examples
#' 1+1
#' 10+1


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

