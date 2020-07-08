#' Fitter Function for Bayesian Linear Models
#'
#' Basic computing engine called to find the posterior mode and a UL decomposition
#' @param mu Prior mean vector of length \code{p}.
#' @param P Prior precision matrix of dimension \code{p * p}.
#' @inheritParams stats::lm.wfit
#' @return a \code{\link{list}} wih components:
#' @example inst/examples/Ex_confint.glmb.R
#' @export 

rNormal_reg.wfit<-function(x,y,P,mu, w,offset=NULL,method="qr",tol=1e-7,singular.ok=TRUE,...){

  ## Handle all four cases of offset and wt here 
  ## Should determine where validity of input checking should be done.
  ## Seems like the lm.wfit function itself has some of this, so not all of these checks are necessary
  
  if(!is.null(offset)&!is.null(w)){  
    Y=matrix((y-offset)*sqrt(w),nrow=length(y))
    X=x*sqrt(w)  }

  if(!is.null(offset)&is.null(w)){  
    Y=matrix((y-offset),nrow=length(y))
    X=x  }
  
  if(is.null(offset)&!is.null(w)){  
    Y=matrix(y*sqrt(w),nrow=length(y))
    X=x*sqrt(w)  }
  
  if(is.null(offset)&is.null(w)){  
    Y=matrix(y,nrow=length(y))
    X=x}
  
  ## For now rename these (modify below to avoid)
  Bbar=mu
  A=P
  
  ### Do dimension checks (may want to do outside of this function)
  
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
  
  RA=chol(A)
  
  # Create modifed design matrix x and modifed observed y matrix
  
  W=rbind(X,RA)    # W should be modified design matrix !
  Z=rbind(Ytemp,matrix(RA%*%Bbar,ncol=1)) ## Z Should be the modified y vector!
  
  ## Call lm.fit
  
  lmf=lm.fit (W, Z,    offset = NULL, method = "qr", tol = 1e-7,
                     singular.ok = TRUE)
  
  ## Also do the IR Decomposition needed by the posterior simulation
  ## might be able to reuse qr output (TBD)
  
  #   note:  Y,X,A,Bbar must be matrices!
  lmf$IR=backsolve(chol(crossprod(W)),diag(k))
  #                      W'W = R'R  &  (W'W)^-1 = IRIR'  -- this is the UL decomp!
  
  lmf$k=k
  
  lmf$Btilde=matrix(lmf$coefficients,ncol=1)
  
  return(lmf)
  
  
}