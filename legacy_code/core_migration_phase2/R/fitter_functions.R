#' Bayesian Weighted Fitting Engines
#'
#' These functions provide the Bayesian analogue of \code{lm.wfit}. They implement the core weighted least squares
#'step used inside Bayesian linear linear models, incorporating prior
#' precision and posterior mode information.
#'
#' \code{rNormal_reg.wfit} performs the Bayesian weighted least squares update
#' for linear models under a Normal prior.
#'
#' \code{glmb.wfit} performs the corresponding update for generalized linear
#' models, reconstructing the weighted least squares step using the posterior
#' mode and the GLM family functions.
#' 
#' @rdname bayes_wfit
#' @param mu Prior mean vector of length \code{p}.
#' @param P Prior precision matrix of dimension \code{p * p}.
#' @inheritParams stats::lm.wfit
#' @return a \code{\link{list}} with components:
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
  
  ## Do cholesky decomposition - Might be inaccurate/inefficient
  
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
  
  #                      E'E
  lmf$S=crossprod(Z-W%*%lmf$Btilde)  # This part likely is used for the Gamma or Wishart simulation - Essentially RSS
  
  # Assign the "lm" class for now to allow the method functions for influence measures to work
  
  class(lmf)="lm"
  
  
  return(lmf)
  
  
}



#' @param Bbar Prior mean vector of length \code{p}.
#' @param P Prior precision matrix of dimension \code{p * p}.
#' @param betastar Posterior mode vector of length \code{p} which has already been estimated.
#' @param weights an optional vector of \emph{prior weights} to be used in the fitting process. 
#' Should be \code{NULL} or a numeric vector.
#' @param family a description of the error distribution and link function to be used in the model.
#' Should be a family function. (see \code{\link{family}} for details of family functions.)
#' @inheritParams stats::lm.wfit
#' @return a \code{\link{list}} wih components:
#' @example inst/examples/Ex_glmb.wfit.R
#' @rdname bayes_wfit
#' @export


glmb.wfit<-function(x,y,weights=rep.int(1, nobs),offset=rep.int(0, nobs),family=gaussian(),Bbar,P,betastar,method="qr",tol=1e-7,singular.ok=TRUE,...){
  
  # Basic checks like in glm.fit
  
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if(is.matrix(y)) rownames(y) else names(y)
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  ## define weights and offset if needed
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  
  # Get needed family functions
  
  
  linkinv<-family$linkinv
  dev.resids<-family$dev.resids
  mu.eta<-family$mu.eta
  variance=family$variance
  
  # Get Cholesky decomposition for Prior Precision
  
  RA=chol(P)
  
  # Update the constants (like in glm.fit)
  
  start <- betastar
  eta <- drop(x %*% start)   # This should be fitted values using posterior mode
  mu <- linkinv(eta <- eta + offset)
  dev <- sum(dev.resids(y, mu, weights))
  
  good <- weights > 0
  varmu <- variance(mu)[good]  # For poisson, this is just lambda =exp(Xbetastar)=mu !
  
  
  if (anyNA(varmu))
    stop("NAs in V(mu)")
  if (any(varmu == 0))
    stop("0s in V(mu)")
  mu.eta.val <- mu.eta(eta)
  if (any(is.na(mu.eta.val[good])))
    stop("NAs in d(mu)/d(eta)")
  ## drop observations for which w will be zero
  good <- (weights > 0) & (mu.eta.val != 0)
  
  if (all(!good)) {
    conv <- FALSE
    warning(gettextf("no observations informative"))
  }
  
  # Update z and w as in glm.fit (after call to fitting function)
  
  z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
  w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])  # These are essentially weights
  
  # Bind values used by glm.fit with the prior components
  
  W2=rbind(x[good, , drop = FALSE] * w,RA)    # W2 should be modified design matrix !
  Z2=rbind(matrix(z * w,ncol=1),matrix(RA%*%Bbar,ncol=1)) ## Z2 Should be the modified y vector!
  
  fit<-lm.fit(W2,Z2)
  
  class(fit)="lm"
  
  ## For now print both for comparison purposes 
  ## Should add comparison measure to see how close these are
  ## Add return a warning if not close [They should essentially match]
  
  ## print(fit$coefficients)
  ##  print(betastar)
  
  return(fit)
}


