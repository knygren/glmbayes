#' Fitter Function for Bayesian Generalized Linear Models
#'
#' Basic computing engine called to replicate the output from the lm.fit function for an already
#' optimized Bayesian Generalized Linear model.
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



