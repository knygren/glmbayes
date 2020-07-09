#' Extract Log-Likelihood
#'
#' This function is a method function for the \code{"glmb"} class used to 
#' Extract the log-Likelihood from a Bayesian Generalized Linear Model.
#' @param object an object of class \code{glmb}, typically the result of a call to \link{glmb}
#' @param ... further arguments to or from other methods
#' @return The function returns a vector, \code{logLikout} with the estimated log-likelihood for each draw. 
#' @example inst/examples/Ex_logLik.glmb.R
#' @export
#' @method logLik glmb


logLik.glmb<-function(object,...){
  
  y<-object$y
  x<-object$x
  wt<-object$prior.weights
  dispersion<-object$dispersion
  
  
  #alpha # to be updated 
  f1<-object$famfunc$f1
  n<-length(object$coefficients[,1])
  
  logLikout<-matrix(0,nrow=n,ncol=1)

  if(length(dispersion)==1){  
    for(i in 1:n){
    logLikout[i,1]<--f1(object$coefficients[i,],y=y,x=x,wt=wt/dispersion)
    }
  }
  
  if(length(dispersion)>1){  
    for(i in 1:n){
      logLikout[i,1]<--f1(object$coefficients[i,],y=y,x=x,wt=wt/dispersion[i])
    }
  }
  
  logLikout
}
