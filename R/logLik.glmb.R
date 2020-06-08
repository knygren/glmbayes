#' Extract Log-Likelihood
#'
#' This function is a method function for the \code{glmb} class used to Extract the
#' Log-Likelihood from a Bayesian Generalized Linear Model.
#' @param object an object of class \code{glmb}, typically the result of a call to \link{glmb}
#' @param ... further arguments to or from other methods
#' @return The sum of \code{x} and \code{y}
#' @examples
#' 1+1
#' 10+1



logLik.glmb<-function(object,...){
  
  y<-object$glm$y
  x<-object$glm$x
  wt<-object$glm$prior.weights
  dispersion<-object$dispersion
  
  
  #alpha # to be updated 
  f1<-object$famfunc$f1
  n<-length(object$coefficients[,1])
  
  logLikout<-matrix(0,nrow=n,ncol=1)
  #f1temp(b=out$coefficients,y=out$y,x=out$x,alpha=0,wt=Claims/dispersion)
  for(i in 1:n){
    logLikout[i,1]<--f1(object$coefficients[i,],y=y,x=x,wt=wt/dispersion)
  }
  logLikout
}
