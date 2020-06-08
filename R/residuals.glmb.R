#' Accessing Bayesian Generalized Linear Model Fits
#'
#' These functions are all \link{methods} for class \code{glmb} or \code{summary.glmb} objects.
#' @param object an object of class \code{glmb}, typically the result of a call to \link{glmb}
#' @param ... further arguments to or from other methods
#' @return The sum of \code{x} and \code{y}
#' @examples
#' 1+1
#' 10+1

residuals.glmb<-function(object,...)
{
  y<-object$glm$y	
  n<-length(object$coefficients[,1])
  fitted.values<-object$fitted.values
  
  dev.residuals<-object$glm$family$dev.resids
  DevRes<-matrix(0,nrow=n,ncol=length(y))
  #	yrep<-matrix(0,nrow=n,ncol=length(y))
  
  
  for(i in 1:n)
  {
    DevRes[i,]<-sign(y-fitted.values[i,])*sqrt(dev.residuals(y,fitted.values[i,],1))
    #	if(fit$family$family=="poisson")yrep[i,]<-rpois(length(y),fitted.values[i,])
    
  }
  
  colnames(DevRes)<-names(y)
  DevRes
}
