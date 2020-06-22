#' Accessing Bayesian Generalized Linear Model Fits
#'
#' These functions are all \link{methods} for class \code{glmb} or \code{summary.glmb} objects.
#' @param object an object of class \code{glmb}, typically the result of a call to \link{glmb}
#' @param ysim Optional simulated data for the data y.
#' @param \ldots further arguments to or from other methods
#' @return A matrix \code{DevRes} of dimension \code{n} times \code{p} containing
#' the Deviance residuals for each draw. If ysim is provided, the residuals are based
#' on a comparison to the simulated data instead. The credible intervals
#' for residuals based on simulated data should be a more approproiate measure of
#' whether individual residuals represent outliers or not  
#' @example inst/examples/Ex_residuals.glmb.R
#' @export
#' @method residuals glmb

residuals.glmb<-function(object,ysim=NULL,...)
{
  y<-object$glm$y	
  n<-length(object$coefficients[,1])
  
  ## Updated to use prior.weights - likely matters for binomial data
  ## Need to verify this performs as expected
  wts <- object$glm$prior.weights

    fitted.values<-object$fitted.values
  dev.residuals<-object$glm$family$dev.resids
  DevRes<-matrix(0,nrow=n,ncol=length(y))

  
  for(i in 1:n)
  {
    if(is.null(ysim))    DevRes[i,]<-sign(y-fitted.values[i,])*sqrt(dev.residuals(y,fitted.values[i,],wts))
    else(DevRes[i,]<-sign(ysim[i,]-fitted.values[i,])*sqrt(dev.residuals(ysim[i,],fitted.values[i,],wts)))
  }
  
  colnames(DevRes)<-names(y)
  DevRes
}
