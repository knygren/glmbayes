#' Calculate Variance-Covariance Matrux for a Fitted Model Object
#'
#' Returns the variance-covariance matrix of the main parameters of
#' a fitted model object.
#' @param object fitted model object, typically the result of a call to \code{"glmb"}.
#' @param \ldots additional arguments for method functions. 
#' @return A matrix of estimated covariances between the parameter estimates
#' in the linear or non-linear predictor of the model. This should have
#' row and column names corresponding to the parameter names given by the
#' \code{\link{coef}} method.
#' @example inst/examples/Ex_extractAIC.glmb.R

vcov.glmb<-function(object,...)
{
  return(cov(object$coefficients))
  
}
