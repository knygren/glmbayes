#' Model Deviance
#'
#' Returns the deviance of a fitted Bayesian Generalized Linear Model
#' @param object an object of class \code{"rglmb"}, typically the result of a call to \link{rglmb}
#' @param ... further arguments to or from other methods
#' @return A vector with the deviance extracted from the \code{object}.
#' @example inst/examples/Ex_extractAIC.glmb.R
#' @export
#' @method deviance rglmb

deviance.rglmb<-function(object,...)
{
   object2=summary(object)
   return(deviance(object2,...))

}
