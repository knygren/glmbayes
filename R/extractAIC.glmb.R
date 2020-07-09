#' Extract DIC from a Fitted Bayesian Model
#'
#' Computes the Deviance Information Criterion for a fitted Bayesian
#' Generalized Linear Model
#' @param fit an object of class \code{"glmb"}, typically the result of a call to \link{glmb}
#' @param ... further arguments to or from other methods
#' @return A list with the Estimated effective number of parameters \code{pD}
#' and the \code{DIC} from the object \code{fit} of class \code{"glmb"}. 
#' @example inst/examples/Ex_extractAIC.glmb.R
#' @export
#' @method extractAIC glmb

extractAIC.glmb<-function(fit,...)
{
  c(pD=fit$pD,DIC=fit$DIC)
  
}


## Alias for extractAIC function

#' @export
#' @rdname extractAIC.glmb
#'
extractDIC<-function(fit,...) UseMethod(extractAIC)
