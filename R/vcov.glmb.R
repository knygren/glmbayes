#' Extract DIC from a Fitted Bayesian Model
#'
#' Computes the Deviance Information Criterion for a fitted Bayesian
#' Generalized Linear Model
#' @param object an object of class \code{"glmb"}, typically the result of a call to \link{glmb}
#' @param ... further arguments to or from other methods
#' @return A list with the Estimated effective number of parameters \code{pD}
#' and the \code{DIC} from the object \code{fit} of class \code{"glmb"}. See \code{\link{glmbdic}}
#' for details on the definition of these objects.
#' @example inst/examples/Ex_extractAIC.glmb.R

vcov.glmb<-function(object,...)
{
  return(cov(object$coefficients))
  
}
