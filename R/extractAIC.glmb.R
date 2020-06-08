#' Extract DIC from a Fitted Bayesian Model
#'
#' Computes the Deviance Information Criterion for a fitted Bayesian
#' Generalized Linear Model
#' @param fit an object of class \code{glmb}, typically the result of a call to \link{glmb}
#' @param ... further arguments to or from other methods
#' @return The sum of \code{x} and \code{y}
#' @examples
#' 1+1
#' 10+1

extractAIC.glmb<-function(fit,...)
{
  c(fit$pD,fit$DIC)
  
}
