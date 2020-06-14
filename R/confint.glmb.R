#' Extract DIC from a Fitted Bayesian Model
#'
#' Computes the Deviance Information Criterion for a fitted Bayesian
#' Generalized Linear Model
#' @param object an object of class \code{"glmb"}, typically the result of a call to \link{glmb}
#' @param parm an object of class \code{"glmb"}, typically the result of a call to \link{glmb}
#' @param level an object of class \code{"glmb"}, typically the result of a call to \link{glmb}
#' @param ... further arguments to or from other methods
#' @return A list with the Estimated effective number of parameters \code{pD}
#' and the \code{DIC} from the object \code{fit} of class \code{"glmb"}. See \code{\link{glmbdic}}
#' for details on the definition of these objects.
#' @example inst/examples/Ex_extractAIC.glmb.R

confint.glmb<-function(object,parm,level=0.95,...)
{
  a <- (1 - level)/2
  a <- c(a, 1 - a)
#  pnames <- colnames(object$coefficients)
#  fac <- qt(a, object$df.residual)
#  ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm, 
#                                                                   pct))
#  ci[] <- cf[parm] + ses[parm] %o% fac
#  ci
#  if (missing(parm))       ## for now force parm to be pnames
#    parm <- pnames 
ci <-t(apply(object$coefficients, 2, FUN = quantile,probs=a))  

#dimnames(ci)=list(parm, a)
return(ci)
}
