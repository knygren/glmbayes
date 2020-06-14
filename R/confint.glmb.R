#' Credible Intervals for Model Parameters
#'
#' Computes credible intervals for Model Parameters
#' @param object a fitted model object of class \code{"glmb"}. Typically the result of a call to \link{glmb}.
#' @param parm a specification (not yet implmented) of which parameters are to be given credible sets,
#' either a vector of numbers or a vector of names. If missing, all parameters are considered. 
#' @param level the credible interval required.
#' @param \ldots additional argument(s) for methods.
#' @return A matrix (or vector) with columnes giving lower and
#' upper credible limits for each parameter. These will be labeled
#' (1-level)/2 and 1-(1-level)/2 in \% (by default 2.5\% and 97.5\%).
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
