#' Extract the Number of Observations from a Fit.
#'
#' 
#' Extract the Number of \sQuote{observations} from a Fit.
#' @param object A fitted model object.
#' @param \ldots further arguments to be passed to other methods
#' @return A single number, normally an integer.  Could be \code{NA}.
#' @example inst/examples/Ex_extractAIC.glmb.R

nobs.glmb<-function(object,...)
{
  return(nobs(object$glm))
  
}
