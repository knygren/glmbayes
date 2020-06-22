#' Extract Log-Likelihood
#'
#' This function is a method function for the \code{"glmb"} class used to 
#' Extract the log-Likelihood from a Bayesian Generalized Linear Model.
#' @param x an object of class \code{glmb}, typically the result of a call to \link{glmb}
#' @param ... further arguments to or from other methods
#' @return The function returns a vector, \code{logLikout} with the estimated log-likelihood for each draw. 
#' @example inst/examples/Ex_logLik.glmb.R
#' @export
#' @method formula summary.rglmb



formula.summary.rglmb<-function(x,...){
  
  z=x
  y=z$y
  x=z$x
  return(formula(glm(y~x-1,family=family(z))))
  
}