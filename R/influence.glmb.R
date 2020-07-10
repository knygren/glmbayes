#' Bayesian Regression Diagnostics
#'
#' This function provides the basic quantities which are used in forming a wide variety of diagnostics for checking
#' the quality of Bayesian regression fits.
#' @inheritParams stats::lm.influence
#' @return a \code{\link{list}} wih components:
#' @example inst/examples/Ex_glmb.wfit.R
#' @method influence glmb
#' @export 


influence.glmb<-function(model,...){
  
  # Just tell function to use fit component
  # necessary because coefficients are draws and not modes
  # so not all items returned by fitting function can be include in glmb returned list
  
  return(influence(model$fit,...)) 
  
}


#' Bayesian Regression Diagnostics
#'
#' This function provides the basic quantities which are used in forming a wide variety of diagnostics for checking
#' the quality of Bayesian regression fits.
#' @param infl influence structure as returned by \code{influence.glmb} 
#' (the latter only for the glm method of \code{rstudent} and \code{cooks.distance}).
#' @inheritParams stats::influence.measures
#' @return a \code{\link{list}} wih components:
#' @example inst/examples/Ex_glmb.wfit.R
#' @export 


glmb.influence.measures<-function (model, infl = influence(model)) 
{

  if(is.null(infl)) infl=influence(model)
  
  return(influence.measures(model$fit,infl))
}

#' @export 
#' @method rstandard glmb
#' @rdname glmb.influence.measures 

rstandard.glmb<-function(model,...,infl=influence(model)){

  if(is.null(infl)) infl=influence(model)
  
  return(rstandard(model$fit,infl))

}
  