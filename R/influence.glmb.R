#' Bayesian Regression Diagnostics
#'
#' This function provides the basic quantities which are used in forming a wide variety of diagnostics for checking
#' the quality of Bayesian regression fits.
#' @inheritParams stats::lm.influence
#' @return a \code{\link{list}} wih components:
#' @example inst/examples/Ex_glmb.wfit.R
#' @method influence glmb
#' @export 



influence.glmb<-function(model,do.coef=TRUE){
  
  # Just tell function to use fit component
  # necessary because coefficients are draws and not modes
  # so not all items returned by fitting function can be include in glmb returned list
    
  return(influence(model$fit,do.coef)) 
  
}


#' Bayesian Regression Diagnostics
#'
#' This function provides the basic quantities which are used in forming a wide variety of diagnostics for checking
#' the quality of Bayesian regression fits.
#' @inheritParams stats::influence.measures
#' @return a \code{\link{list}} wih components:
#' @example inst/examples/Ex_glmb.wfit.R
#' @method influence.measures glmb
#' @export 


influence.measures.glmb<-function (model, infl = influence(model,do.coef=FALSE)) 
{

  if(is.null(infl)) infl=influence(model,do.coef=FALSE)
  
  return(influence.measures(model$fit,infl))
}