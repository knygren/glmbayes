#' Predict Method for Bayesian GLM Fits
#'
#' Obtains predictions and optionally estimates standard errors of those
#' predictions from a fitted Bayesian generalized linear model object.
#' 
#' @param object a fitted object of class inheriting from \code{"glmb"}.
#' @param newdata optionally, a data frame in which to look for variables
#' with which to predict. If omitted, the fitted linear predictors are used.
#' @param type the type of prediction required.  The default is on the scale of the linear predictors; 
#' the alternative \code{"response"} is on the scale of the response variable.  Thus for a default
#' binomial model the default predictions are of log-odds (probabilities
#' on logit scale) and \code{type = "response"} gives the predicted
#' probabilities.  The \code{"terms"} option returns a matrix giving the
#' fitted values of each term in the model formula on the linear predictor
#' scale (not implemented).
#' @param se.fit logical switch indicating if standard errors are required (not implemented).
#' @param dispersion the dispersion of the Bayesian GLM fit to be assumed in 
#' computing the standard errors. If omitted, that returned by the \code{summary} applied 
#' to the object is used.
#' @param terms with \code{type="terms"} by default all terms are returned.
#' A character vector specifies which terms are to be returned.
#' @param na.action function determining what should be done with missing values in
#' \code{newdata}. The default is to predict \code{NA}.
#' @param ... further arguments passed to or from other methods.
#' @return A list with the Estimated effective number of parameters \code{pD}
#' and the \code{DIC} from the object \code{fit} of class \code{"glmb"}. See \code{\link{glmbdic}}
#' for details on the definition of these objects.
#' @details If \code{newdata} is omitted the predictions are based on the data
#' used for the fit.  In that case how cases with missing values in the
#' original fit is determined by the \code{na.action} argument of that
#' fit.  If \code{na.action = na.omit} omitted cases will not appear in 
#' the residuals, whereas if \code{na.action = na.exclude} they will 
#' appear (in predictions and standard errors), with residual value
#' \code{NA}.  See also \code{\link{napredict}}.
#' @note Variables are first looked for in \code{newdata} and then searched for
#' in the usual way (which will include the environment of the formula 
#' used in the fit).  A warning will be given if the 
#' variables found are not of the same length as those in \code{newdata}
#' if it was supplied.
#' @example inst/examples/Ex_predict.glmb.R

predict.glmb<-function(object,newdata=NULL,type="link",
                       se.fit = FALSE, dispersion = NULL, terms = NULL, 
                       na.action = na.pass,...)
{

  ## Note: may need adjustment for case whee dispersion is random !
  ## residual.scalen needs to be adjusted every iteration!
  
  ## If newdata is missing, simply return the model predictions 
  ## on either the scale of the linear predictors or on the scale
  ## of the response variable (the former is the default)
   
  if (missing(newdata)) {
      if(type=="link")   pred <-  object$linear.predictors
      if(type=="response")  pred <-   object$fitted.values
      }
  
  else{
    
    ## if newdata is missing build matrix for storing predictions 
    ## draws come down, predictions for rows in newdata come across
    n=nrow(object$coefficients)
    nvars=ncol(object$coefficients)
    betas=object$coefficients
    npreds=nrow(newdata)
    pred=matrix(0,nrow=n,ncol=nrow(newdata))
    
    ## Initialize lm object to glm object from glmb function
    lm_object=object$glm ## set object to glm object from glmb class

    
    residual.scale <- as.vector(sqrt(dispersion))
    
    ## Loop through draws
    for(i in 1:n)
    {
      # Update object and then return predictions from predict.lm
      lm_object$coefficients=betas[i,1:nvars]
      pred[i,1:npreds]<- predict.lm(object, newdata, se.fit, 
                              scale = residual.scale, type =  "response",
                terms = terms, na.action = na.action)

      ## Rescale predictiosn if type="response"
        if(type=="response") pred[i,1:npreds]<- family(object$glm)$linkinv(pred[i,1:npreds])
         
    }
    
  }
  return(pred)
  
}
