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
#' @param olddata a data frame that should contain all the variables used in the 
#' original model specification. Must currently be provided whenever newdata is 
#' provided. Both olddata and newdata are subsetted to the model variables extracte
#' from the object model formula and rbind(olddata,newdata) must be valid after this step. 
#' A check is also run to verify if the resulting x matrix from olddata is consistent with that 
#' from the original model object.
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
#' @export
#' @method predict glmb

predict.glmb<-function(object,newdata=NULL,type="link",
                       se.fit = FALSE, dispersion = NULL, terms = NULL, 
                       na.action = na.pass,olddata=NULL,...)
{

  ## Additional changes to consider
  
  ### 1) Replacing last step with a step that passes individual draws to glm [might be safer]
  ### 2) QC that offsets are handled properly
  ### 3) Implement version that allows y to be generated
  
  

  ## Note: may need adjustment for case when dispersion is random 
  ## residual.scale be may not be 1 or link-inverse may not be accurate
  ## Possible that this only matters if we want to simulate
  ## from likelihood function, so may be ok. Need to verify.
  
  ## If newdata is missing, simply return the model predictions 
  ## on either the scale of the linear predictors or on the scale
  ## of the response variable (the former is the default)
   
  if (missing(newdata)) {
      if(type=="link")   pred <-  object$linear.predictors
      if(type=="response"){
                          pred <-   object$fitted.values
                          }
      return(pred)
         
            }
  
  else{

    ## For now requires olddata be provided so that data can be compared for consistency
        
    
    ## Pull model frame from original model and olddata
    ## This can likely be remove and moved inside complete_newdata function

    original_frame=object$glm$model
    
    # This could fail so use try
    
    newdata_ok=1
    newdata_frame=try(model.frame(formula(object$glm$model),newdata),silent=TRUE)
    
    #if(class(newdata_frame)=="data.frame") print("newdata built a model frame")    
    if(class(newdata_frame)=="try-error") newdata_ok=0    
    

    ## Attempt comparison of newdata_frame to original_frame
    ## This should pass if all variables (including dependent variables) are contained in 
    ##  newdata_frame and the attributes match

    ## If newdata_frame matches original_frame other than in terms of number of rows, 
    ## produce new_matrix right away (this should be correct)

    if(newdata_ok==1){
        if(isTRUE(Compare_Model_Frames(original_frame,newdata_frame,Check_Rows = FALSE)))
    {
      new_x=model.matrix(formula(object$glm$model),newdata_frame)

      pred=generate_predictions(object,type=type,new_x,newdata_frame)  
      return(pred)
      
      
    }
    else{newdata_ok=0}
      
    }  
    
    ## if newdata_frame does not match, try constructing a valid frame using
    ## information from olddata
    
    if (missing(olddata)){ stop("no olddata available to complete newdata")}
    
    x_matrices=complete_newdata(object,newdata,olddata,type)
    
    new_x=x_matrices$x_new
    new_mod_frame=x_matrices$mod_frame_new
    
    pred=generate_predictions(object,type=type,new_x,new_mod_frame)  
    return(pred)
    
  }
  
}
