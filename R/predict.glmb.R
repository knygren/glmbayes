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
#' @return A matrix of predictions where the rows correspond to the draws from the estimated model, and
#' the columns to the observations in the newdata dataset (or the original data if newdata is missing). 
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

    original_frame=object$model
    
    # This could fail so use try
    
    newdata_ok=1
    newdata_frame=try(model.frame(formula(object$model),newdata),silent=TRUE)
    
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
      new_x=model.matrix(formula(object$model),newdata_frame)

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



#' Compare Model Frames
#'
#' Compares model frames and stops if the dimensions, names, or attributes don't match 
#' @param oldframe The target model frame.
#' @param newframe The model frame being compared to the target
#' @param Check_Rows Logical indicating whether Number of Rows Should be checked (the default). If FALSE,
#' differing number of rows are allowed
#' @return 0 if the two model frames pass the comparison. If not, the function stops.
#' @example inst/examples/Ex_confint.glmb.R
#' @noRd

Compare_Model_Frames<-function(oldframe,newframe,Check_Rows=TRUE){
  
  # Check that dimensions match
  if(Check_Rows==TRUE)  if(!nrow(newframe)==nrow(oldframe)) {
    print("New Model Frame Does Not Have The Same Number of Rows")
    return(FALSE)
  }
  if(!ncol(newframe)==ncol(oldframe))  {
    print("New Model Frame Does Not Have The Same Number of Columns")
    return(FALSE)
  }
  # Check that Column Names Match
  
  for(i in 1:length(names(oldframe))) {
    if(!names(newframe)[i]==names(oldframe)[i])
    {
      print("New Model Frame Has Different Column Names")  
      return(FALSE)
    }
  }
  
  # Check that other attributes match
  
  for(i in 1:length(names(oldframe)))
  {
    if(!is.null(attr.all.equal(newframe[, names(newframe)[i]], oldframe[, names(oldframe)[i]])))
    { 
      not_equal=attr.all.equal(newframe[, names(newframe)[i]], oldframe[, names(oldframe)[i]])
      if(isTRUE(grepl("Attributes", not_equal))) 
      {
        print("Attributes are Different.")
        print(not_equal)
        return(FALSE)
      }
      
    }
  }
  return(TRUE)
}


#' Complete a newdata data frame
#'
#' Takes an incomplete newdata data frame missing some model variables
#' and attempts to complete it. 
#' @param object an object of class \code{glmb}
#' @param newdata A data frame containing newdata that is to be completed (when possible)
#' @param olddata A data frame that should contain the model variables (and values) from the original model
#' @param type Type of response. Either "link" or "response"
#' @return A design matrix x_new and a model frame mod_frame_new
#' @details Attempts to complete a newdata data frame with missing model
#' variables. This will typically fail if variables used in the
#' formula to derive exogenous or offset variables are missing 
#' (unless the environment contains the needed information). It will generally 
#' succeed if the all variables used to derive exogenous or offset
#' variables are present and olddata contains data for 
#' the missing variables that can be used to substitute the missing 
#' data in newdata with means from the olddata  
#' @example inst/examples/Ex_confint.glmb.R
#' @noRd


complete_newdata<-function(object,newdata,olddata,type){
  
  ## Check if olddata available
  if (missing(olddata)){ stop("no olddata available to complete newdata")}
  
  ## Pull original and "olddata" model frames  
  original_frame=object$model
  olddata_frame=model.frame(formula(object$model),olddata)
  
  ## Validate olddata model frame
  if(isFALSE(Compare_Model_Frames(original_frame,olddata_frame))){
    stop("Invalid olddata - generated model frame does not match original")}
  
  
  #############################################################################33
  ## Short internal function to check if model variables are present
  ## in sub-components of formular strings
  
  get_vars=function(all_vars,xstring){
    
    ## Initialize model 
    xvars=all_vars
    xflags=list(rep("TRUE",length(xvars)))
    
    check_data=data.frame(xvars=xvars,xflags=xflags)
    colnames(check_data)=c("xvars","xflags")
    
    for(i in 1:length(xvars)){
      if(isFALSE(grepl(xvars[i], xstring))) check_data$xflags[i]=FALSE
    }
    
    xdata=subset(check_data, xflags==TRUE)
    return(unique(xdata$xvars))
    
    ## Closing brackets for internal function
    
  }
  
  ###################################################################################3333
  
  
  # Step 1: get mod_vars and separate between those used for independent variables/offset calculations
  ## and those only used in calculation of dependent variable
  
  mod_vars=all.vars(formula(object))
  mod_formula=formula(object)
  
  # get strings holding information related to dependent and independent variables 
  
  ystring=as.character(mod_formula[2]) ## Seems to be string holding dependent variables
  xstring=as.character(mod_formula[3]) ## Seems to be string holding other variables
  
  ## Get subset of mod_vars contained in xstring 
  
  if(length(xstring)>0 & length(mod_vars)>0){
    x_vars=get_vars(all_vars=mod_vars,xstring=xstring)
  } 
  else
  {x_vars=NULL}
  
  # Get list of model variables that are not part of xstring
  
  y_vars1=setdiff(mod_vars,x_vars)
  
  ## get list of variables that are confirmed parts of ystring
  
  if(length(ystring)>0 & length(y_vars1)>0){
    y_vars2=get_vars(all_vars=y_vars1,xstring=ystring)
  }
  else{y_vars2=y_vars1}
  
  # Issue warning is no variables appear to be used to construct dependent variables
  
  if(length(y_vars1)==0) {warning("No data variables appear to be used in formula string
                                  for dependent variable. The number of observation set by the formula are likely fixed and recover of x 
                                  matrix will likely fail")}
  
  ## Check for any unclassified and added back to x_vars if the exist 
  ## (in theory there should be none)
  
  miss_vars=setdiff(y_vars1,y_vars2)
  
  if(length(miss_vars)>0){  
    warning("Some Model Variables Not Found - May cause unpredictable behavior")
    x_vars=c(x_vars,miss_vars)
  }
  ## List of y_vars and x_vars to use
  
  y_vars=y_vars2
  x_vars=x_vars
  
  ###############################################################################
  
  ## Check if subsetting olddata and newdata to mod_vars and x_vars respectively would return
  ##  errors
  
  tryCatch(olddata[mod_vars],error=function(e) {stop("olddata does not contain all model variables")})
  tryCatch(newdata[x_vars],error=function(e) {stop("newdata does not contain all explanatory variables")})
  
  ## Get model variables missing from newdata
  
  newvars=names(newdata)
  miss_newvars=setdiff(mod_vars, newvars)  # These should now be any missing dependent variables 
  
  ########################################################################
  
  # Find olddata means for missing variables and use to populate missing
  # variables for newdata
  
  if(length(miss_newvars)>0){
    
    ## Subset olddata to variables missing in newdata and get column means
    
    old_miss_newvars=olddata[miss_newvars]
    m_miss_newvars=round(colMeans(old_miss_newvars))
    
    for(i in 1:length(names(m_miss_newvars)))    {
      Temp1=matrix(round(colMeans(olddata[miss_newvars])[i]),nrow=nrow(newdata),ncol=1)  
      colnames(Temp1)=names(m_miss_newvars)[i]
      if(i==1) Temp2=Temp1
      else(Temp2=cbind(Temp2,Temp1))
    }
    newdata[miss_newvars] <- Temp2
  }
  
  #############################################
  
  ## Verify that newdata now contains all model variables
  
  tryCatch(newdata[mod_vars],error=function(e) {stop("Modified newdata does not contain all model variables")})
  
  ## Subset both olddata and newdata to mod_vars (can check)
  olddata=olddata[mod_vars]
  newdata=newdata[mod_vars]
  
  ########################################################################################
  
  ## Verify that olddata yields a design matrix consistent with original
  
  olddata_mod_frame=model.frame(mod_formula,olddata)
  olddata_mod_matrix=model.matrix(mod_formula,olddata_mod_frame)
  x_old=object$x
  x_new=olddata_mod_matrix
  
  if(isTRUE(all.equal(x_new,x_old))==FALSE) stop("olddata does not yield an x matrix consistent with that 
                                               stored in original model object")
  
  ######################################################################
  
  # get x and model frame matrices and then generate predictions
  
  x_matrices=get_x_matrix(object,olddata,newdata)
  
  return(x_matrices)
  
  
}


#' Design matrix and model frame for newdata
#'
#' Derives the model frame and design matrix for a complete version of the 
#' newdata data frame passed to the predict.glmb function
#' @param object an object of class \code{glmb}
#' @param olddata A data frame containing the model variables (and values) from the original model
#' @param newdata A data frame containing the model variables wih possible new values
#' @return A design matrix x_new and a model frame mod_frame_new
#' @example inst/examples/Ex_confint.glmb.R
#' @noRd

get_x_matrix<-function(object,olddata,newdata){
  
  nrow1=nrow(olddata)
  nrow2=nrow(newdata)
  
  ## Stack olddata and newdata together
  
  combodata=rbind(olddata,newdata)
  combo_mod_frame=model.frame(formula(object$formula),combodata)
  combo_mod_matrix=model.matrix(formula(object$formula),combo_mod_frame)
  
  x_old=combo_mod_matrix[1:(nrow1),]
  x_new=combo_mod_matrix[(nrow1+1):(nrow1+nrow2),]
  
  ## assign attributes of combo matrix to sub-matrices 
  
  attr(x_old,which="assign")<-attr(combo_mod_matrix,which="assign")
  attr(x_old,which="contrasts")<-attr(combo_mod_matrix,which="contrasts")
  
  attr(x_new,which="assign")<-attr(combo_mod_matrix,which="assign")
  attr(x_new,which="contrasts")<-attr(combo_mod_matrix,which="contrasts")
  
  mod_frame_old=combo_mod_frame[1:(nrow1),]
  mod_frame_new=combo_mod_frame[(nrow1+1):(nrow1+nrow2),]
  
  x=object$x
  
  if(!dim(x_old)[1]==dim(x)[1]) stop("Number of rows in final constucted x matrix for olddata does not match original")
  if(!dim(x_old)[2]==dim(x)[2]) stop("Number of columns in final constucted x matrix for olddata does not match original")
  
  for(i in 1:dim(x_old)[2]){
    if(!colnames(x_old)[i]==colnames(x)[i]) stop("Column names in final matrix do not match original")
  }
  
  
  x_new2=model.matrix(formula(object$formula),mod_frame_new)
  
  return(list(x_new=x_new2,mod_frame_new=mod_frame_new))
}



#' Generates predictions used by the glm predict function
#'
#' Simulates New Data tied to the model predictions 
#' @param object an object of class \code{glmb}
#' @param type type of response. Either \code{"link"}, the default, or \code{"response"}.
#' @param new_x new design matrix corresponding to newdata passed to the predict function
#' @param new_model_frame new model frame corresponding to newdata passed to the predict function
#' @return A matrix \code{pred} with model predictions where the number of rows equal the number of 
#' rows in \code{object$coefficients} and the number of columns equal the number of rows in  
#' newdata passed to the predict function.
#' @example inst/examples/Ex_confint.glmb.R
#' @noRd 

generate_predictions<-function(object,type="link",new_x,new_model_frame){
  
  ## Set dimensions
  
  nvars=ncol(object$coefficients)
  n=nrow(object$coefficients)
  npreds=nrow(new_x)    
  
  betas=object$coefficients
  mod_formula=formula(object)
  
  new_terms=terms(formula(object))
  
  ## Initialized offset term to 0
  new_alpha=matrix(0,nrow=npreds,ncol=1)
  
  ## If offsets exists, add together into the new offset term
  if(!is.null(attr(new_terms,"offset"))){
    new_offsets=attr(new_terms,"offset")
    n_offsets=length(attr(new_terms,"offset"))
    for(j in 1:n_offsets){new_alpha[1:npreds,1]=new_alpha[1:npreds,1]+new_model_frame[1:npreds,new_offsets[j]]}        
  }
  
  # Calculate predictions and if needed transform to response scale
  
  #print(new_alpha)
  #print(new_x)
  
  
  pred=matrix(0,nrow=n,ncol=npreds)
  
  for(i in 1:n)
  {
    betas_temp=matrix(betas[i,1:nvars],nrow=nvars,ncol=1)
    
    #print(betas_temp)
    
    pred[i,1:npreds]=t(new_alpha+new_x%*%betas_temp)
    
    ## Rescale predictions if type="response"
    ## May need adjustment for offset 
    
    if(type=="response") {
      pred[i,1:npreds]<- family(object)$linkinv(pred[i,1:npreds])
    }
  }
  
  return(pred=pred)  
}

