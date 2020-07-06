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



