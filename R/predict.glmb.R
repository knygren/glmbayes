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

predict.glmb<-function(object,newdata=NULL,type="link",
                       se.fit = FALSE, dispersion = NULL, terms = NULL, 
                       na.action = na.pass,olddata=NULL,...)
{


  ## Note: may need adjustment for case when dispersion is random 
  ## residual.scale be may not be 1 or link-inverse may not be accurate
  ## Possible that this only matters if we want to simulate
  ## from likelihood function, so may be ok. Need to verify.
  
  ## If newdata is missing, simply return the model predictions 
  ## on either the scale of the linear predictors or on the scale
  ## of the response variable (the former is the default)
   
  if (missing(newdata)) {
      if(type=="link")   pred <-  object$linear.predictors
      if(type=="response")  pred <-   object$fitted.values
      }
  
  else{

    ## Function assumes columns the same so data frames can be combined
    ## also assumes no critical variables are missing
    
    if (missing(olddata)) stop("olddata must be provided whenever newdata is provided")
    
    # Short function used to separate model 
    # variables into those used in derivation of explanatory and dependent variables
    
    
    get_vars=function(all_vars,xstring){
      
      ## Initialize model t
      xvars=all_vars
      xflags=list(rep("TRUE",length(xvars)))
      
      check_data=data.frame(xvars=xvars,xflags=xflags)
      colnames(check_data)=c("xvars","xflags")
      
      for(i in 1:length(xvars)){
        if(isFALSE(grepl(xvars[i], xstring))) check_data$xflags[i]=FALSE
      }
      
      xdata=subset(check_data, xflags==TRUE)
      return(unique(xdata$xvars))
      
    }
    
    mod_vars=all.vars(formula(object))
    mod_formula=formula(object)
    
    ystring=as.character(mod_formula[2]) ## Seems to be formula for dependent variables
    xstring=as.character(mod_formula[3]) ## Seems to be formula for explanatory variables

    
    if(length(xstring)>0 & length(mod_vars)>0){  x_vars=get_vars(all_vars=mod_vars,xstring=xstring)}
    else{x_vars=NULL    }
    
    
    y_vars1=setdiff(mod_vars,x_vars)
    
    if(length(ystring)>0 & length(y_vars1)>0){
    y_vars2=get_vars(all_vars=y_vars1,xstring=ystring)
    }
    else{y_vars2=y_vars1}
    
    if(length(y_vars1)==0) {warning("No data variables appear to be used in formula 
for dependent variable. The number of observation set by formula are likely fixed and 
predictions may currently only work if the number of observations in newdata match those in olddata")}

    miss_vars=setdiff(y_vars1,y_vars2)
    
    if(length(miss_vars)>0){  
    warning("Some Model Variables Not Found - May cause unpredictable behavior")
      x_vars=c(x_vars,miss_vars)
            }

    
    
    y_vars=y_vars2
    x_vars=x_vars
    
    ## Subset olddata and newdata to model variables
    
    ## Try with error handling

    
    tryCatch(olddata[mod_vars],error=function(e) {stop("olddata does not contain all model variables")})
    tryCatch(newdata[x_vars],error=function(e) {stop("newdata does not contain all explanatory variables")})

    newvars=names(newdata)
    miss_newvars=setdiff(mod_vars, newvars)  # These should now be any missing dependent variables 

    
    if(length(miss_newvars)>2) stop("newdata has more than two missing dependent variables")

  
    if(length(miss_newvars)>0){
    old_miss_newvars=olddata[miss_newvars]
    m_miss_newvars=round(colMeans(old_miss_newvars))

    for(i in 1:length(names(m_miss_newvars)))
      
    {
      Temp1=matrix(round(colMeans(olddata[miss_newvars])[i]),nrow=nrow(newdata),ncol=1)  
      colnames(Temp1)=names(m_miss_newvars)[i]
      if(i==1) Temp2=Temp1
      else(Temp2=cbind(Temp2,Temp1))
      
    }
    
    newdata[miss_newvars] <- Temp2

    }
    
    ## Look to see if this is producing the desired result
    

    tryCatch(newdata[mod_vars],error=function(e) {stop("Modified newdata does not contain all model variables")})

    olddata=olddata[mod_vars]
    newdata=newdata[mod_vars]
    
    #return(list(olddata=olddata,newdata=newdata))
    
    ## Run check to see if olddata returns same x_data a stored in object
    temp_glm1=glm(object$glm$formula, family = object$glm$family,x=TRUE,data=olddata)

    x_old=object$glm$x
    x_new=temp_glm1$x
    
    if(isFALSE(all.equal(x_new,x_old))) stop("olddata does not yield an x matrix consistent with that 
                                             stored in the original model object")

    get_x_matrix<-function(object,olddata,newdata){
      nrow1=nrow(olddata)
      nrow2=nrow(newdata)
      combodata=rbind(olddata,newdata)
    
      temp_glm=glm(object$formula, family = object$family,x=TRUE,data=combodata)
      
      return(temp_glm$x[(nrow1+1):(nrow1+nrow2),])
    }

    
    x_matrix=get_x_matrix(object$glm,olddata,newdata)
    nvars=ncol(object$coefficients)
    n=nrow(object$coefficients)
    betas=object$coefficients
    npreds=nrow(newdata)
    pred=matrix(0,nrow=n,ncol=npreds)
    
    for(i in 1:n)
    {## Making sure names match
      betas_temp=betas[i,1:nvars]
      pred[i,1:npreds]=t(x_matrix%*%as.matrix(betas_temp,ncol=1))
      
      ## Rescale predictions if type="response"
      if(type=="response") pred[i,1:npreds]<- family(object$glm)$linkinv(pred[i,1:npreds])
      
    }
    

  }
  return(pred)
  
}
