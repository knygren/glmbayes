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

  