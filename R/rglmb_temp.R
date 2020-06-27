#' The Bayesian Generalized Linear Model Distribution
#'
#' \code{rglmb} is used to generate iid samples for Bayesian Generalized Linear Models.
#' The model is specified by providing a data vector, a design matrix, 
#' the family (determining the likelihood function) and the pfamily (determining the 
#' prior distribution).
#' @param pfamily the prior family to use for the model (including the constants passed to prior). 
#' @param offset an offset parameter
#' @param weights a weighting variable
#' @param control2 a list of parameters for controlling the Bayesian fitting process.
#' @inheritParams glmb
#' @return Currently mainly the draws for the dispersion and the regression coefficients
#' will be updated to return outputs consistent with other function
#' @example inst/examples/Ex_rindep_norm_gamma_reg.R
#' @export 



rglmb_temp<-function(n=1,y,x,family=gaussian(),pfamily,offset=NULL,weights=1,control2=list(...)){
  
  ## Pull in information from the pfamily  
  pf=pfamily$pfamily
  okfamilies=pfamily$okfamilies  
  plinks=pfamily$plinks
  prior_list=pfamily$prior_list 
  simfun=pfamily$simfun
  
  ## Pull in information on families  
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  ## Check that the family is implemented for the pfamily
  
  if(family$family %in% okfamilies){
    oklinks=plinks(family)
    if(!family$link %in% oklinks){      
      stop(gettextf("link \"%s\" not available for selected pfamily/family combination; available links are %s", 
                    family$link , paste(sQuote(oklinks), collapse = ", ")), domain = NA)
    }
  }
  else{
    stop(gettextf("family \"%s\" not available for current pfamily; available families are %s", 
                  family$family , paste(sQuote(okfamilies), collapse = ", ")), 
         domain = NA)
    
  }
  
  
  ## Call relevant simulation function (for now without control2 list)
  
  outlist=simfun(n=n,y=y,x=x,prior_list=prior_list,offset=offset,weights=weights,family=family)
  
  return(outlist)
  
}
