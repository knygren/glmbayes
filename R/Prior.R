#' Setup Prior Objects
#'
#' Sets up the structure for the Prior mean and Variance Matrices using information from a classical model.
#' @param na.action how \code{NAs} are treated. The default is first, any \code{\link{na.action}} attribute of 
#' data, second a \code{na.action} setting of \link{options}, and third \code{na.fail} if that is unset. 
#' The \code{factory-fresh} default is \code{na.omit}. Another possible value is \code{NULL}.
#' @inheritParams stats::model.frame
#' @return A list with items related to the prior.
#' \item{mu}{An initial version of the prior mean vector, populated with all zeros}
#' \item{Sigma}{An initial version of the prior Variance-Covariance vector, populated as the diagonal identity matrix}
#' \item{model}{The model frame from \code{object} if it exists}
#' \item{x}{The design matrix from \code{object} if it exists}
#' @family prior
#' @seealso 

#' @example inst/examples/Ex_Prior_Setup.R
#' @export
#' @rdname Prior_Setup
#' @order 1

## Note arguments outside of first two are currently not used

Prior_Setup<-function(formula,data=NULL, subset = NULL, na.action = na.fail, 
                         drop.unused.levels = FALSE, xlev = NULL, ...){
  
  #mf<-model.frame(formula,data,subset=subset,na.action=na.action,
  #                drop.unused.levels=drop.unused.levels,xlev=xlev)
  
  mf<-model.frame(formula,data)
  x<-model.matrix(formula,mf)
  
  nvar=ncol(x)
  var_names=colnames(x)
  #nvar=length(object$coefficients)
  mu=matrix(0,nrow=nvar,ncol=1)
  
  if(var_names[1]=='(Intercept)'){
    lm_out=lm(formula,data=mf,y=TRUE)
    y=lm_out$y
    mu[1,1]=mean(y)
    
  } 
  
  
  Sigma=as.matrix(diag(nvar))
  
  rownames(mu)=var_names
  colnames(mu)=c("mu")
  rownames(Sigma)=var_names
  colnames(Sigma)=var_names
  
  print("Variable names are:")
  print(var_names)
  return(list(mu=mu,Sigma=Sigma,model=mf,x=x))    
  
}



#' Checks for Prior-data conflicts
#'
#' Checks if the credible intervals for the prior overlaps with the implied confidence intervals from 
#' the classical model that comes from a call to the glm function 
#' @param level the confidence level at which the Prior-data conflict should be checked.
#' @inheritParams glmb
#' @return A vector where each item provided the ratio of the absolue value for the difference between the 
#' prior and maximum likelihood estimate divided by the length of the sum of half of the two intervals 
#' (where normality is assumed)
#' @family prior
#' @seealso 
#' @example inst/examples/Ex_Prior_Check.R
#' @export
#' @rdname Prior_Check
#' @order 1

Prior_Check<-function(formula,family,pfamily,level=0.95,data=NULL, weights, subset,na.action, 
                      start = NULL, etastart, mustart, offset ,control = list(...) , model = TRUE, 
                      method = "glm.fit",x = FALSE, y = TRUE, contrasts = NULL, ...){
  
  pf=pfamily
  prior_list=pfamily$prior_list
  
  ## For now, the below is really only correct for the dNormal pfamily
  
  mu=prior_list$mu
  Sigma=prior_list$Sigma
  
  
  object=glm(formula=formula,family=family,data=data)
  
  Like_est=object$coefficients
  Like_std=summary(object)$coefficients[,2]
  
  if(is.null(mu)){
    print("No Prior mean vector provided. Variables with needed Priors are:")
    print(names(Like_est))
    names(Like_est)    
    
  }
  
  
  if(level<=0.5) stop("level must be greater than 0.5")
  
  Sigma=as.matrix(Sigma)
  Prior_std=sqrt(diag(Sigma))
  
  print("Variables in the Model Are:")
  print(names(Like_est))
  std_dev_sum=qnorm(level)*(Prior_std+Like_std)
  
  abs_ratio=matrix(rep(0,length(Like_est),nrow=length(Like_est),ncol=1))
  abs_diff=abs(mu-Like_est)
  abs_ratio[1:length(Like_est),]=abs_diff/std_dev_sum
  
  rownames(abs_ratio)=names(Like_est)
  colnames(abs_ratio)=c("abs_ratio")
  max_abs_ratio=max(abs_ratio)
  
  if(max_abs_ratio>1) {
    print("At least one of the maximum likelihood estimates appears to be inconsistent with the prior")
  }
  
  else{
    print("The maximum likelihood estimates for all coefficients appear to be roughly consistent with the prior.")
  }
  return(abs_ratio)
  
}

