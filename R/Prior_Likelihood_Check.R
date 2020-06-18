#' Checks for Prior-data conflicts
#'
#' Checks if the credible intervals for the prior overlaps with the implied confidence intervals from 
#' the classical model that comes from a call to the glm function 
#' @param object a fitted model object of class \code{"glm"}. Typically the result of a call to \link{glm}.
#' @param mu a vector of length \code{p} giving the prior means of the variables in the design matrix.
#' @param Sigma a positive-definite symmetric matrix of dimension \code{p * p} specifying the prior covariance matrix of the variables.
#' @param level the confidence level at whihc the Prior-data conflict should be checked.
#' @return A vector where each item provided the ratio of the absolue value for the difference between the 
#' prior and maximum likelihood estimate divided by the length of the sum of half of the two intervals 
#' (where normality is assumed)
#' @example inst/examples/Ex_Prior_Check.R


Prior_Check<-function(object,mu=NULL, Sigma=NULL,level=0.95){

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

