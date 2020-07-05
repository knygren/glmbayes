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
#' @example inst/examples/Ex_Prior_Check.R
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
  Sigma=as.matrix(diag(nvar))
  
  rownames(mu)=var_names
  colnames(mu)=c("mu")
  rownames(Sigma)=var_names
  colnames(Sigma)=var_names
  
  print("Variable names are:")
  print(var_names)
  return(list(mu=mu,Sigma=Sigma,model=mf,x=x))    
  
}
