#' Setup Prior Objects
#'
#' Sets up the structure for the Prior mean and Variance Matrices using information from a classical model.
#' @param object a fitted model object of class \code{"glm"}. Typically the result of a call to \link{glm}.
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


Prior_Setup<-function(object){

    nvar=length(object$coefficients)
    var_names=names(object$coefficients)
    mu=matrix(0,nrow=nvar,ncol=1)
    Sigma=as.matrix(diag(nvar))
    
    rownames(mu)=var_names
    colnames(mu)=c("mu")
    rownames(Sigma)=var_names
    colnames(Sigma)=var_names
    
    return(list(mu=mu,Sigma=Sigma,model=as.data.frame(object$model),x=as.matrix(object$x)))    
    
  }

#' @export
#' @rdname Prior_Setup
#' @order 2

Prior_Setup_v2<-function(formula,data, subset = NULL, na.action = na.fail, 
                         drop.unused.levels = FALSE, xlev = NULL, ...){
  
  
  mf<-model.frame(formula)
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
  
  return(list(mu=mu,Sigma=Sigma,model=as.data.frame(object$model),x=as.matrix(object$x)))    
  
}

