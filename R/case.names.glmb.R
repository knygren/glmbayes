#' Case and Variable Names of Fitted Models
#'
#' Simple utilities returning (non-missing) canse names, and (non-eliminated) variable names.
#' @param object a fitted model object, typically of class \code{"glmb"}. 
#' @param full logical; if \code{TRUE}, all names (including zero weights, ...) are returned.
#' @param \ldots further arguments passed to or from other methods.
#' @return A character vector
#' @example inst/examples/Ex_confint.glmb.R
#' @rdname case.names.glmb
#' @export
#' @method case.names glmb


case.names.glmb<-function(object,full=FALSE,...){
  w=weights(object)
  dn=rownames(object$x)
  if(full || is.null(w)) dn else dn[w!=0]
}

#' @rdname case.names.glmb
#' @export
#' @method variable.names glmb

variable.names.glmb <- function(object, full = FALSE, ...)
{
  # Currently just returns results based on glm object
  return(colnames(object$coefficients))
}