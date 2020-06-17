#' Model Formulae
#'
#' The generic function \code{formula} and its specific methods provide a
#' way of extracting formulae which have been included in other objects.
#' @param x \R object.
#' @param \ldots further arguments passed to or from other methods.
#' @return Produce an object of class \code{"formula"} which
#' contains a symbolic model formula.
#' @example inst/examples/Ex_formula.glmb.R

formula.glmb<-function(x,...)
{
  return(formula(x$glm))
  
}
