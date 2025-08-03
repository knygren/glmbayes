#' Derives the inverse of the gradient vector.
#'
#' Computes the inverse of the gradient vector for the gaussian model. Typically done to find the set 
#' of tangency points that yield the same gradient as an initial set of gradients used for an envelope 
#' positioned at the posterior mode for a specific dispersion.
#' @param cbars Gradient vectors desired for new thetabars. Typically the output of an initial call to 
#' Envelopebuild.
#' @param mu Prior mean 
#' @param P Prior Precision matrix
#' @param alpha offset vector
#' @param wt weighting vector
#' @return Refer to Nygren and Nygren (2006) for details. The first set of items refers to Example 2 in section 
#' 3.1. All except the last item in this list of returned items has a number of rows equaling the number of 
#' components of the grid and a number of columns equaling the number of coefficients in the model. All quantities 
#' refer to the respective coefficient for each of the components of the grid.
#' \item{Down}{The lower bounds for the interval to be evaluated. Either negative infinity or a real number.}
#' \item{Up}{The upper bounds for the interval to be evaluated. Either positive infinity or a real number.}
#' \item{lglt}{The log of the density between negative infinity and the upper bound}
#' \item{lgrt}{The log of the density between the lower bound and infinity}
#' \item{lgct}{The log of the density between the lower and upper bounds}
#' \item{logU}{The one of the 3 above that is relevant for the component of the grid}
#' \item{logP}{A two column matrix, the first of which holds sum of logU across the components. The second column is 
#' 0 and is later populated by the Set_logP function}
#' @inheritParams  glmb
#' @example inst/examples/Ex_extractAIC.glmb.R
#' @export

Inv_f3_gaussian<-function(cbars, y, x, mu, P, alpha, wt){
  
  return(.Inv_f3_gaussian(cbars, y, x, mu, P, alpha, wt))
  
}

