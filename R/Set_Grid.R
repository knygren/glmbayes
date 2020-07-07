#' Calculate Log-densities for Grid Components
#'
#' Computes un-normalized log of the density associated with each component of a Grid
#' used to sample using the likelihood-subgradient density approach
#' @param GridIndex A matrix containing information for each component in the grid related
#' to whether the components is to the left, in the center, or to the right of the density. 
#' Each row corresponds to a component in the grid, while the columns correspond to the 
#' transformed (standardized) variables.
#' @param cbars A matrix containing the subgradient for the (adjusted) negative log-likelihood function
#' for each the component in the grid. 
#' @param Lint A matrix storing information on where the upper and lower bounds are, depending
#' on whether the sampling is from the left, the center, or the right.
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
#' @example inst/examples/Ex_extractAIC.glmb.R
#' @export
#' @keywords internal

Set_Grid<-function(GridIndex,cbars,Lint){
  
  return(.Set_Grid_cpp(GridIndex,cbars,Lint))
  
}

