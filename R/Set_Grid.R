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
#' @return A list with the following components (add).
#' @example inst/examples/Ex_extractAIC.glmb.R
#' @export

Set_Grid<-function(GridIndex,cbars,Lint){
  
  return(.Set_Grid_cpp(GridIndex,cbars,Lint))
  
}