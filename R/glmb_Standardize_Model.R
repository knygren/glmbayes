#' Standardize A Non-Gaussian Model
#'
#' Standardizes a Non-Gaussian Model prior to Envelope Creation
#' @param y a vector of observations of length m
#' @param x a design matrix of dimension m*p
#' @param P a positive-definite symmetric matrix specifying the prior precision matrix of the variables.
#' @param bstar a matrix containing the posterior mode from an optimization step
#' @param A1 a matrix containing the posterior precision matrix at the posterior mode
#' @return A list with the following components
#' \item{bstar2}{Standardized Posterior Mode}
#' \item{A}{Standardized Data Precision Matrix}
#' \item{x2}{Standardized Design Matrix}
#' \item{mu2}{Standardized Prior Mean vector}
#' \item{P2}{Standardized Precision Matrix Added to log-likelihood}
#' \item{L2Inv}{A matrix used when undoing the first step in standardization described below}
#' \item{L3Inv}{A matrix used when undoing the second step in standardization described below}
#' @details This functions starts with basic information about the model in the argument list and then
#' uses the following steps to further standardize the model (the model is already assumed to have a 0 prior mean vector
#' when this step is applied).
#' 
#' 1) An eigenvalue composition is applied to the posterior precision matrix, and the model is (as an interim step)
#' standardized to have a posterior precision matrix equal to the identity matrix. Please note that this means
#' that the prior precision matrix after this step is \code{"smaller"} than the identity matrix.
#' 
#' 
#' 2) A diagonal matrix epsilon is pulled out from the standardized prior precision matrix so that the remaining
#' part of the prior precision matrix still is positive definite. That part is then treated as part of the posterior
#' for the rest of the standardization and simulation and only the part connected to epsilon is treated as part of the prior. 
#' Note that the exact epsilon chosen seems not to matter. Hence there are many possible ways of doing this 
#' standardization and future versions of this package may tweak the current approach 
#' if it helps improve numerical accuracy or acceptance rates.
#' 
#' 
#' 3) The model is next standardized (using a second eigenvalue decomposition) so that the prior (i.e., the portion connected to epsilon) is the identity 
#' matrix. The standardized model then simutaneously has the feature that the prior precision matrix is the 
#' identity matrix and that the data precision A (at the posterior mode) is a diagonal matrix. Hence the variables
#' in the standardized model are approximately independent at the posterior mode.
#' @example inst/examples/Ex_glmb_Standardize_Model.R
#' @export


glmb_Standardize_Model<-function(y, x, P, bstar, A1){
  
  return(.glmb_Standardize_Model_cpp(y, x, P, bstar, A1))

  }
  
