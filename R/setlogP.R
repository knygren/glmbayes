#' Calculate constants used during sampling from Likelihood subgradient densities
#'
#' This function computes two vectors of constants needed during sampling. The first is used to determine
#' the probability with which each component of the grid should be visited, while the second 
#' is used as a constants when calculating acceptance rates.
#' @param logP A matrix that typically contains two columns and information for each component of the grid. The 
#' first column will typically hold the final output from the Set_Grid function, which is the density
#' associated with the related restricted normal.
#' @param NegLL A vector with evaluations of the Negative Log-likelihood for each of the components of the grid. 
#' @param cbars A matix holding the gradients for the Negative of the Log-likelihood for each of the 
#' componentsof the grid.
#' @param G3 A matrix containing the set of tangent points used in the grid.
#' @return Refer to Nygren and Nygren (2006) for details. The first .
#' \item{logP}{The first column holds the value passed into the function while the second contains the log of the 
#' (un-normalized) probabilities with which each of the components of the grid should be visited. This corresponds 
#' to the log of the denominator components used to compute p_i in remark 6 in Nygren and Nygren (2006)}
#' \item{LLconst}{This holds a vector of constants used as upper bounds when deriving acceptance rates during the 
#' accept-reject sampling process. This constant corresponds to the denominator for the function h() in Theorem 1 
#' of Nygren and Nygren (2006). During the sampling, the log of the numerator of the same function if evaluated for
#' each candidate and the difference between the candidate value and this constant is used to determine the acceptance
#' rate to use when evaluating acceptance of the candidate. If the evaluated value for the candidate is close to 
#' this constant, then the chance of acceptance rate is high. If it is much smaller, then the chance of acceptance 
#' is low.}
#' @example inst/examples/Ex_extractAIC.glmb.R
#' @export


setlogP<-function(logP, NegLL, cbars, G3){
  
return(.setlogP_cpp(logP, NegLL, cbars, G3))
  
}