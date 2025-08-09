#' Builds Envelope function for simulation
#'
#' Builds an Enveloping function for simulation using a grid and tangencies for 
#' the posterior density
#' @param bStar Point at which envelope should be centered (typically posterior mode
#' @param A Diagonal precision matrix for log-likelihood function associated with model in standard form
#' @param y a vector of observations of length \code{m}
#' @param x a design matrix of dimension \code{m * p}
#' @param mu a vector giving the prior means of the variables
#' @param P a positive-definite symmetric matrix specifying the prior precision matrix of the variables
#' @param alpha offset vector 
#' @param wt a vector of weights
#' @param family Family for the envelope. Can take on values binomial, quasibinomial, poisson, quasibinomial, and Gamma
#' @param link Link function for the envelope. Valid links for the binomial and quasibinomial families are logit, probit, and cloglog. The valid value for the poisson, quasipoisson, and the Gamma families is the log
#' @param Gridtype an optional argument specifying the method used to determine number of likelihood subgradient 
#' densities used to construct the enveloping function
#' @param n number of draws to generate from posterior density. Used here to help determine the size of the grid (Gridtype 1 and 2 only)
#' @param sortgrid an optional Logical argument determining whether the final envelope should be sorted in descending order based on probability for each part of the envelope
#' @return The function returns a list consisting of the following components (
#' the first six of which are matrices with the number of rows equal to the 
#' number of components in the Grid and columns equal to the number of parameters):
#' \item{GridIndex}{A matrix containing information on how each dimension should be 
#' sampled (1 means left tail of a restricted normal, 2 center, 3 right tail, and 4 the 
#' entire line)}
#' \item{thetabars}{A matrix containing the points of tangencies associated with each component of the grid}
#' \item{cbars}{A matrix containing the gradients for the negative log-likelihood at each tangency}
#' \item{logU}{A matrix containing the log of the cummulative probability associated with each dimension}
#' \item{logrt}{A matrix containing the log of the probability associated with the right tail (i.e. 
#' that to the right of the lower bound)}
#' \item{loglt}{A matrix containing the log of the probability associated with the left tail (i.e.,
#' that to the left of the upper bound)}
#' \item{LLconst}{A vector containing constant for each component of the grid used during the accept-reject procedure}
#' \item{logP}{A matrix containing log-probabilities related to the components of the grid}
#' \item{PLSD}{A vector containing the probability of each component in the Grid}
#' @details To construct an enveloping function, we follow the approach in Nygren and Nygren (2006)
#' which involves the following steps when a maximally sized grid is constructed 
#' (if the prior for some dimensions is relatively strong, this may not be needed)
#' 
#'
#' 1) For each dimension, a constant \code{omega_i} is found that depends on the 
#' corresponding diagonal element in the precision matrix.
#' 
#'  
#' 2) Corresponding intervals \code{(thetastar_i-0.5 *omega,thetastar_i-0.5 *omega)} 
#'  are constructed around the posterior mode thetastar for each dimension
#'  
#'  
#' 3) The mode as well as the points \code{thetastar_i-omega_i} and \code{thetastar_i+omega_i}
#' are selected as the components of the points at which tangencies will be found for each of 
#' the dimensions.  
#' 
#' 
#' 4) A Grid is constructed with all possible combinations of points and 
#' negative log-likelihood and gradient for the negative log-likelihood are evaluated (see the
#' EnvelopeBuild_c.cpp function source code for details)
#' 
#' 
#' 5) The \code{\link{Set_Grid}} function is called in order to evaluate the log of the density associated
#' with each of the resulting restricted multivariate normals by evaluating the differences between the cummulative 
#' density for each dimension between its lower and upper bound.
#' 
#' 
#' 6)  The Set_LogP is called in order to help set the probabilities with which each of the components
#' of the grid should be sampled (see remark 6 in Nygren and Nygren (2006)).
#' 
#' 
#' Any constants needed by the sampling are added to a list and returned by the function. 
#'     
#' @example inst/examples/Ex_EnvelopeBuild.R
#' @export
#' @keywords internal



EnvelopeBuild<-function(bStar, A, y, x, mu, P, alpha, wt, family = "binomial",
                         link = "logit", Gridtype = 2L, n = 1L, sortgrid = FALSE,
                        use_opencl = FALSE,        
                        verbose = FALSE 
                        ){
  
  if(family=="gaussian"){
    return(.EnvelopeBuild_Ind_Normal_Gamma(bStar, A, y, x, mu, P, alpha, wt, family = family,
                         link = link, Gridtype = Gridtype, n = n, sortgrid))}
  
  else{
  return(.EnvelopeBuild_cpp(bStar, A, y, x, mu, P, alpha, wt, family = family,
                            link = link, Gridtype = Gridtype, n = n, sortgrid,use_opencl=use_opencl,verbose=verbose ))}
}
  
