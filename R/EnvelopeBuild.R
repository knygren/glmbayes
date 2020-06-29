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
#' @example inst/examples/Ex_EnvelopeBuild.R
#' @export



EnvelopeBuild<-function(bStar, A, y, x, mu, P, alpha, wt, family = "binomial",
                         link = "logit", Gridtype = 2L, n = 1L, sortgrid = FALSE){
  
  if(family=="gaussian"){
    return(.EnvelopeBuild_Ind_Normal_Gamma(bStar, A, y, x, mu, P, alpha, wt, family = family,
                         link = link, Gridtype = Gridtype, n = n, sortgrid))}
  
  else{
  return(.EnvelopeBuild_cpp(bStar, A, y, x, mu, P, alpha, wt, family = family,
                            link = link, Gridtype = Gridtype, n = n, sortgrid))}
}
  
