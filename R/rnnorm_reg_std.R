#' The Bayesian Generalized Linear Model Distribution in Standard Form
#'
#' \code{rnnorm_reg_std} is used to generate iid samplers from Non-Gaussian Generalized
#' Linear Models in standard form. The function should onlybe called after standardization 
#' of a Generalized Linear Model.
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x a design matrix of dimension \code{m * p}.
#' @param mu a vector of length \code{p} giving the prior means of the variables in the design matrix.
#' @param P a positive-definite symmetric matrix of dimension \code{p * p} specifying the prior precision
#'  matrix of the variable.
#' @param alpha this can be used to specify an \emph{a priori} known component 
#' to be included in the linear predictor during fitting. This should be 
#' \code{NULL} or a numeric vector of length equal to the number of cases. 
#' One or more offset terms can be included in the formula instead or as well,
#' and if more than one is specified their sum is used. See 
#' \code{\link{model.offset}}.
#' @param wt an optional vector of \sQuote{prior weights} to be used in the fitting process. 
#' Should be NULL or a numeric vector.
#' @param f2 function used to calculate the negative of the log-posterior function
#' @param Envelope an object of type \code{glmbenvelope}.
#' @param family family used for simulation. Used that this is different from the 
#' family used in other functions.
#' @param link link function used for simulation.
#' @param progbar dummy for flagging if a progressbar should be produced during the call
#' @return A list consisting of the following:
#' \item{out}{A matrix with simulated draws from a model in standard form. Each row represents 
#' one draw from the density}
#' \item{draws}{A vector with the number of candidates required before 
#' acceptance for each draw}
#' @details Simulates from the posterior density of a model in standard form
#' @example inst/examples/Ex_rnnorm_reg_std.R
#' @export


rnnorm_reg_std<-function(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link, progbar = 1L){
  
  return(.rnnorm_reg_std_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link, progbar = progbar))

  }
  
