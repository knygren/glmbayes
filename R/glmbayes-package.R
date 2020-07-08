#' Bayesian Generalized Linear Models (iid Samples)
#'
#' Generates iid samples for Bayesian Generalized Linear Models.
#' @aliases
#' glmbayes-package
#' glmbayes
#' @details 
#' The \code{glmbayes} package produces iid samples for Bayesian Genereralized Linear Models and is
#' intended as a Bayesian version of the \code{\link{glm}} function for classical models. Estimation can be performed
#' using three main functions. For models with fixed dispersion parameters, the \code{\link{rglmb}} 
#' function is the workhorse function and comes with a minimialistic interface for the input and output.
#' It is also suitable for use as part of block Gibbs sampling procedures. The \code{\link{glmb}} function is 
#' essentially a wrapper function for the rglmb function that provides an interface closer to that of the \code{\link{glm}} 
#' function. The \code{\link{rGamma_reg}} function can be leveraged in order to produce samples for the 
#' dispersion parameters associated with the gaussian and Gamma link functions. Most methods 
#' defined for the output of the \code{\link{glm}} function are also defined for the \code{\link{glmb}}, 
#' \code{\link{rglmb}}, and \code{\link{rGamma_reg}} functions (see their respective documentation 
#' for details).
#' 
#' For the regression parameters, multivariate normal priors are assumed. Simulation for the 
#' gaussian family with the identify link function is performed using standard procedures for 
#' multivariate normal densities.  For all other families and link functions, simulation is performed using 
#' the likelihood subgradient approach of Nygren and Nygren (2006). This approach involves the 
#' construction of an enveloping function for the full posterior density followed by accept-reject 
#' sampling. For models that are approximately multivariate normal, the expected number of draws 
#' required per acceptance are bounded from above as noted in Nygren and Nygren (2006). 
#' 
#' Currently implemented models include the gaussian (identity link), poisson/quasipoisson (log link), 
#' binomial/quasibinomial (logit, probit, and cloglog links), and Gamma (log link) families. These 
#' models all have log-concave likelihood functions that allow us to leverage the likelihood-subgradient 
#' approach for the iid sampling. Models that fail to have log-concave likelihood functions are not 
#' implemented.  Our demos (viewable by entering the \code{demo()} command) provides examples of each 
#' of these families and links. 
#' 
#' The current implementation requires separate use of the \code{\link{rGamma_reg}} function in order 
#' to generate samples for dispersion parameters (gaussian, Gammma, quasipoisson, quasi-binomial 
#' families). Our demos include examples of the joint use of the \code{\link{rglmb}} and \code{\link{rGamma_reg}} to 
#' produce samples for both regression and dispersion parameters using two-block Gibbs samplers. 
#' As these two-block Gibbs samplers likely are geometrically ergodic, future implementations may 
#' incorporate these two-block Gibbs samplers into the \code{\link{rglmb}} and \code{\link{glmb}} functions by leveraging 
#' theoretical bounds om convergence rates derived using Rosenthal (1996) type drift and 
#' minorization conditions.
#' 
#' The \code{\link{rglmb}} function can also be used in Block-Gibbs sampling implementations for Hierarchical 
#' Bayesian models. The demos associated with this package contains examples of such models.
#'  
#'
#' @references 
#' 
#' Dobson, A. J. (1990)
#' \emph{An Introduction to Generalized Linear Models.}
#' London: Chapman and Hall.
#' 
#' Hastie, T. J. and Pregibon, D. (1992)
#' \emph{Generalized linear models.}
#' Chapter 6 of \emph{Statistical Models in S}
#' eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#' McCullagh P. and Nelder, J. A. (1989)
#' \emph{Generalized Linear Models.}
#' London: Chapman and Hall.
#' 
#' Nygren, K.N. and Nygren, L.M (2006)
#' Likelihood Subgradient Densities. \emph{Journal of the American Statistical Association}.
#' vol.101, no.475, pp 1144-1156.
#' 
#' Raiffa, Howard and Schlaifer, R (1961)
#' \emph{Applied Statistical Decision Theory.}
#' Boston: Clinton Press, Inc.
#' 
#' Venables, W. N. and Ripley, B. D. (2002)
#' \emph{Modern Applied Statistics with S.}
#' New York: Springer.
#' 
#' 
#' @example inst/examples/Ex_glmbayes-package.R
#' @docType package
#' @author Kjell Nygren
#' @import stats Rcpp RcppArmadillo
#' @importFrom Rcpp evalCpp
#' @importFrom MASS mvrnorm
#' @useDynLib glmbayes
#' @name glmbayes-package
#' 
NULL