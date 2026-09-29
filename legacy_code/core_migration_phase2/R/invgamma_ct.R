#' The Central Inverse-Gamma Distribution
#'
#' Distribution function, quantile function, and random generation for the
#' inverse-Gamma distribution on the dispersion scale. These functions provide
#' numerically stable evaluation and sampling when the dispersion parameter is
#' restricted to lie between a lower and an upper bound.
#'
#' @name InvGamma_ct
#' @aliases
#' InvGamma_ct
#' pinvgamma_ct
#' qinvgamma_ct
#' rinvgamma_ct
#'
#' @param dispersion Value(s) at which the inverse-Gamma distribution function
#'   is evaluated.
#' @param p Probability value(s) for the quantile function.
#' @param n Number of random draws to generate. If \code{length(n) > 1}, the
#'   length is taken to be the number required.
#' @param shape Shape parameter of the inverse-Gamma distribution.
#' @param rate Rate parameter of the inverse-Gamma distribution.
#' @param disp_lower Lower bound of the dispersion parameter.
#' @param disp_upper Upper bound of the dispersion parameter.
#'
#' @return
#' For \code{pinvgamma_ct}, a vector of distribution function values evaluated at
#' \code{dispersion}.  
#' For \code{qinvgamma_ct}, a vector of quantiles corresponding to the
#' probabilities \code{p}.  
#' For \code{rinvgamma_ct}, a vector of length \code{n} containing random draws
#' from the inverse-Gamma distribution restricted to the interval
#' \code{[disp_lower, disp_upper]}.
#'
#' @details
#' The inverse-Gamma distribution is defined by the transformation
#' \eqn{D = 1 / X}, where \eqn{X} follows a Gamma distribution with the same
#' shape and rate parameters. The functions \code{pinvgamma_ct} and
#' \code{qinvgamma_ct} therefore compute probabilities and quantiles by mapping
#' the dispersion value \code{D} to the corresponding Gamma scale and applying
#' the Gamma CDF or quantile function.
#'
#' The function \code{rinvgamma_ct} generates random draws from a truncated
#' inverse-Gamma distribution by sampling a uniform probability and inverting
#' the truncated CDF on the Gamma scale. This approach avoids numerical
#' instability when the truncation interval is narrow or when the dispersion
#' parameter is close to zero.
#'
#' These functions are primarily intended for hierarchical Bayesian models in
#' which dispersion parameters are updated under tight truncation constraints.
#' They provide a stable alternative to direct manipulation of the Gamma
#' distribution when working on the dispersion scale is more natural or more
#' numerically robust. They are used in envelope-based dispersion sampling
#' \insertCite{Nygren2006}{glmbayes}.
#'
#' @seealso \code{\link{Gamma_ct}}, \code{\link{Normal_ct}}, \code{\link{EnvelopeDispersionBuild}}
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_InvGamma_ct.R
#' @rdname InvGamma_ct
#' @export
pinvgamma_ct <- function(dispersion, shape, rate) {
  1 - pgamma(1 / dispersion, shape = shape, rate = rate)
}

#' @rdname InvGamma_ct
#' @export
qinvgamma_ct <- function(p, shape, rate, disp_upper, disp_lower) {
  p_upp <- pinvgamma_ct(disp_upper, shape = shape, rate = rate)
  p_low <- pinvgamma_ct(disp_lower, shape = shape, rate = rate)
  p1    <- p_low + p * (p_upp - p_low)
  p2    <- 1 - p1
  1 / qgamma(p2, shape, rate)
}

#' @rdname InvGamma_ct
#' @export
rinvgamma_ct <- function(n, shape, rate, disp_upper, disp_lower) {
  p <- runif(n)
  qinvgamma_ct(p = p, shape = shape, rate = rate,
               disp_upper = disp_upper, disp_lower = disp_lower)
}


