#' Summarizing Bayesian Generalized Linear Model Distribution Functions
#'
#' @inherit glmbayesCore::summary.rglmb return details title description references examples format note
#' @param object an object of class \code{"rglmb"} or \code{"rlmb"} for which a
#'   summary is desired.
#' @param x an object of class \code{"summary.rglmb"} for which a printed output
#'   is desired.
#' @param digits the number of significant digits to use when printing.
#' @param \ldots additional optional arguments.
#' @seealso \code{\link{rglmb}}, \code{\link{rlmb}}, \code{\link{summary.glmb}},
#'   \code{\link{summary.rGamma_reg}}, \code{\link[stats]{summary.glm}}.
#' @aliases summary.rglmb summary.rlmb print.summary.rglmb
#' @name summary.rglmb
NULL

## \code{summary.rglmb} and \code{summary.rlmb} are implemented and registered
## in glmbayesCore; loaded via \code{Imports} when \code{glmbayes} is attached.
## \code{\link{summary.glmb}} remains glmbayes-specific for \code{glmb} fits.
