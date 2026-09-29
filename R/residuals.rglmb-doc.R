#' Deviance Residuals for \code{rglmb} and \code{summary.rglmb} Objects
#'
#' @inherit glmbayesCore::residuals.rglmb return title description
#' @param object an object of class \code{rglmb}, \code{rlmb}, or
#'   \code{summary.rglmb}.
#' @param ysim optional matrix of simulated responses (one row per draw), as
#'   produced by a posterior-predictive simulation.
#' @param ... further arguments (currently unused).
#' @seealso \code{\link{residuals.glmb}}, \code{\link{rglmb}}, \code{\link{rlmb}},
#'   \code{\link{summary.rglmb}}, \code{\link[stats]{residuals.glm}}.
#' @name residuals.rglmb
NULL

## \code{residuals} methods for \code{rglmb}, \code{rlmb}, and
## \code{summary.rglmb} are registered from glmbayesCore when \code{glmbayes}
## is loaded. \code{\link{residuals.glmb}} remains in this package.
