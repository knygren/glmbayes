#' Accessing Bayesian Generalized Linear Model Fits
#'
#' Extract deviance residuals from fitted Bayesian GLM objects. The residuals
#' use the family's deviance residuals function as in \code{\link[stats]{residuals.glm}}
#' \insertCite{McCullagh1989}{glmbayes}.
#'
#' These functions are \link{methods} for class \code{glmb} and \code{lmb}
#' objects. Methods for \code{rglmb}, \code{rlmb}, and \code{summary.rglmb}
#' are registered from \code{glmbayesCore} when this package is loaded; see
#' \code{\link{residuals.rglmb}}.
#' @param object an object of class \code{glmb} or \code{lmb}, typically the
#'   result of a call to \link{glmb} or \link{lmb}.
#' @param ysim Optional matrix of simulated responses (one row per draw), as
#'   produced by a posterior-predictive simulation (e.g. \code{\link{simulate.glmb}}).
#'   When supplied, \code{ysim} substitutes for the observed response \code{y}
#'   while the fitted value for each draw is held fixed at that draw's own
#'   fit; this builds a reference ("what would a typical residual look like
#'   under the model") distribution for posterior-predictive residual checks,
#'   comparable to the residuals obtained from the actual data
#'   (\code{ysim = NULL}).
#' @param \ldots further arguments to or from other methods
#' @return A matrix \code{DevRes} of dimension \code{n} times \code{p} containing
#' the Deviance residuals for each draw. If ysim is provided, the residuals are based
#' on a comparison to the simulated data instead. The credible intervals
#' for residuals based on simulated data should be a more appropriate measure of
#' whether individual residuals represent outliers or not.
#' @seealso \code{\link{predict.glmb}}, \code{\link{summary.glmb}}, \code{\link{glmb}},
#'   \code{\link{glmbayes-package}}; \code{\link{residuals.rglmb}};
#'   \code{\link{rglmb}}, \code{\link{rlmb}}, \code{\link{lmb}};
#'   \code{\link[stats]{residuals.glm}}
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_residuals.glmb.R
#' @export 
#' @method residuals glmb

## This method follows stats::residuals.glm() conventions while computing
## residuals across posterior draws. See inst/COPYRIGHTS.
residuals.glmb<-function(object,ysim=NULL,...)
{
  .residuals_rglmb_draws(object, ysim = ysim)
}


#' @rdname residuals.glmb
#' @export 
#' @method residuals lmb

residuals.lmb<-function(object,ysim=NULL,...)
{
  return(residuals.lm(object,ysim,...))
  }


## Shared draw-wise deviance-residual computation for glmb/rglmb/rlmb/
## summary.rglmb objects. `glmb` and `summary.rglmb` objects carry a
## precomputed `fitted.values` matrix; plain `rglmb`/`rlmb` objects do not, so
## the linear predictor and fitted values are recomputed from `x` and
## `coefficients` in that case.
##
## `ysim` (when supplied) substitutes for the observed response y, with the
## fitted value mu held fixed at each draw's own fit -- not the reverse. This
## is the posterior-predictive-check convention demonstrated in
## vignette("Chapter-05", package = "glmbayes") (simulate new data via
## `simulate.glmb()`, recompute residuals against the same fitted values, and
## compare to the actual residuals to assess whether they look unusual under
## the model). Matches glmbayesCore's `.residuals_rglmb_draws()`.
.residuals_rglmb_draws <- function(object, ysim = NULL) {
  y <- object$y
  n <- nrow(object$coefficients)
  wts <- object$prior.weights

  if (!is.null(object$fitted.values)) {
    fv_mat <- object$fitted.values
  } else {
    lp_mat <- t(object$x %*% t(object$coefficients))
    fv_mat <- object$family$linkinv(lp_mat)
  }

  devfun <- object$family$dev.resids
  DevRes <- matrix(0, nrow = n, ncol = length(y))

  for (i in seq_len(n)) {
    mu_vec <- fv_mat[i, ]
    y_vec <- if (is.null(ysim)) y else ysim[i, ]
    DevRes[i, ] <- sign(y_vec - mu_vec) * sqrt(devfun(y_vec, mu_vec, wts))
  }

  colnames(DevRes) <- names(y)
  DevRes
}
