#' bayestestR Prior-Checking Methods for glmb/lmb Model Objects
#'
#' @description
#' Methods implementing three of \pkg{bayestestR}'s prior-related generics for
#' objects of class \code{"glmb"} (which includes \code{\link{lmb}},
#' \code{\link{rglmb}}, \code{\link{rlmb}}, and \code{\link{rGamma_reg}} fits,
#' since these all inherit \code{"glmb"} in their class vector).
#'
#' \code{simulate_prior.glmb()} delegates to the \code{pfun} stored in
#' \code{x$pfamily} (see \code{\link{pfamily}} and \code{\link{prior_simfuncs}}),
#' drawing \code{n} independent samples directly from the \emph{actual} prior
#' specification in \code{prior_list}, rather than reconstructing an approximation
#' \code{\link{glmbayes_insight_methods}}'s \code{get_priors.glmb()} flattened
#' \code{Parameter}/\code{Distribution}/\code{Location}/\code{Scale} table.
#' This matters when the coefficient prior is a genuine \emph{multivariate}
#' normal (\code{\link{dNormal}}, \code{\link{dNormal_Gamma}},
#' \code{\link{dIndependent_Normal_Gamma}}) with a non-diagonal covariance
#' matrix \code{Sigma}: the flattened table can only report each parameter's
#' marginal mean/SD, but \code{simulate_prior.glmb()} draws from the true
#' joint \eqn{N(\mu, \Sigma)} (via a Cholesky factor of the complete
#' \code{Sigma}), preserving any correlation between coefficients.
#' \code{\link{dBeta}}/\code{\link{dGamma}} priors are drawn from their exact
#' Beta/Gamma form rather than their Normal-moment surrogate. \code{dGamma}
#' with \code{Inv_Dispersion = TRUE} is a special case: that prior is on the
#' inverse dispersion (precision) \emph{only} -- \code{beta} is a fixed,
#' known input, not a modeled quantity -- so \code{simulate_prior.glmb()}
#' returns \emph{only} simulated dispersion draws for it, with no
#' coefficient columns at all (mirroring
#' \code{\link{glmbayes_insight_methods}}'s \code{get_parameters.glmb()}
#' treatment of the posterior side).
#'
#' \code{\link{dNormal_Gamma}} and \code{\link{dIndependent_Normal_Gamma}}
#' both place priors on \emph{both} the coefficients (joint normal) and the
#' dispersion (inverse-gamma), and \code{simulate_prior.glmb()} simulates
#' both parts for either. For \code{\link{dIndependent_Normal_Gamma}},
#' however, the dispersion prior is \emph{two-sided truncated} to
#' \code{x$pfamily$prior_list$disp_lower}/\code{disp_upper} -- bounds that
#' \code{\link{rindepNormalGamma_reg}} (and \code{\link{rglmb}}/
#' \code{\link{rlmb}}/\code{\link{glmb}}/\code{\link{lmb}}, which call it)
#' always compute and write back into the fitted object's \code{pfamily}
#' even when the user leaves \code{disp_lower}/\code{disp_upper} at their
#' default \code{NULL} at prior-specification time, because the
#' accept-reject envelope requires a bounded dispersion domain to construct
#' (lack of conjugacy between the independent coefficient and precision
#' priors rules out the closed-form, unbounded posterior available to
#' \code{\link{dNormal_Gamma}}). \code{simulate_prior.glmb()} draws from
#' this same truncated Inverse-Gamma via the exact inverse-CDF method
#' \file{src/invgamma_ct.cpp} uses internally, so it matches what the fit
#' actually assumed rather than the untruncated Inverse-Gamma(\code{shape},
#' \code{rate}). \code{\link{dNormal_Gamma}}'s dispersion prior has no
#' truncation concept at all (its \code{prior_list} has no
#' \code{disp_lower}/\code{disp_upper} fields; it is fully conjugate), and
#' \code{\link{dGamma}}'s truncation is opt-in only (via explicit
#' \code{disp_lower}/\code{disp_upper} arguments at construction; the
#' default \code{NULL} is never overwritten after fitting) -- both are
#' handled by the same helper, which falls back to the plain untruncated
#' Inverse-Gamma whenever \code{disp_lower}/\code{disp_upper} are \code{NULL}.
#'
#' \code{check_prior.glmb()} compares those prior draws against the fit's
#' posterior draws (\code{\link[insight]{get_parameters}}), one parameter at a
#' time, using the same \code{"gelman"} (posterior SD vs. prior SD ratio) and
#' \code{"lakeland"} (fraction of posterior mass inside the prior's 95\% HDI)
#' rules \pkg{bayestestR} uses for other model classes -- reimplemented
#' locally so \code{glmbayes} does not depend on \pkg{bayestestR}'s unexported
#' internal helpers. Because both rules only ever compare one parameter's
#' prior draws to that same parameter's posterior draws, the joint/marginal
#' distinction above never changes \code{check_prior()}'s verdicts; it only
#' matters if you use \code{simulate_prior()}'s draws directly (e.g. for
#' prior-predictive visualization).
#'
#' \code{describe_prior.glmb()} does \emph{not} route through
#' \code{get_priors()}'s flattened, marginal-only table (unlike
#' \pkg{bayestestR}'s own \code{describe_prior.stanreg()}, which is a thin
#' wrapper around \code{\link[insight]{get_priors}}). Instead it returns
#' \code{\link{pfamily}(model)} directly, so it prints the exact same
#' \code{Call} / \code{Prior Family} / \code{Prior List} report as
#' \code{pfamily(model)} (via the package's existing
#' \code{print.pfamily()} method) -- including the complete covariance
#' matrix for a multivariate normal coefficient prior, not just its
#' diagonal.
#'
#' @param model An object of class \code{"glmb"}.
#' @param n Number of prior draws to simulate.
#' @param method Either \code{"gelman"} or \code{"lakeland"}; see
#'   \code{\link[bayestestR]{check_prior}} for the general rule definitions.
#' @param simulate_priors Ignored (prior draws are always simulated via
#'   \code{simulate_prior.glmb()}); included for signature compatibility with
#'   \code{\link[bayestestR]{check_prior}}. \code{glmbayes} fits have no
#'   \code{unupdate()}-style refit mechanism (unlike \code{stanreg}/\code{brmsfit}),
#'   so \code{simulate_priors = FALSE} is not supported.
#' @param parameters Not used; included for signature compatibility with
#'   \code{\link[bayestestR]{describe_prior}}.
#' @param ... Not used; included for S3 signature compatibility with the
#'   \pkg{bayestestR} generics.
#' @return
#' \code{simulate_prior.glmb()} returns a \code{\link{data.frame}} with
#' \code{n} rows and one column per parameter (matching
#' \code{\link[insight]{get_parameters}}'s shape).
#'
#' \code{check_prior.glmb()} returns a \code{\link{data.frame}} with columns
#' \code{Parameter} and \code{Prior_Quality}.
#'
#' \code{describe_prior.glmb()} returns the fit's \code{"pfamily"} object
#' (see \code{\link{pfamily}}).
#' @seealso \code{\link{glmbayes_insight_methods}}, \code{\link{pfamily}},
#'   \code{\link{prior_simfuncs}}; \code{\link[bayestestR]{simulate_prior}},
#'   \code{\link[bayestestR]{check_prior}}, \code{\link[bayestestR]{describe_prior}}.
#' @example inst/examples/Ex_glmbayes_bayestestR_prior_methods.R
#' @name glmbayes_bayestestR_prior_methods
#' @aliases simulate_prior.glmb check_prior.glmb describe_prior.glmb
NULL

#' @rdname glmbayes_bayestestR_prior_methods
#' @method simulate_prior glmb
#' @export
simulate_prior.glmb <- function(model, n = 1000, ...) {
  pf <- model$pfamily
  if (is.null(pf) || is.null(pf$pfun)) {
    stop("simulate_prior() requires model$pfamily$pfun (see ?pfamily).", call. = FALSE)
  }

  params <- if (.glmbayes_is_dispersion_only_prior(model)) {
    NULL
  } else {
    colnames(model$coefficients)
  }

  pf$pfun(n = n, prior_list = pf$prior_list, params = params, ...)
}

## Column-by-column "gelman"/"lakeland" prior-informativeness rules, matching
## bayestestR's own (unexported) .check_prior() helper; reimplemented here
## rather than called via bayestestR::: to avoid depending on unexported API.
.glmbayes_check_prior_gelman <- function(prior, posterior) {
  if (all(is.na(prior))) {
    return("not determinable")
  }
  prior_sd <- stats::sd(prior, na.rm = TRUE)
  if (isTRUE(prior_sd == 0)) {
    ## Defensive fallback for a fixed/point-mass prior column (none of the
    ## five pfamily constructors currently produce one in simulate_prior.glmb()
    ## -- dGamma(Inv_Dispersion = TRUE)'s fixed beta is excluded entirely,
    ## see above -- but a zero-SD prior has an undefined posterior-SD/prior-SD
    ## ratio regardless of how it arose). Such a prior trivially and maximally
    ## determines the parameter, so it is reported as "informative" rather
    ## than falling through to the (here always FALSE)
    ## sd(posterior) > 0.1 * 0 comparison.
    return("informative")
  }
  if (stats::sd(posterior, na.rm = TRUE) > 0.1 * prior_sd) {
    "informative"
  } else {
    "uninformative"
  }
}

.glmbayes_check_prior_lakeland <- function(prior, posterior) {
  if (all(is.na(prior))) {
    return("not determinable")
  }
  hdi95 <- bayestestR::hdi(prior, ci = 0.95)
  r <- bayestestR::rope(posterior, ci = 1, range = c(hdi95$CI_low, hdi95$CI_high))
  if (as.numeric(r) > 0.99) "informative" else "misinformative"
}

#' @rdname glmbayes_bayestestR_prior_methods
#' @method check_prior glmb
#' @export
check_prior.glmb <- function(model, method = "gelman", simulate_priors = TRUE, ...) {
  method <- match.arg(method, c("gelman", "lakeland"))
  posteriors <- insight::get_parameters(model)
  priors <- simulate_prior(model, n = max(1000L, nrow(posteriors)))

  common <- intersect(colnames(priors), colnames(posteriors))
  check_fn <- if (method == "gelman") .glmbayes_check_prior_gelman else .glmbayes_check_prior_lakeland
  result <- mapply(check_fn, priors[common], posteriors[common])

  data.frame(Parameter = common, Prior_Quality = unname(result), stringsAsFactors = FALSE)
}

#' @rdname glmbayes_bayestestR_prior_methods
#' @method describe_prior glmb
#' @export
describe_prior.glmb <- function(model, parameters = NULL, ...) {
  pfamily(model)
}
