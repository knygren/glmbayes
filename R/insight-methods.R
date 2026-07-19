#' insight Accessor Methods for glmb/lmb Model Objects
#'
#' @description
#' Methods implementing the \pkg{insight} accessor generics for objects of
#' class \code{"glmb"} (which includes \code{\link{lmb}}, \code{\link{rglmb}},
#' \code{\link{rlmb}}, and \code{\link{rGamma_reg}} fits, since these all
#' inherit \code{"glmb"} in their class vector). This lets \pkg{insight}-based
#' tooling (including \pkg{bayestestR}) recognize \code{glmbayes} fits
#' as Bayesian models and extract posterior draws in the standard
#' \pkg{insight} shape (one row per draw, one column per parameter).
#'
#' \code{model_info.glmb()} always reports \code{is_bayesian = TRUE} and
#' derives family/link-specific flags (\code{is_linear}, \code{is_count},
#' \code{is_poisson}, \code{is_binomial}, \code{is_trial}, \code{is_logit},
#' \code{is_probit}, \code{is_exponential}) from \code{x$family} (defaulting
#' to \code{\link[stats]{gaussian}} when absent, which is the case for
#' \code{\link{lmb}} fits). All other boolean fields reported by
#' \code{\link[insight]{model_info}} (e.g. \code{is_mixed},
#' \code{is_zero_inflated}, \code{is_ordinal}, \code{is_gam},
#' \code{is_multivariate}) are explicitly set to \code{FALSE} rather than
#' omitted, since some \pkg{insight}/\pkg{bayestestR} internals test these
#' fields with \code{if()} and error on \code{NULL}.
#'
#' \code{get_parameters.glmb()} returns \code{x$coefficients} (the posterior
#' draws matrix) as a data frame, with a \code{dispersion} column appended
#' only when the dispersion parameter was actually sampled (i.e. varies
#' across draws) rather than held fixed. \code{\link{dGamma}} with
#' \code{Inv_Dispersion = TRUE} (an \code{"rGamma_reg"} fit) is a special
#' case: that prior is on the inverse dispersion (precision) \emph{only}
#' -- \code{beta} is a fixed, known input (typically held fixed in a Gibbs
#' step), not itself a modeled/estimated quantity -- so
#' \code{get_parameters.glmb()} reports \emph{only} a \code{dispersion}
#' column and omits the coefficients entirely, rather than recycling the
#' fixed row as if it were a (degenerate) posterior draw.
#' \code{\link{dGamma}} with \code{Inv_Dispersion = FALSE} is the mirror
#' image: there, the prior \emph{is} on the coefficient (the conjugate
#' Gamma-Poisson/Gamma-Gamma rate), which is genuinely sampled, and
#' dispersion is fixed -- so that case reports coefficients as usual and
#' omits \code{dispersion}.
#'
#' \code{find_parameters.glmb()} returns the same parameter grouping as a
#' named list of parameter \emph{names} (not values): \code{conditional} for
#' the regression coefficients, plus \code{dispersion} only when sampled --
#' or, for \code{dGamma(Inv_Dispersion = TRUE)}, \code{dispersion} only,
#' with no \code{conditional} entry at all (mirroring
#' \code{get_parameters.glmb()} above).
#'
#' \code{find_algorithm.glmb()} reports the estimation method as
#' \code{"iid rejection sampling"} rather than falling through to the
#' inherited \code{"glm"}/\code{"lm"} methods (which would otherwise report
#' \code{"ML"}/\code{"OLS"} — incorrect for a Bayesian iid-sampling fit).
#' There is no \code{chains}/\code{warmup} concept (draws are independent,
#' not an MCMC chain), so only \code{algorithm} and \code{iterations} are
#' reported, mirroring \code{\link[insight]{find_algorithm}}'s treatment of
#' other non-MCMC Bayesian samplers (e.g. \code{bayesQR}).
#'
#' \code{get_data.glmb()} returns \code{x$data} (or \code{x$model}, the model
#' frame, when \code{x$data} is not itself a data frame -- e.g. when
#' \code{glmb()}/\code{lmb()} was called without an explicit \code{data}
#' argument, in which case \code{x$data} is the calling environment, as for
#' \code{\link[stats]{glm}}). This avoids a spurious \code{"Could not recover
#' model data from environment"} warning that \code{\link[insight]{get_data}}'s
#' default method emits for \code{glmb} fits before falling back to
#' (successfully) reconstructing the same data from the model frame.
#'
#' \code{get_priors.glmb()} does \emph{not} translate the \code{pfamily}
#' prior specification into \code{insight}'s usual
#' \code{Parameter}/\code{Distribution}/\code{Location}/\code{Scale} shape
#' (unlike, say, \code{get_priors.stanreg}/\code{get_priors.brmsfit}). That
#' shape has one row per parameter with only a scalar
#' \code{Location}/\code{Scale}, so it cannot represent off-diagonal
#' covariance for a genuinely multivariate normal coefficient prior.
#' Instead, \code{get_priors.glmb()} returns \code{\link{pfamily}(x)}
#' directly, so it prints the exact same \code{Call}/\code{Prior
#' Family}/\code{Prior List} report as \code{pfamily(x)} (via the package's
#' existing \code{print.pfamily()} method) -- including the complete
#' covariance matrix, not just its diagonal. See
#' \code{\link{glmbayes_bayestestR_prior_methods}}'s
#' \code{simulate_prior.glmb()} for a consumer that draws from this true
#' joint distribution, and \code{describe_prior.glmb()}, which returns the
#' same thing under \pkg{bayestestR}'s naming.
#'
#' @param x An object of class \code{"glmb"}, typically returned by
#'   \code{\link{glmb}}, \code{\link{lmb}}, \code{\link{rglmb}}, or
#'   \code{\link{rlmb}}.
#' @param verbose Not used; included for S3 signature compatibility with
#'   \code{\link[insight]{get_priors}}.
#' @param ... Not used; included for S3 signature compatibility with the
#'   \pkg{insight} generics.
#' @return
#' \code{model_info.glmb()} returns a named \code{\link{list}} of model
#' metadata matching the shape of \code{\link[insight]{model_info}}'s output
#' for other model classes.
#'
#' \code{get_parameters.glmb()} returns a \code{\link{data.frame}} with one
#' row per posterior draw and one column per parameter (plus a
#' \code{dispersion} column when sampled).
#'
#' \code{find_parameters.glmb()} returns a named \code{\link{list}} of
#' character vectors (parameter names), with components \code{conditional}
#' and (when sampled) \code{dispersion}.
#'
#' \code{find_algorithm.glmb()} returns a named \code{\link{list}} with
#' components \code{algorithm} and \code{iterations}.
#'
#' \code{get_data.glmb()} returns a \code{\link{data.frame}} (the model data).
#'
#' \code{get_priors.glmb()} returns the fit's \code{"pfamily"} object (see
#' \code{\link{pfamily}}).
#' @seealso \code{\link{glmb}}, \code{\link{lmb}}, \code{\link{rglmb}},
#'   \code{\link{rlmb}}, \code{\link{glmbayes-package}}, \code{\link{pfamily}};
#'   \code{\link[insight]{model_info}}, \code{\link[insight]{get_parameters}},
#'   \code{\link[insight]{find_parameters}}, \code{\link[insight]{find_algorithm}},
#'   \code{\link[insight]{get_data}}, \code{\link[insight]{get_priors}}.
#' @example inst/examples/Ex_glmbayes_insight_methods.R
#' @name glmbayes_insight_methods
#' @aliases model_info.glmb get_parameters.glmb find_parameters.glmb
#'   find_algorithm.glmb get_data.glmb get_priors.glmb
NULL

## Full field template mirroring insight::model_info()'s output shape (checked
## against installed insight 1.5.1 for stats::lm/glm objects); defaults to
## FALSE/is_bayesian=TRUE so every field is present with a concrete value,
## since some insight/bayestestR internals test these with if() and error on
## NULL rather than treating a missing field as falsy.
.glmbayes_model_info_template <- function() {
  list(
    is_binomial = FALSE, is_bernoulli = FALSE, is_count = FALSE, is_poisson = FALSE,
    is_negbin = FALSE, is_beta = FALSE, is_betabinomial = FALSE, is_orderedbeta = FALSE,
    is_dirichlet = FALSE, is_exponential = FALSE, is_logit = FALSE, is_probit = FALSE,
    is_censored = FALSE, is_truncated = FALSE, is_survival = FALSE, is_linear = FALSE,
    is_tweedie = FALSE, is_zeroinf = FALSE, is_zero_inflated = FALSE, is_dispersion = FALSE,
    is_hurdle = FALSE, is_ordinal = FALSE, is_cumulative = FALSE, is_multinomial = FALSE,
    is_categorical = FALSE, is_mixed = FALSE, is_multivariate = FALSE, is_trial = FALSE,
    is_bayesian = TRUE, is_gam = FALSE, is_anova = FALSE, is_timeseries = FALSE,
    is_ttest = FALSE, is_correlation = FALSE, is_onewaytest = FALSE, is_chi2test = FALSE,
    is_ranktest = FALSE, is_levenetest = FALSE, is_variancetest = FALSE, is_ztest = FALSE,
    is_xtab = FALSE, is_proptest = FALSE, is_binomtest = FALSE, is_ftest = FALSE,
    is_meta = FALSE, is_wiener = FALSE, is_rtchoice = FALSE, is_mixture = FALSE
  )
}

## TRUE only when the dispersion parameter was actually sampled (varies across
## draws) rather than held fixed at a constant value; see data-raw/
## glmbayes-bayestestR-integration.md for the underlying rationale.
.glmbayes_dispersion_sampled <- function(x) {
  disp <- x$dispersion
  !is.null(disp) && length(unique(disp)) > 1L
}

## TRUE only for dGamma(Inv_Dispersion = TRUE) fits ("rGamma_reg"), where the
## prior is on the inverse dispersion (precision) *only*: beta is a fixed,
## known input (typically held fixed in a Gibbs step), not itself a
## modeled/estimated quantity. Shared by get_parameters.glmb(),
## find_parameters.glmb(), and simulate_prior.glmb() (R/bayestestR-methods.R)
## so all three consistently omit the coefficients rather than treating the
## fixed beta as if it were a (degenerate) posterior/prior draw.
.glmbayes_is_dispersion_only_prior <- function(x) {
  pf <- x$pfamily
  identical(pf$pfamily, "dGamma") && isTRUE(pf$prior_list$Inv_Dispersion)
}

#' @rdname glmbayes_insight_methods
#' @method model_info glmb
#' @export
model_info.glmb <- function(x, ...) {
  family   <- if (!is.null(x$family)) x$family else stats::gaussian()
  fam_name  <- family$family
  link_name <- family$link

  out <- .glmbayes_model_info_template()

  if (identical(fam_name, "gaussian")) {
    out$is_linear <- TRUE
  } else if (fam_name %in% c("poisson", "quasipoisson")) {
    out$is_count   <- TRUE
    out$is_poisson <- TRUE
  } else if (fam_name %in% c("binomial", "quasibinomial")) {
    out$is_binomial <- TRUE
    out$is_trial    <- TRUE
    if (identical(link_name, "logit"))  out$is_logit  <- TRUE
    if (identical(link_name, "probit")) out$is_probit <- TRUE
  } else if (identical(fam_name, "Gamma")) {
    out$is_exponential <- TRUE
  }

  out$link_function <- link_name
  out$family         <- fam_name
  out$n_obs          <- NROW(x$y)
  out$n_grouplevels   <- NULL

  out
}

#' @rdname glmbayes_insight_methods
#' @method get_parameters glmb
#' @export
get_parameters.glmb <- function(x, ...) {
  ## dGamma(Inv_Dispersion = TRUE): beta is a fixed, known input, not a
  ## modeled quantity -- report only the (actually sampled) dispersion,
  ## and omit the coefficients entirely.
  if (.glmbayes_is_dispersion_only_prior(x)) {
    return(data.frame(dispersion = x$dispersion))
  }

  out <- as.data.frame(x$coefficients)

  if (.glmbayes_dispersion_sampled(x)) {
    out$dispersion <- x$dispersion
  }

  out
}

#' @rdname glmbayes_insight_methods
#' @method find_parameters glmb
#' @export
find_parameters.glmb <- function(x, ...) {
  if (.glmbayes_is_dispersion_only_prior(x)) {
    return(list(dispersion = "dispersion"))
  }

  out <- list(conditional = colnames(x$coefficients))

  if (.glmbayes_dispersion_sampled(x)) {
    out$dispersion <- "dispersion"
  }

  out
}

#' @rdname glmbayes_insight_methods
#' @method find_algorithm glmb
#' @export
find_algorithm.glmb <- function(x, ...) {
  list(algorithm = "iid rejection sampling", iterations = nrow(x$coefficients))
}

#' @rdname glmbayes_insight_methods
#' @method get_data glmb
#' @export
get_data.glmb <- function(x, ...) {
  ## x$data is only a proper data frame when glmb()/lmb() was called with an
  ## explicit data= argument; otherwise it is the calling environment (as for
  ## stats::glm()). x$model (the model frame) is always a data frame.
  if (is.data.frame(x$data)) return(x$data)
  if (is.data.frame(x$model)) return(x$model)
  NextMethod()
}

#' @rdname glmbayes_insight_methods
#' @method get_priors glmb
#' @export
get_priors.glmb <- function(x, verbose = TRUE, ...) {
  pfamily(x)
}
