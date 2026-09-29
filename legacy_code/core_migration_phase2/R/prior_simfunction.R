#' Prior Simulation Functions for glmbayes pfamily Objects
#'
#' @description
#' Prior-simulation functions provide a unified interface for generating iid
#' draws directly from the \emph{prior} distribution stored in a
#' \code{\link{pfamily}} object's \code{prior_list}, mirroring the existing
#' \code{simfuncs} (\code{\link{rNormal_reg}}, \code{\link{rGamma_reg}},
#' \code{\link{rGamma_Conjugate_reg}}, \code{\link{rNormalGamma_reg}},
#' \code{\link{rindepNormalGamma_reg}}, \code{\link{rBeta_reg}}) that already
#' generate iid draws from the \emph{posterior}.
#'
#' Every \code{pfamily} constructor (\code{\link{dNormal}}, \code{\link{dGamma}},
#' \code{\link{dNormal_Gamma}}, \code{\link{dIndependent_Normal_Gamma}},
#' \code{\link{dBeta}}) has a corresponding prior-simulation function below,
#' named by dropping the leading \code{"d"} from the constructor's name,
#' prepending \code{"r"}, and appending \code{"_prior"} (e.g.
#' \code{dNormal()} -> \code{rNormal_prior()}), matching the existing
#' \code{simfuncs}' own \code{"r..._reg"} naming. \code{dGamma()} maps to one
#' of two functions depending on \code{Inv_Dispersion}, exactly as it already
#' selects between \code{\link{rGamma_reg}} and
#' \code{\link{rGamma_Conjugate_reg}} for its \code{simfun}:
#' \code{rGamma_prior()} (\code{Inv_Dispersion = TRUE}: prior on the inverse
#' dispersion only) and \code{rGamma_Conjugate_prior()}
#' (\code{Inv_Dispersion = FALSE}: conjugate prior on the rate directly).
#'
#' All six functions share an \strong{identical signature}
#' \code{function(n, prior_list, params = NULL, ...)}, so a caller holding
#' only \code{n} and a \code{prior_list} (plus, optionally, a vector of
#' parameter names) can invoke any of them the same way -- exactly as
#' \code{\link{rglmb}}/\code{\link{rlmb}} call whichever \code{simfun} a
#' \code{pfamily} object carries without needing to know which one it is.
#' Each \code{\link{pfamily}} constructor stores the matching function as
#' \code{pfun}, alongside \code{simfun}, so
#' \code{\link[bayestestR]{simulate_prior}} methods can extract the
#' \code{pfamily} and call \code{pfun(n, prior_list, params)} directly --
#' mirroring how \code{\link{rglmb}}/\code{\link{rlmb}} call \code{simfun}.
#'
#' Coefficient draws for a genuinely multivariate normal prior
#' (\code{\link{dNormal}}, \code{\link{dNormal_Gamma}},
#' \code{\link{dIndependent_Normal_Gamma}}) are drawn from the \emph{true}
#' joint \eqn{N(\mu, \Sigma)} via a Cholesky factor of the complete
#' (possibly non-diagonal) \code{Sigma}, preserving any correlation between
#' coefficients -- not reconstructed from per-parameter marginal moments.
#' \code{\link{dBeta}}/\code{\link{dGamma}} (\code{Inv_Dispersion = FALSE})
#' coefficients are drawn independently, one \code{\link[stats]{rbeta}}/
#' \code{\link[stats]{rgamma}} draw per parameter, from their exact
#' \code{shape1}/\code{shape2} or \code{shape}/\code{rate} form, rather than
#' their normal-moment surrogate.
#'
#' \code{\link{dNormal_Gamma}} and \code{\link{dIndependent_Normal_Gamma}}
#' both place priors on \emph{both} the coefficients (joint normal) and the
#' dispersion (inverse-gamma), so \code{rNormal_Gamma_prior()} and
#' \code{rIndependent_Normal_Gamma_prior()} return both a coefficient block
#' and a \code{dispersion} column. For \code{\link{dIndependent_Normal_Gamma}}
#' the dispersion prior is two-sided truncated to
#' \code{prior_list$disp_lower}/\code{disp_upper} (see
#' \code{\link{rindepNormalGamma_reg}} and the \code{dIndependent_Normal_Gamma}
#' branches in \code{\link{rglmb}}/\code{\link{rlmb}}/\code{\link{glmb}}/
#' \code{\link{lmb}}, which always populate these -- even when left at their
#' default \code{NULL} at prior-specification time -- with the bounds the
#' accept-reject envelope actually used); \code{rIndependent_Normal_Gamma_prior()}
#' samples that exact truncated Inverse-Gamma via
#' \code{.glmbayes_rinvgamma_prior()} (\file{R/prior_simfunction.R}), the
#' same inverse-CDF method \file{src/invgamma_ct.cpp} uses internally.
#' \code{\link{dNormal_Gamma}}'s dispersion prior has no truncation concept
#' at all (fully conjugate; \code{disp_lower}/\code{disp_upper} are absent
#' from its \code{prior_list}), and \code{\link{dGamma}}'s truncation is
#' opt-in only (via explicit \code{disp_lower}/\code{disp_upper} arguments;
#' the default \code{NULL} is never overwritten after fitting) -- both fall
#' back to the plain untruncated Inverse-Gamma via the same helper.
#'
#' @param n Number of prior draws to generate.
#' @param prior_list A list with prior parameters (\code{mu}, \code{Sigma},
#'   \code{shape}, \code{rate}, \code{shape1}, \code{shape2}, \code{beta},
#'   \code{disp_lower}, \code{disp_upper}, etc., as relevant), of the same
#'   form stored in a \code{\link{pfamily}} object's \code{prior_list}.
#' @param params Optional character vector of coefficient names for the
#'   output columns. Defaults to \code{rownames(prior_list$mu)} when
#'   \code{NULL} (as set by \code{\link{Prior_Setup}} and the
#'   \code{\link{pfamily}} constructors when given named input); falls back
#'   to \code{"V1"}, \code{"V2"}, ... if neither is available.
#' @param \ldots Additional arguments; currently unused, present so all six
#'   functions share an identical signature and so extra arguments passed by
#'   a generic caller are silently ignored.
#'
#' @return
#' A \code{\link{data.frame}} with \code{n} rows. \code{rGamma_prior()}
#' (\code{dGamma(Inv_Dispersion = TRUE)}) returns a single \code{dispersion}
#' column and no coefficient columns, since that prior is on the inverse
#' dispersion only. The other five return one column per coefficient, named
#' from \code{params}; \code{rNormal_Gamma_prior()} and
#' \code{rIndependent_Normal_Gamma_prior()} additionally append a
#' \code{dispersion} column.
#'
#' @seealso \code{\link{pfamily}}, \code{\link{dNormal}}, \code{\link{dGamma}},
#'   \code{\link{dNormal_Gamma}}, \code{\link{dIndependent_Normal_Gamma}},
#'   \code{\link{dBeta}}; \code{\link{simfuncs}} for the analogous
#'   posterior-simulation functions; \code{\link{glmbayes_bayestestR_prior_methods}}
#'   for \code{\link[bayestestR]{simulate_prior}} on fitted \code{glmb} objects.
#' @name prior_simfuncs
#' @aliases rNormal_prior rGamma_prior rGamma_Conjugate_prior
#'   rNormal_Gamma_prior rIndependent_Normal_Gamma_prior rBeta_prior
NULL

## Exact inverse-CDF sampler for the Inverse-Gamma(shape, rate) dispersion
## prior, mirroring src/invgamma_ct.cpp. When lower/upper are both NULL this
## is the plain (untruncated) Inverse-Gamma; dIndependent_Normal_Gamma() fits
## always carry concrete disp_lower/disp_upper in prior_list after sampling.
.glmbayes_rinvgamma_prior <- function(n, shape, rate, lower = NULL, upper = NULL) {
  if (is.null(lower) || is.null(upper)) {
    return(1 / stats::rgamma(n, shape = shape, rate = rate))
  }
  cdf <- function(d) 1 - stats::pgamma(1 / d, shape = shape, rate = rate)
  p_low <- cdf(lower)
  p_upp <- cdf(upper)
  p1 <- p_low + stats::runif(n) * (p_upp - p_low)
  1 / stats::qgamma(1 - p1, shape = shape, rate = rate)
}

## Shared helper: n draws from the joint N(mu, Sigma) via a Cholesky factor
## of the complete (possibly non-diagonal) Sigma, so any correlation between
## coefficients is preserved (see prior_simfuncs). Returns an n x length(mu)
## data.frame with column names `params` (or "V1", "V2", ... if params and
## rownames(mu) are both NULL).
.glmbayes_prior_mvn_draws <- function(n, mu, Sigma, params = NULL) {
  mu_vec <- as.numeric(mu)
  p <- length(mu_vec)
  nm <- params
  if (is.null(nm)) nm <- rownames(as.matrix(mu, ncol = 1L))
  if (is.null(nm)) nm <- paste0("V", seq_len(p))

  R <- chol(as.matrix(Sigma))     ## upper-triangular, t(R) %*% R == Sigma
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  draws <- sweep(z %*% R, 2, mu_vec, "+")
  colnames(draws) <- nm
  as.data.frame(draws)
}

## Shared helper: n independent draws per parameter for priors whose
## coefficients are each marginally Beta(shape1, shape2)- or
## Gamma(shape, rate)-distributed (dBeta(), dGamma(Inv_Dispersion = FALSE)),
## rather than jointly normal. `rfun` is stats::rbeta or stats::rgamma;
## `...` supplies its shape parameters.
.glmbayes_prior_indep_draws <- function(rfun, n, p, params = NULL, ...) {
  nm <- params
  if (is.null(nm)) nm <- paste0("V", seq_len(p))
  draws <- matrix(rfun(n * p, ...), nrow = n, ncol = p)
  colnames(draws) <- nm
  as.data.frame(draws)
}

#' @rdname prior_simfuncs
#' @export
rNormal_prior <- function(n, prior_list, params = NULL, ...) {
  .glmbayes_prior_mvn_draws(n, prior_list$mu, prior_list$Sigma, params)
}

#' @rdname prior_simfuncs
#' @export
rGamma_prior <- function(n, prior_list, params = NULL, ...) {
  data.frame(
    dispersion = .glmbayes_rinvgamma_prior(
      n, prior_list$shape, prior_list$rate,
      prior_list$disp_lower, prior_list$disp_upper
    )
  )
}

#' @rdname prior_simfuncs
#' @export
rGamma_Conjugate_prior <- function(n, prior_list, params = NULL, ...) {
  p <- length(as.numeric(prior_list$mu))
  nm <- params
  if (is.null(nm)) nm <- rownames(as.matrix(prior_list$mu, ncol = 1L))
  .glmbayes_prior_indep_draws(
    stats::rgamma, n, p, nm,
    shape = prior_list$shape, rate = prior_list$rate
  )
}

#' @rdname prior_simfuncs
#' @export
rNormal_Gamma_prior <- function(n, prior_list, params = NULL, ...) {
  out <- .glmbayes_prior_mvn_draws(n, prior_list$mu, prior_list$Sigma, params)
  out$dispersion <- .glmbayes_rinvgamma_prior(n, prior_list$shape, prior_list$rate)
  out
}

#' @rdname prior_simfuncs
#' @export
rIndependent_Normal_Gamma_prior <- function(n, prior_list, params = NULL, ...) {
  out <- .glmbayes_prior_mvn_draws(n, prior_list$mu, prior_list$Sigma, params)
  out$dispersion <- .glmbayes_rinvgamma_prior(
    n, prior_list$shape, prior_list$rate,
    prior_list$disp_lower, prior_list$disp_upper
  )
  out
}

#' @rdname prior_simfuncs
#' @export
rBeta_prior <- function(n, prior_list, params = NULL, ...) {
  p <- length(as.numeric(prior_list$mu))
  nm <- params
  if (is.null(nm)) nm <- rownames(as.matrix(prior_list$mu, ncol = 1L))
  .glmbayes_prior_indep_draws(
    stats::rbeta, n, p, nm,
    shape1 = prior_list$shape1, shape2 = prior_list$shape2
  )
}
