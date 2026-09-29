#' Compute Calibrated Gaussian Normal–Gamma Prior Components
#'
#' Internal Gaussian calibration routine used by \code{\link{Prior_Setup}}.
#' Given weighted Gaussian regression inputs and a dispersion–independent
#' coefficient–scale prior covariance \eqn{\Sigma_0}, this function computes
#' all Normal–Gamma quantities required by the Gaussian prior families:
#' \code{dispersion}, \code{shape}, \code{shape_ING}, \code{rate},
#' \code{rate_gamma}, and the calibrated coefficient–scale covariance
#' \code{Sigma}. The input \code{Sigma_0} is returned unchanged.
#'
#' The computation follows a structured pipeline:
#' \enumerate{
#'   \item Validate dimensions and numeric inputs.
#'   \item Compute the weighted residual sum of squares at \code{bhat}.
#'   \item Form the weighted Gram matrix \eqn{X^\top W X}, invert it, and
#'         construct the marginal quadratic term \eqn{S_{marg}}.
#'   \item Map \code{n_prior} and \code{k} to the Normal–Gamma shape and rate.
#'   \item Calibrate the implied dispersion and coefficient covariance.
#'   \item Return all calibrated prior components.
#' }
#'
#' @details
#' The function assumes the Chapter 11 convention that \eqn{\Sigma_0} is
#' *dispersion–free*: it encodes prior structure on the precision–weighted
#' coefficient scale. The returned \code{Sigma} is the corresponding
#' coefficient–scale covariance after calibration.
#'
#' A common choice is the Zellner–type form
#' \deqn{
#'   \Sigma_0 = \frac{1 - \mathrm{pwt}}{\mathrm{pwt}} (X^\top W X)^{-1},
#' }
#' where \eqn{\mathrm{pwt}} is a scalar prior weight. More generally,
#' \code{Sigma_0} may be any positive–definite matrix.
#'
#' The function computes:
#' \itemize{
#'   \item The marginal quadratic term
#'     \deqn{
#'       S_{marg}
#'       = RSS_w
#'       + (\hat\beta - \mu)^\top
#'         \left(\Sigma_0 + (X^\top W X)^{-1}\right)^{-1}
#'         (\hat\beta - \mu).
#'     }
#'
#'   \item Prior Gamma shape:
#'     \eqn{a_0 = (n_{\mathrm{prior}} + k)/2}.
#'
#'   \item Posterior Gamma shape:
#'     \eqn{a_n = (n_{\mathrm{prior}} + n_w + k)/2},
#'     where \eqn{n_w = n_{\mathrm{effective}}}.
#'
#'   \item Calibrated dispersion:
#'     \deqn{
#'       E[\sigma^2 \mid y]
#'       = \frac{S_{marg}}{n_w - p},
#'     }
#'     the usual weighted residual–df estimator.
#'
#'   \item Prior Gamma rate:
#'     \deqn{
#'       b_0
#'       = \frac{1}{2} S_{marg}
#'         \frac{n_{\mathrm{prior}} + k + p - 2}{n_w - p},
#'     }
#'     ensuring \eqn{E[\sigma^2 \mid y] = S_{marg}/(n_w - p)}.
#'
#'   \item Calibrated coefficient covariance:
#'     \deqn{
#'       \Sigma
#'       = \frac{n_w}{n_{\mathrm{prior}}}
#'         E[\sigma^2 \mid y]
#'         (X^\top W X)^{-1}.
#'     }
#' }
#'
#' \strong{Limiting behavior.}
#' \itemize{
#'   \item As \eqn{n_{\mathrm{prior}} \to \infty}, the prior becomes increasingly
#'         concentrated and dominates the likelihood.
#'
#'   \item As \eqn{n_{\mathrm{prior}} \to 0^+} with \eqn{n_w > p},
#'         \eqn{S_{marg} \to RSS_w} and
#'         \eqn{E[\sigma^2 \mid y] \to RSS_w/(n_w - p)},
#'         matching the classical weighted Gaussian estimator.
#'
#'   \item Strict positivity of the prior rate requires
#'         \eqn{n_{\mathrm{prior}} + k + p > 2}.
#'
#'   \item The prior mean \eqn{\mu} is never altered; this function calibrates
#'         only scale parameters.
#' }
#'
#' @param X Numeric model matrix with \code{nrow(X) == length(Y)}.
#' @param Y Numeric response vector.
#' @param weights Numeric vector of case weights.
#' @param offset Numeric offset vector.
#' @param dispersion Optional scalar dispersion. If supplied, overrides the
#'   calibrated value.
#' @param n_effective Effective sample size (typically \code{sum(weights)}).
#' @param bhat Numeric coefficient vector, usually the weighted least–squares
#'   estimate.
#' @param mu Numeric prior mean vector.
#' @param Sigma_0 Dispersion–independent prior covariance matrix \eqn{[p \times p]}.
#' @param Sigma Optional coefficient–scale covariance matrix. If supplied,
#'   overrides the calibrated \code{Sigma}.
#' @param n_prior Effective prior sample size.
#' @param k Non–negative scalar with \eqn{k + p \ge 2}. Controls tail behavior
#'   of the variance prior; does not affect posterior means.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{dispersion} — calibrated Gaussian dispersion.
#'   \item \code{shape} — Gamma shape for residual precision.
#'   \item \code{shape_ING} — shape for the independent Normal–Gamma prior.
#'   \item \code{rate} — Gamma rate for residual precision.
#'   \item \code{rate_gamma} — Gamma rate for the fixed–\eqn{\beta} path
#'         (\code{\link{dGamma}}).
#'   \item \code{Sigma} — calibrated coefficient–scale covariance.
#'   \item \code{Sigma_0} — the input dispersion–free covariance.
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom Rdpack reprompt
#' @export
compute_gaussian_prior <- function(
    X,
    Y,
    weights,
    offset,
    dispersion = NULL,
    n_effective,
    bhat,
    mu,
    Sigma_0,
    Sigma = NULL,
    n_prior,
    k = 1
) {
  ## ---------------------------------------------------------------------------
  ## Gaussian calibration pipeline:
  ## Step A: Validate inputs and dimensions.
  ## Step B: Compute weighted RSS from (Y, X, bhat, offset, weights).
  ## Step C: Build Gram terms and S_marg using Sigma_0 + (X'WX)^{-1}.
  ## Step D: Gamma shape = (n_prior + k) / 2; rate from n_prior, k, p.
  ## Step E: Calibrate dispersion/rate and map to coefficient Sigma.
  ## Step F: Return calibrated terms.
  ## ---------------------------------------------------------------------------
  if (!is.null(dispersion)) {
    if (!is.numeric(dispersion) || length(dispersion) != 1L ||
        !is.finite(dispersion) || dispersion <= 0) {
      stop("compute_gaussian_prior: dispersion must be NULL or a single positive finite numeric value.", call. = FALSE)
    }
  }
  dispersion_input <- dispersion
  Sigma_input <- Sigma

  ## Step A: validate all required Gaussian inputs.
  n_obs <- NROW(Y)
  if (!is.matrix(X) || NROW(X) != n_obs) {
    stop("compute_gaussian_prior: X must be a matrix with nrow(X) == length(Y).", call. = FALSE)
  }
  if (!is.numeric(Y) || length(Y) != n_obs) {
    stop("compute_gaussian_prior: Y must be a numeric vector with length equal to nrow(X).", call. = FALSE)
  }
  if (!is.numeric(weights) || length(weights) != n_obs) {
    stop("compute_gaussian_prior: weights must be a numeric vector with length equal to nrow(X).", call. = FALSE)
  }
  if (!is.numeric(offset) || length(offset) != n_obs) {
    stop("compute_gaussian_prior: offset must be a numeric vector with length equal to nrow(X).", call. = FALSE)
  }
  p <- NCOL(X)
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 0) {
    stop("compute_gaussian_prior: k must be a single non-negative finite numeric value.", call. = FALSE)
  }
  if (k + p < 2) {
    stop(
      "compute_gaussian_prior: require k + p >= 2, where p = ncol(X). Got k = ", k, ", p = ", p, ".",
      call. = FALSE
    )
  }
  if (!is.numeric(bhat) || length(bhat) != p || any(!is.finite(bhat))) {
    stop("compute_gaussian_prior: bhat must be a finite numeric vector with length ncol(X).", call. = FALSE)
  }
  mu_num <- as.numeric(mu)
  if (length(mu_num) != p || any(!is.finite(mu_num))) {
    stop("compute_gaussian_prior: mu must be a finite numeric vector with length ncol(X).", call. = FALSE)
  }
  if (!is.matrix(Sigma_0) || nrow(Sigma_0) != p || ncol(Sigma_0) != p || anyNA(Sigma_0)) {
    stop("compute_gaussian_prior: Sigma_0 must be a numeric [p x p] matrix with no missing values.", call. = FALSE)
  }
  if (!is.numeric(n_prior) || length(n_prior) != 1L || !is.finite(n_prior) || n_prior <= 0) {
    stop("compute_gaussian_prior: n_prior must be a single positive finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(n_effective) || length(n_effective) != 1L || !is.finite(n_effective) || n_effective <= 0) {
    stop("compute_gaussian_prior: n_effective must be a single positive finite numeric value.", call. = FALSE)
  }

  ## Step B: weighted residual sum of squares at bhat.
  res <- as.numeric(Y) - as.numeric(X %*% bhat) - as.numeric(offset)
  rss_weighted <- sum(as.numeric(weights) * res^2)
  if (!is.finite(rss_weighted) || rss_weighted <= 0) {
    stop("compute_gaussian_prior: weighted RSS must be strictly positive.", call. = FALSE)
  }
  if (!is.finite(p) || p < 1L) {
    stop("compute_gaussian_prior: require ncol(X) >= 1.", call. = FALSE)
  }
  if (n_effective <= p) {
    stop(
      "compute_gaussian_prior: require n_effective > p (number of coefficients) for Gaussian dispersion (denominator n_effective - p). ",
      "Got n_effective = ", n_effective, " and p = ", p, ".",
      call. = FALSE
    )
  }

  ## Step C: weighted Gram inverse and S_marg quadratic augmentation.
  XtW <- sweep(X, 1, as.numeric(weights), `*`)
  Gm <- crossprod(XtW, X)
  Ginv <- tryCatch(
    solve(Gm),
    error = function(e) {
      stop("compute_gaussian_prior: cannot invert weighted Gram matrix X'WX. ", conditionMessage(e), call. = FALSE)
    }
  )
  dlt <- matrix(bhat, ncol = 1L) - matrix(mu_num, ncol = 1L)
  M <- Sigma_0 + Ginv
  Mi <- tryCatch(
    solve(M),
    error = function(e) {
      stop("compute_gaussian_prior: cannot invert Sigma_0 + (X'WX)^{-1}. ", conditionMessage(e), call. = FALSE)
    }
  )
  quad <- as.numeric(crossprod(dlt, Mi %*% dlt))
  if (!is.finite(quad) || quad < 0) {
    stop("compute_gaussian_prior: S_marg quadratic form is not finite or nonnegative.", call. = FALSE)
  }
  S_marg <- rss_weighted + quad

  ## Step D: prior Gamma shape (precision) uses n_prior and k.
  shape <- (n_prior + k) / 2
  if (!is.finite(shape) || shape <= 0) {
    stop("compute_gaussian_prior: computed shape must be strictly positive.", call. = FALSE)
  }

  ## Step E: calibrate Gaussian dispersion/rate and implied Sigma.
  ## Prior shape a_0 = (n_prior + k)/2; posterior a_n = a_0 + n_effective/2
  ##   = (n_prior + k + n_effective)/2. With b_n = b_0 + S_marg/2 and
  ##   b_0 = (S_marg/2) * (n_prior + k + p - 2) / (n_effective - p),
  ##   E[sigma^2|y] = b_n/(a_n - 1) = S_marg/(n_effective - p) = dispersion_cal
  ##   (k cancels between a_n and b_n).
  den_resid_df <- n_effective - p
  dispersion_cal <- S_marg / den_resid_df
  b_0_S_marg_formula <- 0.5 * S_marg * (n_prior +k+ p - 2L) / den_resid_df

  if (!is.finite(dispersion_cal) || dispersion_cal <= 0) {
    stop("compute_gaussian_prior: calibrated dispersion (S_marg/(n_effective-p)) is missing or not positive.", call. = FALSE)
  }
  if (!is.finite(b_0_S_marg_formula) || b_0_S_marg_formula <= 0) {
    stop(
      "compute_gaussian_prior: prior rate b_0 is missing or not positive. ",
      "Require n_prior + k + p > 2. Got n_prior = ", n_prior, ", k = ", k, ", p = ", p, ".",
      call. = FALSE
    )
  }
  rate <- b_0_S_marg_formula

  ## Prior rate for dGamma / fixed-beta path: same scaling as \code{rate}, but with
  ## RSS at Zellner blend \eqn{\beta_\star=(1-\mathrm{pwt})\hat\beta+\mathrm{pwt}\mu},
  ## \eqn{\mathrm{pwt}=n_{\mathrm{prior}}/(n_{\mathrm{prior}}+n_{\mathrm{effective}})}.
  pwt_scalar <- n_prior / (n_prior + n_effective)
  beta_star <- (1 - pwt_scalar) * bhat + pwt_scalar * mu_num
  res_star <- as.numeric(Y) - as.numeric(X %*% beta_star) - as.numeric(offset)
  rss_star <- sum(as.numeric(weights) * res_star^2)
  if (!is.finite(rss_star) || rss_star <= 0) {
    stop(
      "compute_gaussian_prior: weighted RSS at default coefficient blend must be strictly positive.",
      call. = FALSE
    )
  }
  rate_gamma <- 0.5 * rss_star * (n_prior + k + p - 2L) / den_resid_df
  if (!is.finite(rate_gamma) || rate_gamma <= 0) {
    stop("compute_gaussian_prior: computed rate_gamma must be strictly positive.", call. = FALSE)
  }
  shape_ING <- shape + p / 2

  Sigma_calibrated <- (n_effective / n_prior) * dispersion_cal * Ginv
  dimnames(Sigma_calibrated) <- list(colnames(X), colnames(X))

  dispersion <- dispersion_cal
  Sigma <- Sigma_calibrated
  Sigma_0_out <- Sigma_0
  dimnames(Sigma_0_out) <- list(colnames(X), colnames(X))

  if (!is.null(Sigma_input)) {
    if (!is.matrix(Sigma_input) || nrow(Sigma_input) != p || ncol(Sigma_input) != p || anyNA(Sigma_input)) {
      stop("compute_gaussian_prior: Sigma must be NULL or a numeric [p x p] matrix with no missing values.", call. = FALSE)
    }
    Sigma <- Sigma_input
    if (!is.null(dispersion_input)) {
      dispersion <- dispersion_input
    }
    dimnames(Sigma) <- list(colnames(X), colnames(X))
    Sigma_0_out <- Sigma / dispersion
    dimnames(Sigma_0_out) <- list(colnames(X), colnames(X))
  } else if (!is.null(dispersion_input)) {
    dispersion <- dispersion_input
    Sigma <- Sigma_0_out * dispersion
    dimnames(Sigma) <- list(colnames(X), colnames(X))
  }

  ## Step F: return calibrated outputs.
  list(
    dispersion = dispersion,
    shape = shape,
    shape_ING = shape_ING,
    rate = rate,
    rate_gamma = rate_gamma,
    Sigma = Sigma,
    Sigma_0 = Sigma_0_out
  )
}
