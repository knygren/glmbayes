#' Parallel Truncated‐Normal Regression (Standard Worker)
#'
#' Wraps the C++ function `rnnorm_reg_std_cpp_parallel` for fast, threaded sampling.
#'
#' @param n integer; number of observations.
#' @param y numeric vector of length n; response values.
#' @param x numeric matrix with n rows; predictors.
#' @param mu numeric matrix (p × …); prior means for each coefficient.
#' @param P numeric matrix or block‐matrix; prior precision/covariance.
#' @param alpha numeric vector of length p; regression coefficients.
#' @param wt numeric vector of length n; observation weights.
#' @param f2 R function; auxiliary update callback used internally.
#' @param Envelope list with the following elements:
#'   - PLSD: numeric vector of slice‐sampling probabilities.  
#'   - loglt: numeric matrix of lower‐tail log‐densities.  
#'   - logrt: numeric matrix of upper‐tail log‐densities.  
#'   - cbars: numeric matrix of truncation bounds.  
#'   - LLconst: numeric vector of log‐likelihood offsets.  
#' @param family character string; distribution family (e.g. "gaussian").
#' @param link character string; link function (e.g. "identity").
#' @param progbar integer (0 or 1); whether to display a progress bar.
#'
#' @return A list with components:
#'   - out: numeric matrix (n × p) of sampled latent values.  
#'   - draws: numeric vector of length n of final draws.  
#'
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib glmbayes, .registration = TRUE
#' @export
rnnorm_reg_std_cpp_parallel <- function(n,
                                        y,
                                        x,
                                        mu,
                                        P,
                                        alpha,
                                        wt,
                                        f2,
                                        Envelope,
                                        family,
                                        link,
                                        progbar = 1L) {
  # insure RcppParallel symbols are registered
  RcppParallelLibs()
  
  .Call(
    "_glmbayes_rnnorm_reg_std_cpp_parallel",
    as.integer(n),
    as.numeric(y),
    as.matrix(x),
    as.matrix(mu),
    as.matrix(P),
    as.numeric(alpha),
    as.numeric(wt),
    f2,
    Envelope,
    as.character(family),
    as.character(link),
    as.integer(progbar)
  )
}