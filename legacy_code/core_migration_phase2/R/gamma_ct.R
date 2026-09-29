#' The Central Gamma Distribution
#'
#' Distribution function and random generation for the center (between a lower
#' and an upper bound) of the Gamma distribution with shape and rate parameters.
#' These functions provide numerically stable evaluation and sampling when the
#' truncation interval is narrow or when the Gamma density is highly skewed.
#'
#' @name Gamma_ct
#' @aliases
#' Gamma_ct
#' pgamma_ct
#' rgamma_ct
#' ctrgamma
#'
#' @param n Number of draws to generate. If \code{length(n) > 1}, the length
#'   is taken to be the number required.
#' @param shape Shape parameter of the Gamma distribution.
#' @param rate Rate parameter of the Gamma distribution.
#' @param lower_prec Lower truncation point on the precision scale. If
#'   \code{NULL}, no lower truncation is applied.
#' @param upper_prec Upper truncation point on the precision scale. If
#'   \code{NULL}, no upper truncation is applied.
#'   
#' @return
#' For \code{pgamma_ct}, a vector of probabilities corresponding to the mass of
#' the Gamma distribution between \code{a} and \code{b}.  
#' For \code{rgamma_ct} or \code{ctrgamma}, a vector of length \code{nn}
#' containing random draws from the Gamma distribution restricted to the
#' interval \code{[a, b]} (or \code{[lower_prec, upper_prec]} on the precision
#' scale).
#'
#' @details
#' The function \code{pgamma_ct} computes the probability mass between a lower
#' bound \code{a} and an upper bound \code{b} under a Gamma density with the
#' specified shape and rate parameters. This is particularly useful when the
#' interval \code{b - a} is small, where the naive computation
#' \code{pgamma(b) - pgamma(a)} may underflow to zero even when the true
#' probability is positive.
#'
#' The function \code{ctrgamma} provides a numerically robust sampler for the
#' Gamma distribution under one-sided or two-sided truncation. It handles:
#'
#' \itemize{
#'   \item no truncation (reducing to \code{rgamma})
#'   \item lower truncation only
#'   \item upper truncation only
#'   \item two-sided truncation with \code{lower_prec < upper_prec}
#'   \item exact degeneracy when the truncation interval collapses
#'   \item numerical degeneracy when the Gamma CDF collapses in floating point
#' }
#'
#' All computations are performed on the log scale using stable log–CDF and
#' log–sum–exp transformations. This avoids the catastrophic cancellation that
#' occurs when the Gamma CDF values at the truncation points are extremely close.
#'
#' These functions are primarily intended for use in hierarchical Bayesian
#' models where precision parameters are updated under tight truncation
#' constraints, and where numerical stability is essential for reliable sampling
#' performance. They are used in envelope-based dispersion sampling
#' \insertCite{Nygren2006}{glmbayes}.
#'
#' @seealso \code{\link{Normal_ct}}, \code{\link{InvGamma_ct}}, \code{\link{EnvelopeDispersionBuild}}
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_Gamma_ct.R
#' @rdname Gamma_ct
#' @order 1
#' @export


# Constrained truncated gamma sampler on precision scale (numerically robust)
rgamma_ct <- function(n, shape, rate, lower_prec = NULL, upper_prec = NULL) {
  # ---- Case 0: no truncation ----
  if (is.null(lower_prec) && is.null(upper_prec)) {
    return(stats::rgamma(n, shape = shape, rate = rate))
  }
  
  # ---- Extract numeric bounds (support-aware) ----
  # Gamma support is [0, infinity), so clamp at 0 if lower is NULL
  L <- if (!is.null(lower_prec)) lower_prec else 0
  U <- if (!is.null(upper_prec)) upper_prec else Inf
  
  # ---- Exact numeric degeneracy in the bounds ----
  if (is.finite(L) && is.finite(U) && L == U) {
    # Interval truly collapsed
    return(rep(L, n))
  }
  
  # ---- One-sided truncation: lower only ----
  if (!is.null(lower_prec) && is.null(upper_prec)) {
    # Sample from Gamma conditioned on X >= L.
    # Work on the upper-tail CDF: T(x) = P(X >= x) = 1 - F(x)
    log_tail_L <- stats::pgamma(L, shape = shape, rate = rate,
                                lower.tail = FALSE, log.p = TRUE)
    # U ~ Unif(0,1), target tail probability: T_star = U * T(L)
    u      <- stats::runif(n)
    log_u  <- log(u)
    log_T  <- log_u + log_tail_L  # log T_star
    
    # Invert using qgamma on the upper tail with log.p
    out <- stats::qgamma(p = log_T, shape = shape, rate = rate,
                         lower.tail = FALSE, log.p = TRUE)
    return(out)
  }
  
  # ---- One-sided truncation: upper only ----
  if (is.null(lower_prec) && !is.null(upper_prec)) {
    # Sample from Gamma conditioned on X <= U.
    # Work on the lower-tail CDF: F(x) = P(X <= x)
    log_F_U <- stats::pgamma(U, shape = shape, rate = rate,
                             lower.tail = TRUE, log.p = TRUE)
    
    u      <- stats::runif(n)
    log_u  <- log(u)
    log_F  <- log_u + log_F_U  # log F_star
    
    out <- stats::qgamma(p = log_F, shape = shape, rate = rate,
                         lower.tail = TRUE, log.p = TRUE)
    return(out)
  }
  
  # ---- Two-sided truncation: L and U both finite ----
  # Now both lower_prec and upper_prec are non-NULL, and L < U (checked above).
  
  # Work on the lower-tail CDF: F(x) = P(X <= x)
  log_F_L <- stats::pgamma(L, shape = shape, rate = rate,
                           lower.tail = TRUE, log.p = TRUE)
  log_F_U <- stats::pgamma(U, shape = shape, rate = rate,
                           lower.tail = TRUE, log.p = TRUE)
  
  # Enforce ordering (should normally hold, but guard anyway)
  if (log_F_U < log_F_L) {
    tmp      <- log_F_L
    log_F_L  <- log_F_U
    log_F_U  <- tmp
    tmp      <- L
    L        <- U
    U        <- tmp
  }
  
  # If pgamma collapses the CDF in floating point, the conditional law
  # is essentially degenerate at L. We do not use an arbitrary epsilon;
  # we only test exact equality on the log scale.
  if (log_F_U == log_F_L) {
    return(rep(L, n))
  }
  
  # Mass of the interval (L, U] on the CDF scale, in log:
  # delta = F(U) - F(L), log_delta = logdiffexp(log_F_L, log_F_U)
  log_delta <- logdiffexp(log_F_L, log_F_U)
  
  # U ~ Unif(0,1), target unconditional CDF value:
  # F_star = F(L) + U * (F(U) - F(L))
  # Work in log-space via log-sum-exp of two terms:
  #   term1 = log_F_L
  #   term2 = log(U) + log_delta
  u      <- stats::runif(n)
  log_u  <- log(u)
  t1     <- log_F_L
  t2     <- log_u + log_delta
  
  m      <- pmax(t1, t2)
  log_F  <- m + log(exp(t1 - m) + exp(t2 - m))  # log(F_star)
  
  # Invert using qgamma on the lower tail with log.p
  out <- stats::qgamma(p = log_F, shape = shape, rate = rate,
                       lower.tail = TRUE, log.p = TRUE)
  out
}


