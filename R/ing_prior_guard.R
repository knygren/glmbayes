#' Implied effective prior sample size for an ING Gamma component
#'
#' Under the \code{Prior_Setup()} / Block~2 calibration
#' \code{shape = (n_prior + 1 + p) / 2}, invert to
#' \code{n_prior = 2 * shape - 1 - p}.
#'
#' @param shape Gamma shape (scalar numeric).
#' @param p Number of hyper-parameters (\code{length(mu)} / \code{ncol(x)}).
#' @return Scalar implied \code{n_prior}.
#' @noRd
.ing_n_prior_from_shape <- function(shape, p) {
  2 * shape - 1 - p
}

#' Stop when an ING dispersion prior carries more information than the data
#'
#' The dispersion envelope caps its log-tilt at \code{n_w/2} (Remark 4.1.3 of
#' the ING vignette), which requires \code{n_prior <= n_w} (equivalently
#' \code{pwt <= 0.5} under \code{Prior_Setup()} calibration).
#'
#' @param shape Gamma shape (scalar).
#' @param p Number of hyper-parameters.
#' @param n_w Effective observation count (\code{sum(weights)} or Block~2
#'   \code{J = length(group_levels)}).
#' @param detail Clause after \code{but} describing the data limit (no
#'   trailing period).
#' @param limit_label Symbol used in the requirement clause (e.g.
#'   \code{"n_w"} or \code{"J"}).
#' @param prefix Optional message prefix (component / function name).
#' @noRd
.ing_stop_if_prior_exceeds_data <- function(
    shape,
    p,
    n_w,
    detail,
    limit_label = "n_w",
    prefix = NULL
) {
  n_prior_implied <- .ing_n_prior_from_shape(shape, p)
  if (n_prior_implied > n_w) {
    msg <- paste0(
      prefix,
      "dIndependent_Normal_Gamma prior implies n_prior = ",
      signif(n_prior_implied, 4),
      " effective prior observations, but ",
      detail,
      ". The dispersion envelope requires n_prior <= ",
      limit_label,
      " (prior weight pwt <= 0.5); weaken the prior (smaller shape) ",
      "or supply more data."
    )
    stop(msg, call. = FALSE)
  }
  invisible(n_prior_implied)
}
