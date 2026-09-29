## Re-exports of insight generics extended by glmbayes (see R/insight-methods.R)
## so they can be called unqualified after library(glmbayes), without also
## attaching insight. insight is a lightweight, low-dependency package and is
## listed in Imports (not Suggests) specifically to support this re-export.

#' @importFrom insight model_info
#' @export
insight::model_info

#' @importFrom insight get_parameters
#' @export
insight::get_parameters

#' @importFrom insight find_parameters
#' @export
insight::find_parameters

#' @importFrom insight find_algorithm
#' @export
insight::find_algorithm

#' @importFrom insight get_data
#' @export
insight::get_data

#' @importFrom insight get_priors
#' @export
insight::get_priors

## bayestestR generics extended by glmbayes (see R/bayestestR-methods.R).
## bayestestR's own Imports are just insight + datawizard + base packages
## (all lightweight), and insight is already a hard dependency above, so
## moving bayestestR from Suggests to Imports here adds only one real
## transitive dependency (datawizard).

#' @importFrom bayestestR simulate_prior
#' @export
bayestestR::simulate_prior

#' @importFrom bayestestR check_prior
#' @export
bayestestR::check_prior

#' @importFrom bayestestR describe_prior
#' @export
bayestestR::describe_prior

## Type (2): canonical implementation in glmbayesCore; help from @inherit (local roxygen in R/prior.R is commented out).
#' @inherit glmbayesCore::Prior_Setup params return details title description references examples format note
#' @family prior
#' @seealso
#' \code{\link{pfamily}} for prior-family objects and the constructors
#' \code{\link{dNormal}}, \code{\link{dNormal_Gamma}}, \code{\link{dGamma}},
#' and \code{\link{dIndependent_Normal_Gamma}}.
#'
#' \code{\link{glmb}}, \code{\link{lmb}} for formula-based fits with a
#' \code{pfamily} built from \code{Prior_Setup()} output; \code{\link{rglmb}},
#' \code{\link{rlmb}} for matrix-based sampling that consumes the same prior
#' structure; \code{\link[glmbayesCore]{simfuncs}} for functions that take a \code{prior_list}
#' assembled from those components (including
#' \code{\link{rindepNormalGamma_reg}} for
#' \code{\link{dIndependent_Normal_Gamma}()}).
#' \code{\link{multi_prior_setup}} for a matrix/cbind response with Gaussian;
#' use with \code{\link{lmb}}
#' \code{Prior_Setup} per column.
#'
#' \insertCite{zellner1986gprior}{glmbayes};
#' \insertCite{Raiffa1961}{glmbayes};
#' \insertCite{Gelman2013}{glmbayes};
#' \insertCite{McCullagh1989}{glmbayes};
#' \insertCite{glmbayesChapter03}{glmbayes};
#' \insertCite{glmbayesChapterA12}{glmbayes}.
#' @export
Prior_Setup <- glmbayesCore::Prior_Setup
