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
