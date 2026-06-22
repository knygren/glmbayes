#' @aliases glmbayes
#'
#' @title glmbayes: Bayesian Generalized Linear Models with iid Sampling
#'
#' @description
#' `glmbayes` provides independent and identically distributed (iid) samples for
#' Bayesian generalized linear models (GLMs), serving as a Bayesian analogue to
#' the base `glm()` function. Supported likelihood families include Gaussian,
#' Poisson, Binomial, and Gamma models with log-concave likelihoods.
#'
#' @details
#' The main user-facing interface is `glmb()`, which mirrors the structure of
#' `glm()` and supports prior specification through `pfamily` objects. Lower-level
#' functions such as `rglmb()` and `rGamma_reg()` provide direct access to the
#' underlying samplers and can be used in block Gibbs sampling or hierarchical
#' model implementations.
#'
#' For an introduction to the package, examples, and a complete set of vignettes,
#' see:
#'
#' - README: <https://github.com/knygren/glmbayes#readme>
#' - All vignettes: `browseVignettes("glmbayes")`
#'
#' The package includes extensive documentation on model fitting, prior
#' construction, diagnostics, and optional GPU acceleration using OpenCL.
#'
#' **Releases:** This source tree is **0.9.6.9000** (in development). The current
#' **CRAN** release is **0.9.6** (`install.packages("glmbayes")`).
#' Source is available from GitHub; R-Universe (\url{https://knygren.r-universe.dev/glmbayes})
#' also builds binaries from that source.
#' Prebuilt CRAN and R-Universe binaries do not include OpenCL; GPU support
#' requires a source install once the host OpenCL environment is ready
#' (see \code{vignette("Chapter-16", "glmbayes")} for the three-step process).
#'
#' IID posterior simulation for non-Gaussian GLMs and several non-conjugate
#' linear-model setups uses the likelihood-subgradient envelope method of
#' \insertCite{Nygren2006}{glmbayes}. Introductory material and worked
#' examples are in \insertCite{glmbayesChapter00,glmbayesChapterA01}{glmbayes};
#' estimation and simulation background in
#' \insertCite{glmbayesChapterA02,glmbayesSimmethods,glmbayesChapterA08}{glmbayes};
#' prior derivations for \code{Prior_Setup()} in
#' \insertCite{glmbayesChapterA12}{glmbayes};
#' GPU/OpenCL topics in
#' \insertCite{glmbayesChapter12,glmbayesChapterA10}{glmbayes}.
#'
#' @section OpenCL startup checks:
#' In interactive sessions, attaching the package with \code{library(glmbayes)}
#' may emit a short \code{\link{packageStartupMessage}}
#' when \code{has_opencl()} is \code{FALSE} (typical for CRAN binaries) but a
#' GPU or OpenCL stack appears available on the host. OpenCL modelling paths
#' require a source install of \pkg{glmbayes} with OpenCL at compile time;
#' \code{has_opencl()} then reports whether that build succeeded. The note
#' confirms full CPU use and points to \code{vignette("Chapter-16")}. Machines
#' without a detectable GPU stack stay silent.
#' Set \code{options(glmbayes.quiet_opencl_startup = TRUE)} to suppress attach
#' notes (recommended for CI and \command{R CMD check}).
#'
#' @example inst/examples/Ex_glmbayes-package.R
#'
#' @seealso
#' Main interfaces: \code{\link{glmb}}, \code{\link{lmb}},
#' \code{\link{rglmb}}, \code{\link{rlmb}}; low-level simulation API
#' \code{\link{simfuncs}}; envelope construction \code{\link{EnvelopeBuild}}.
#'
#' Useful links:
#' \itemize{
#'   \item CRAN: <https://CRAN.R-project.org/package=glmbayes>
#'   \item GitHub: <https://github.com/knygren/glmbayes>
#'   \item R-Universe: <https://knygren.r-universe.dev/glmbayes>
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @author
#' Kjell Nygren
#'
#' @import stats Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom MASS mvrnorm
#' @importFrom Rdpack reprompt
#' @importFrom RcppParallel RcppParallelLibs
#' @import nmathopencl
#' @import opencltools
#' @useDynLib glmbayes, .registration = TRUE
"_PACKAGE"