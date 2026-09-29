#' Envelope Centering for Bayesian Gaussian Regression
#'
#' @description
#' `EnvelopeCentering()` computes an initial dispersion and the expected
#' posterior weighted RSS (closed form under the Normal posterior for
#' coefficients) for use in envelope construction when the dispersion is
#' unknown. The dispersion-anchoring loop updates dispersion from the Gamma
#' posterior using that expected RSS each iteration.
#' This step is typically called inside \code{rIndepNormalGammaReg()} before
#' \code{\link{EnvelopeOrchestrator}}, but may be used directly for diagnostics
#' or custom workflows.
#'
#' @param y Numeric response vector of length \code{m}.
#' @param x Numeric design matrix of dimension \code{m * p}.
#' @param mu Numeric vector of prior means (length \code{p}).
#' @param P Numeric matrix of prior precision (\code{p * p}).
#' @param offset Numeric vector of length \code{m}. Use \code{rep(0, m)} for none.
#' @param wt Numeric vector of prior weights.
#' @param shape Numeric. Shape parameter of the Gamma prior for the dispersion.
#' @param rate Numeric. Rate parameter of the Gamma prior for the dispersion.
#' @param Gridtype Integer. Grid construction method (default \code{2}).
#' @param verbose Logical. Reserved for API compatibility; currently unused in C++.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{dispersion}}{Numeric. Anchored dispersion value.}
#'   \item{\code{RSS_post}}{Numeric. Expected posterior weighted RSS (closed form;
#'     last iteration).}
#' }
#'
#' @details
#' The function first obtains an initial dispersion via \code{lm.wfit} residual
#' variance, then iteratively: (1) computes the expected weighted RSS under the
#' Normal posterior (closed form), (2) updates the
#' dispersion via the Gamma posterior using the expected RSS. The result is used
#' as \code{dispersion2} and \code{RSS_Post2} in downstream envelope construction
#' (e.g., \code{\link{EnvelopeOrchestrator}}).
#'
#' This anchors the joint Normal--Gamma accept--reject construction in
#' \insertCite{Nygren2006}{glmbayes}; see vignettes \code{Chapter-A07},
#' \code{Chapter-A11}, and \insertCite{glmbayesChapterA08,glmbayesIndNormGammaVignette}{glmbayes}.
#'
#' @seealso
#' \code{\link{EnvelopeOrchestrator}} for envelope construction;
#' \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeDispersionBuild}};
#' \code{\link{rindepNormalGamma_reg}} for the full simulation routine;
#' \code{\link{rlmb}} for the user-facing linear-model interface.
#'
#' @references
#' \insertAllCited{}
#'
#' @example inst/examples/Ex_EnvelopeCentering.R
#'
#' @export
EnvelopeCentering <- function(y, x, mu, P, offset, wt, shape, rate,
                             Gridtype = 2L, verbose = FALSE) {
  .EnvelopeCentering_cpp(y, x, mu, P, offset, wt, shape, rate,
                         Gridtype, verbose)
}


#' Envelope Construction Orchestrator for Bayesian Gaussian Regression
#'
#' @description
#' `EnvelopeOrchestrator()` provides a unified interface for constructing the
#' fixed‑dispersion and dispersion‑aware envelopes used in likelihood‑subgradient
#' simulation for Bayesian Gaussian regression with Normal–Gamma priors.
#'
#' This function coordinates:
#'
#' * fixed‑dispersion envelope construction via \link[glmbayes]{EnvelopeBuild},
#' * dispersion‑refined envelope construction via \link[glmbayes]{EnvelopeDispersionBuild},
#' * envelope sorting and reindexing via \link[glmbayes]{EnvelopeSort}, and
#' * UB‑list alignment (reordered `lg_prob_factor` and `UB2min`).
#'
#' It is typically used inside *.cpp routines such as
#' \code{rIndepNormalGammaReg()}, but may also be called directly for
#' diagnostics, envelope visualization, or custom simulation workflows.
#'
#' @param bstar2 Numeric vector. Posterior mode of the standardized regression
#'   coefficients (from the standardized model).
#' @param A Numeric matrix. Posterior precision matrix (Hessian) at the mode.
#' @param y Numeric response vector of length \code{m}.
#' @param x2 Numeric matrix of standardized predictors (\code{m × p}).
#' @param mu2 Numeric vector. Standardized prior mean (typically a zero vector).
#' @param P2 Numeric matrix. Standardized prior precision component moved into
#'   the log‑likelihood.
#' @param alpha Numeric vector. Offset‑adjusted mean component.
#' @param wt Numeric vector of prior weights.
#' @param n Integer. Number of envelope grid points or simulation draws.
#' @param Gridtype Integer specifying the envelope grid construction method for
#'   API compatibility. The C++ orchestrator \strong{overrides} this to \code{3L}
#'   (full \eqn{3^{p}} grid): unknown dispersion does not use smaller grids.
#' @param n_envopt Optional integer. Effective sample size passed to
#'   `EnvelopeOpt` during grid construction. Larger values encourage tighter
#'   envelopes.
#' @param shape Numeric. Shape parameter of the Gamma prior for the dispersion.
#' @param rate Numeric. Rate parameter of the Gamma prior for the dispersion.
#' @param RSS_Post2 Numeric. Expected posterior weighted RSS used to anchor the
#'   dispersion axis (typically \code{centering_out$RSS_post} from
#'   \code{\link{EnvelopeCentering}} inside \code{\link{rindepNormalGamma_reg}};
#'   see vignette \code{Chapter-A11}).
#' @param RSS_ML Numeric. Maximum‑likelihood residual sum of squares.
#' @param max_disp_perc Numeric in \code{(0,1)}. Tail probability used to
#'   determine dispersion bounds when not explicitly supplied.
#' @param disp_lower Optional numeric. Lower bound for the dispersion
#'   (\eqn{\sigma^2}). If supplied, overrides quantile‑based bounds.
#' @param disp_upper Optional numeric. Upper bound for the dispersion
#'   (\eqn{\sigma^2}). Must be strictly greater than \code{disp_lower}.
#' @param use_parallel Logical. Whether to allow parallel computation inside
#'   \link[glmbayes]{EnvelopeDispersionBuild}.
#' @param use_opencl Logical. Whether to allow OpenCL acceleration inside
#'   \link[glmbayes]{EnvelopeBuild}.
#' @param verbose Logical. Whether to print detailed progress and timing
#'   messages.
#'
#' @return
#' A list with components:
#'
#' \describe{
#'   \item{\code{Env}}{The fully constructed and sorted envelope, including the
#'     PLSD component inserted by the dispersion‑aware refinement step.}
#'   \item{\code{gamma_list}}{Updated Gamma‑prior parameters for the dispersion
#'     (shape, rate, and dispersion bounds).}
#'   \item{\code{UB_list}}{Updated UB‑list including reordered
#'     \code{lg_prob_factor} and \code{UB2min}.}
#'   \item{\code{diagnostics}}{Diagnostic quantities returned by
#'     \link[glmbayes]{EnvelopeDispersionBuild}, useful for debugging or envelope
#'     visualization.}
#'   \item{\code{low}}{Lower dispersion bound used.}
#'   \item{\code{upp}}{Upper dispersion bound used.}
#' }
#'
#' @details
#' `EnvelopeOrchestrator()` is the **envelope-construction stage** for Bayesian
#' **Gaussian regression with an independent Normal--Gamma prior** on
#' \eqn{(\beta, \phi)} (dispersion \eqn{\phi}; precision \eqn{\tau = 1/\phi} in
#' much of the theory). It is implemented in \file{src/EnvelopeOrchestrator.cpp}
#' and composes direct C++ calls to `EnvelopeBuild` and
#' `EnvelopeDispersionBuild` with an R call to \code{\link{EnvelopeSort}}.
#'
#' **What this function does not do.** It does **not** run the iterative
#' dispersion centering loop (\code{\link{EnvelopeCentering}}), **not** optimize
#' the posterior mode or Hessian, **not** standardize the model
#' (\code{\link{glmb_Standardize_Model}}), and **not** draw posterior samples.
#' Those steps are performed by \code{\link{rindepNormalGamma_reg}} (see vignette
#' \code{Chapter-A11}) before and after the orchestrator. Inputs such as
#' \code{bstar2}, \code{A}, \code{x2}, \code{mu2}, and \code{P2} must therefore
#' already be in **standard form** for the coefficient subproblem, exactly as
#' passed from that workflow.
#'
#' **What the return value is for.** The returned \code{Env}, \code{gamma_list},
#' and \code{UB_list} are consumed by the internal standardized samplers
#' `rIndepNormalGammaReg_std` and `rIndepNormalGammaReg_std_parallel` in
#' \file{src/rIndepNormalGammaReg.cpp}, which implement the joint
#' accept--reject procedure over \eqn{(\beta, \phi)}. Theory for the dispersion
#' envelope and bounding arguments is in vignette \code{Chapter-A07}; the
#' end-to-end implementation map is in \code{Chapter-A11}. The coefficient-only
#' likelihood-subgradient envelope (\insertCite{Nygren2006}{glmbayes}) is
#' documented under \code{\link{EnvelopeBuild}} and vignette \code{Chapter-A08}.
#'
#' The function does **not** perform simulation. Simulation is carried out
#' afterward via \code{.rIndepNormalGammaReg_std_cpp()} or
#' \code{.rIndepNormalGammaReg_std_parallel_cpp()}, depending on
#' \code{use_parallel}.
#'
#' @section Use of the envelope during sampling:
#'
#' After \code{EnvelopeOrchestrator()} returns, \code{\link{rindepNormalGamma_reg}}
#' delegates iid simulation to `rIndepNormalGammaReg_std` (serial) or
#' `rIndepNormalGammaReg_std_parallel` (parallel). These routines are **not
#' exported**; they are the direct analogues of the fixed-dispersion path
#' \code{.rNormalGLM_std_cpp()} for GLMs, but for the **joint** posterior
#' \eqn{\pi(\beta, \phi \mid y)} under the independent Normal--Gamma prior.
#'
#' **Dominating proposal (conceptual).** The envelope list \code{Env} still
#' describes a **mixture of restricted multivariate Normal** proposal pieces for
#' the **standardized** regression coefficients, with mixture weights
#' \eqn{\tilde{p}_j} stored in \code{PLSD}. After
#' \code{\link{EnvelopeDispersionBuild}}, those weights and the per-face
#' constants are adjusted so that, together with a **truncated inverse-Gamma**
#' (dispersion) proposal derived from \code{gamma_list}, the joint proposal
#' dominates the target posterior on the truncated dispersion interval
#' \code{[low, upp]}. \code{vignette("Chapter-A07", package = "glmbayes")} derives
#' the dispersion-related bounds; \code{vignette("Chapter-A11", package = "glmbayes")}
#' records how \code{UB_list} entries enter the code.
#'
#' **One accept--reject iteration** (standardized coordinates) proceeds as follows:
#'
#' \enumerate{
#'   \item **Draw a mixture component (face)** \eqn{J}. An index \eqn{J} is
#'         drawn from the discrete distribution with probabilities \code{PLSD}.
#'   \item **Propose coefficients** \eqn{\beta^\star}. Conditional on \eqn{J},
#'         each coordinate is drawn from the **restricted Normal** used in the
#'         fixed-dispersion construction: cumulative-normal tail probabilities
#'         \code{loglt[J, ]}, \code{logrt[J, ]}, and subgradient shift
#'         \code{-cbars[J, ]} (internal \code{rnorm_ct} truncated Normal
#'         sampling, same structural role as \code{ctrnorm_cpp()} in the GLM path).
#'   \item **Propose dispersion** \eqn{\phi}. A draw is taken from the truncated
#'         inverse-Gamma / Gamma piece defined by \code{shape3}, \code{rate2},
#'         \code{disp_lower}, and \code{disp_upper} in \code{gamma_list}
#'         (\code{rinvgamma_ct_safe}).
#'   \item **Re-weight the likelihood for** \eqn{\phi}. Observation weights in
#'         the Gaussian log-likelihood are scaled by \eqn{1/\phi}
#'         (\code{wt2 = wt / dispersion} in the C++ sources).
#'   \item **Dispersion-adjusted tangency.** Because the tangency point for the
#'         linear upper bound depends on dispersion, the code recomputes a
#'         **face-specific** \eqn{\bar{\theta}_J(\phi)} via
#'         \code{Inv_f3_with_disp} (using a one-time cache from
#'         \code{Inv_f3_precompute_disp} built from \code{cbars} and the data).
#'         The negative log-likelihood at that point feeds the **UB1** tangent term.
#'   \item **Log-likelihood at the proposal.** Compute
#'         \eqn{-\log f(y \mid \beta^\star, \phi)} with the same \eqn{\phi} and
#'         scaled weights (output \code{LL_Test} in the serial implementation).
#' }
#'
#' **Acceptance inequality (structure).** Write \eqn{\ell(\beta,\phi)} for the
#' Gaussian log-likelihood (weighted, with offset). The serial sampler forms
#' \eqn{\mathrm{UB1}} from the tangent to \eqn{-\ell} at
#' \eqn{\bar{\theta}_J(\phi)} along subgradient \eqn{c_J = }\code{cbars[J, ]}:
#' \deqn{
#'   \mathrm{UB1} =
#'   -\ell\!\big(\bar{\theta}_J(\phi), \phi\big)
#'   - c_J^\top \big(\beta^\star - \bar{\theta}_J(\phi)\big).
#' }
#' Additional terms bound RSS variation along \eqn{\phi} (\strong{UB2}, using
#' \code{RSS_Min} and \code{UB2min} from \code{UB_list}) and dispersion-axis
#' majorization (\strong{UB3A}, \strong{UB3B}) built from \code{lg_prob_factor},
#' \code{lmc1}, \code{lmc2}, \code{lm_log1}, \code{lm_log2},
#' \code{max_New_LL_UB}, and \code{max_LL_log_disp}. With
#' \eqn{L_{\mathrm{test}} = -\ell(\beta^\star,\phi)}, define
#' \eqn{T_1 = L_{\mathrm{test}} - \mathrm{UB1}} and
#' \eqn{T = T_1 - (\mathrm{UB2} + \mathrm{UB3A} + \mathrm{UB3B})}. The code draws
#' \eqn{U_2 \sim \mathrm{Unif}(0,1)} and **accepts** \eqn{(\beta^\star, \phi)} when
#' \deqn{
#'   T - \log(U_2) \ge 0.
#' }
#' Serial and parallel workers use the same logical decomposition up to
#' implementation detail. Under the construction in
#' \code{\link{EnvelopeDispersionBuild}}, the terms are arranged so that
#' \eqn{T_1 \le 0} and \eqn{\mathrm{UB2}}, \eqn{\mathrm{UB3A}}, \eqn{\mathrm{UB3B}}
#' are nonnegative up to controlled numerical slack. Iteration counts are stored
#' in \code{iters_out}.
#'
#' **Mapping orchestrator outputs to the sampler.**
#' \itemize{
#'   \item \code{Env$PLSD}: mixture probabilities over envelope faces for Step 1.
#'   \item \code{Env$loglt}, \code{Env$logrt}, \code{Env$cbars}: restricted
#'         Normal proposal for \eqn{\beta^\star} in Step 2.
#'   \item \code{Env$GridIndex}, \code{Env$thetabars}, \code{Env$logU},
#'         \code{Env$logP}: same role as in \code{\link{EnvelopeBuild}} for the
#'         coefficient mixture; dispersion refinement may update \code{PLSD}
#'         before sorting.
#'   \item \code{gamma_list}: truncated dispersion proposal parameters
#'         (\code{shape3}, \code{rate2}, bounds) for Step 3.
#'   \item \code{UB_list}: global and per-face constants (\code{RSS_Min},
#'         \code{UB2min}, \code{lg_prob_factor}, linear \code{lmc}/\code{lm_log}
#'         pieces) for \eqn{\mathrm{UB2}}, \eqn{\mathrm{UB3A}}, \eqn{\mathrm{UB3B}}.
#'   \item \code{low}, \code{upp}: dispersion interval endpoints (duplicated
#'         from \code{gamma_list} for convenience).
#' }
#' Unlike the fixed-dispersion GLM sampler, this path does **not** apply the
#' stored \code{LLconst} vector directly in the acceptance test; the tangent
#' piece is recomputed as \eqn{\mathrm{UB1}} once \eqn{\phi} and
#' \eqn{\bar{\theta}_J(\phi)} are known.
#'
#' @section Algorithmic steps:
#'
#' The orchestrator implements the **independent Normal--Gamma** envelope
#' pipeline: first a **coefficient** envelope at a **dispersion anchor**
#' (\insertCite{Nygren2006}{glmbayes}; vignette \code{Chapter-A08}), then
#' **dispersion-aware** refinement (\code{Chapter-A07}), then sorting. Steps 3--8
#' repeat the internal logic of \code{\link{EnvelopeBuild}} (same formulas on that
#' help page); here the likelihood is **Gaussian** with **identity** link, weights
#' are \eqn{w_i / d_\star} with \eqn{d_\star} from the anchor below, and the first
#' pass uses \code{sortgrid = FALSE} so sorting runs **after** dispersion
#' refinement.
#'
#' 1. **Force full grid for unknown dispersion.** The argument \code{Gridtype}
#'    is overridden to \code{3L} so the coefficient grid always uses the full
#'    \eqn{3^{p}} partition (implementation policy in
#'    \file{src/EnvelopeOrchestrator.cpp}).
#'
#' 2. **Anchor dispersion and rescale weights for** \code{EnvelopeBuild}.
#'    Let \eqn{n_w = \sum_i w_i}. With prior hyperparameters \code{shape}
#'    (\eqn{a_0}) and \code{rate} (\eqn{b_0}) and centered RSS \code{RSS_Post2},
#'    define \eqn{s = a_0 + n_w/2} and \eqn{r = b_0 + \mathrm{RSS}_{\mathrm{post}}/2}
#'    where \eqn{\mathrm{RSS}_{\mathrm{post}}} denotes \code{RSS_Post2} (the C++
#'    code names the scalars \code{shape2} and \code{rate3}). The dispersion anchor
#'    is \eqn{d_\star = r/(s - 1)}, and observation weights \eqn{w_i} passed into
#'    the embedded \code{EnvelopeBuild} call are scaled by \eqn{1/d_\star}. This ties
#'    the coefficient envelope to the Gamma posterior for the precision conditional
#'    on the centered RSS (Chapters A07, A11).
#'
#' 3. **Compute width parameters** \eqn{\omega_i} **from the diagonal precision
#'    matrix.** Let \eqn{\theta^{\ast}} be the standardized posterior mode. For
#'    each dimension \eqn{i},
#'    \deqn{
#'      \omega_{i} :=
#'      \frac{\sqrt{2} - \exp\!\big(-1.20491 - 0.7321\,\sqrt{0.5 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2}}\big)}
#'           {\sqrt{1 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2}}}.
#'    }
#'    Here \eqn{f} is the weighted Gaussian log-posterior for \eqn{\beta} at the
#'    anchored dispersion.
#'
#' 4. **Construct intervals and the** \eqn{3^{p}} **partition** around
#'    \eqn{\theta^\star}. Set
#'    \deqn{
#'      \ell_{i,1} = \theta^{\ast}_{i} - 0.5\,\omega_{i}, \quad
#'      \ell_{i,2} = \theta^{\ast}_{i} + 0.5\,\omega_{i},
#'    }
#'    and
#'    \deqn{
#'      A_{i,1} = (-\infty,\ell_{i,1}), \quad
#'      A_{i,2} = [\ell_{i,1},\ell_{i,2}], \quad
#'      A_{i,3} = (\ell_{i,2},\infty).
#'    }
#'    With \eqn{J = \prod_{i=1}^{p} \{1,2,3\}} and \eqn{j = (j_1,\ldots,j_p)},
#'    \eqn{A^{\ast}_{j} = \prod_{i=1}^{p} A_{i,j_i}} partitions standardized
#'    coefficient space.
#'
#' 5. **Select tangency points** \eqn{\bar{\theta}_j} **per cell** (left / mode /
#'    right of each interval). For index sets \eqn{C_{j1},C_{j2},C_{j3}} by
#'    coordinate,
#'    \deqn{
#'      \bar{\theta}_{j,i} =
#'      \begin{cases}
#'        \theta^{\ast}_{i} - \omega_{i}, & i \in C_{j1}, \\
#'        \theta^{\ast}_{i},              & i \in C_{j2}, \\
#'        \theta^{\ast}_{i} + \omega_{i}, & i \in C_{j3}.
#'      \end{cases}
#'    }
#'
#' 6. **Evaluate negative log-likelihood and gradient at each grid point.**
#'    Subgradients \eqn{c(\bar{\theta}_j)} and negative log-likelihoods define
#'    the likelihood-subgradient envelope pieces (\insertCite{Nygren2006}{glmbayes};
#'    \code{Chapter-A08}). CPU: \code{f2_f3_non_opencl}; GPU (optional):
#'    \code{f2_f3_opencl}.
#'
#' 7. **Call** \code{EnvelopeSet_Grid_C2_pointwise} **and**
#'    \code{EnvelopeSet_LogP_C2} **(C++ pipeline)** to obtain restricted Normal
#'    log-densities, mixture log-probabilities, and constants as in Remarks 5--6
#'    of the JASA paper (same as \code{\link{EnvelopeBuild}}).
#'
#' 8. **Normalize to** \code{PLSD} **without sorting.** The embedded
#'    \code{EnvelopeBuild} call sets \code{sortgrid = FALSE} so an intermediate
#'    sort is not wasted before \code{\link{EnvelopeDispersionBuild}} revises
#'    mixture weights for the joint \eqn{(\beta,\phi)} target and
#'    \code{\link{EnvelopeSort}} runs once at the end.
#'
#' 9. **Call** \code{EnvelopeDispersionBuild} **(C++).** Pass the coefficient
#'    envelope list, prior \code{shape}, \code{rate}, standardized \code{P2},
#'    data \code{y}, \code{x2}, \code{alpha}, \code{mu2}, \code{wt},
#'    \code{RSS_Post2}, \code{RSS_ML}, dispersion controls
#'    (\code{max_disp_perc}, optional bounds), and \code{use_parallel}. This
#'    constructs the dispersion truncation interval, updates Gamma proposal
#'    parameters (\code{gamma_list}), computes \code{UB_list}, and returns
#'    \code{Env_out} with adjusted \code{PLSD}. See \code{Chapter-A07} and
#'    \code{Chapter-A11}, Section 3.3.
#'
#' 10. **Call** \code{\link{EnvelopeSort}} **(R).** Reorder envelope components and
#'     align \code{lg_prob_factor} and \code{UB2min} with the sorted indexing. If
#'     sorting cannot allocate safely, the implementation falls back to unsorted
#'     \code{Env_out} with UB fields patched (\file{src/EnvelopeOrchestrator.cpp}).
#'
#' 11. **Return** \code{Env}, \code{gamma_list}, \code{UB_list},
#'     \code{diagnostics}, \code{low}, and \code{upp} for the standardized
#'     samplers.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' * \link[glmbayes]{EnvelopeBuild} – fixed‑dispersion envelope construction
#' * \link[glmbayes]{EnvelopeDispersionBuild} – dispersion‑aware envelope refinement
#' * \link[glmbayes]{EnvelopeSort} – envelope sorting and reindexing
#' * \link[glmbayes]{EnvelopeCentering} – \code{RSS_Post2} and dispersion anchor
#' * \link[glmbayes]{glmb_Standardize_Model} – standardized inputs for the orchestrator
#' * \code{\link{rindepNormalGamma_reg}} – full Normal–Gamma workflow (R + C++)
#' * \code{\link{rlmb}}, \code{\link{rglmb}}, \code{\link{simfuncs}} – higher-level sampling entry points
#' * Vignettes \code{Chapter-A07}, \code{Chapter-A08}, \code{Chapter-A11}; cited as
#'   \insertCite{Nygren2006,glmbayesChapterA08,glmbayesIndNormGammaVignette}{glmbayes}
#'
#' @example inst/examples/Ex_EnvelopeOrchestrator.R
#'
#' @export

EnvelopeOrchestrator <- function(bstar2,
                                 A,
                                 y,
                                 x2,
                                 mu2,
                                 P2,
                                 alpha,
                                 wt,
                                 n,
                                 Gridtype,
                                 n_envopt,
                                 shape,
                                 rate,
                                 RSS_Post2,
                                 RSS_ML,
                                 max_disp_perc,
                                 disp_lower,
                                 disp_upper,
                                 use_parallel = TRUE,
                                 use_opencl  = FALSE,
                                 verbose     = FALSE) {
  
  # --- NEW: Call the C++ orchestrator directly ---
  out_cpp <- .EnvelopeOrchestrator_cpp(
    bstar2      = bstar2,
    A           = A,
    y           = y,
    x2          = x2,
    mu2         = as.matrix(mu2, ncol = 1),
    P2          = P2,
    alpha       = alpha,
    wt          = wt,
    n           = n,
    Gridtype    = Gridtype,
    n_envopt    = n_envopt,
    shape       = shape,
    rate        = rate,
    RSS_Post2   = RSS_Post2,
    RSS_ML      = RSS_ML,
    max_disp_perc = max_disp_perc,
    disp_lower  = disp_lower,
    disp_upper  = disp_upper,
    use_parallel = use_parallel,
    use_opencl   = use_opencl,
    verbose      = verbose
  )
  
  if (verbose) {
    cat("[EnvelopeOrchestrator] Using C++ orchestrator output\n")
  }
  
  return(out_cpp)
}

