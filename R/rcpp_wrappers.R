# -------------------------------------------------------------------------
#  Rcpp Interface Wrappers for glmbayes
#
#  These functions provide the minimal, strictly positional R → C++ bridges
#  required by the package.  Each wrapper mirrors the exact argument order
#  expected by the corresponding C++ routine and performs no preprocessing,
#  validation, or postprocessing.  Their sole purpose is to ensure that
#  high‑level R code calls the correct compiled symbol with the correct
#  signature.
#
#  All wrappers are internal:
#    - They are not part of the public API.
#    - They exist only to guarantee stable, explicit R–C++ boundaries.
#    - They prevent accidental reliance on .Call() with named arguments,
#      which R ignores, and which can silently break when signatures change.
#
#  Any future C++ interface changes must be reflected here to maintain
#  positional consistency and avoid NULL → double coercion errors.
#
#  Wrappers are organized by tier:
#    Tier 1: Core Simulation   - Main sampling entry points (rNormal_reg, etc.)
#    Tier 2: Envelope          - Envelope build/eval, EnvelopeCentering,
#                                rNormalGLM_std, rIndepNormalGammaReg_std
#    Tier 3: Model Utilities   - Standardization
#    Tier 4: OpenCL/GPU        - Kernel loading, GPU diagnostics
# -------------------------------------------------------------------------


# =============================================================================
#  Tier 1: Core Simulation
#  Callers: rNormal_reg, rNormalGamma_reg, rindepNormalGamma_reg, rGamma_reg,
#           rNormalGLM_reg_block
#  User:    All users – primary paths via rglmb, rlmb, glmb, pfamily
# =============================================================================

#' @noRd
#' @keywords internal
.rNormalGLM_cpp <- function(n, y, x, mu, P, offset, wt, dispersion, f2, f3, start, family = "binomial", link = "logit", Gridtype = 2L, n_envopt = -1L, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE) {
  .Call(`_glmbayes_rNormalGLM_cpp_export`, n, y, x, mu, P, offset, wt, dispersion, f2, f3, start, family, link, Gridtype, n_envopt, use_parallel, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.rNormalGLMBlocks_cpp <- function(n, y, x, offset, wt, dispersion, mu, P_blocks, prior_by_block, row_blocks, f2, f3, family = "binomial", link = "logit", Gridtype = 2L, n_envopt = -1L, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE) {
  .Call(`_glmbayes_rNormalGLMBlocks_cpp_export`, n, y, x, offset, wt, dispersion, mu, P_blocks, prior_by_block, row_blocks, f2, f3, family, link, Gridtype, n_envopt, use_parallel, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.rNormalReg_cpp <- function(
    n, y, x, mu, P, offset, wt, dispersion,
    f2, f3, start,
    family = "gaussian",
    link = "identity",
    Gridtype = 2
) {
  .Call(
    "_glmbayes_rNormalReg_cpp_export",
    n, y, x, mu, P, offset, wt, dispersion,
    f2, f3, start,
    family, link, Gridtype
  )
}

#' @noRd
#' @keywords internal
.rIndepNormalGammaReg_cpp <- function(n, y, x, mu, P, offset, wt, shape, rate, max_disp_perc, disp_lower, disp_upper, Gridtype, n_envopt, use_parallel, use_opencl, verbose, progbar) {
  .Call(`_glmbayes_rIndepNormalGammaReg_cpp_export`, n, y, x, mu, P, offset, wt, shape, rate, max_disp_perc, disp_lower, disp_upper, Gridtype, n_envopt, use_parallel, use_opencl, verbose, progbar)
}

#' @noRd
#' @keywords internal
.rNormalGammaReg_cpp <- function(n, y, x, mu, P, offset, wt, shape, rate,
                                 max_disp_perc, disp_lower, disp_upper,
                                 verbose = FALSE) {
  .Call(`_glmbayes_rNormalGammaReg_cpp_export`,
        n, y, x, mu, P, offset, wt, shape, rate,
        max_disp_perc, disp_lower, disp_upper, verbose)
}

#' @noRd
#' @keywords internal
.rGammaGaussian_cpp <- function(n, y, x, beta, wt, alpha, shape, rate,
                                disp_lower = NULL, disp_upper = NULL,
                                verbose = FALSE) {
  .Call(`_glmbayes_rGammaGaussian_cpp_export`,
        n, y, x, beta, wt, alpha, shape, rate,
        disp_lower, disp_upper, verbose)
}

#' @noRd
#' @keywords internal
.rGammaGamma_cpp <- function(n, y, x, beta, wt, alpha, shape, rate,
                             max_disp_perc, disp_lower = NULL,
                             disp_upper = NULL, verbose = FALSE) {
  .Call(`_glmbayes_rGammaGamma_cpp_export`,
        n, y, x, beta, wt, alpha, shape, rate,
        max_disp_perc, disp_lower, disp_upper, verbose)
}


# =============================================================================
#  Tier 2: Envelope & Standardization
#  Callers: EnvelopeSize, EnvelopeBuild, EnvelopeEval, EnvelopeDispersionBuild,
#           EnvelopeOrchestrator, EnvelopeCentering, rNormalGLM_std,
#           rIndepNormalGammaReg_std; EnvelopeSet_* are internal
#  User:    Advanced users – understanding algorithm, custom envelope workflows
# =============================================================================

#' @noRd
#' @keywords internal
.rNormalGLM_std_cpp <- function(n, y, x, mu, P, alpha, wt,
                                f2, Envelope,
                                family, link,
                                progbar = 1L,
                                verbose = FALSE) {
  .Call(`_glmbayes_rNormalGLM_std_cpp_export`,
        n, y, x, mu, P, alpha, wt,
        f2, Envelope,
        family, link,
        progbar, verbose)
}

#' @noRd
#' @keywords internal
.rIndepNormalGammaReg_std_cpp <- function(n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose) {
  .Call(`_glmbayes_rIndepNormalGammaReg_std_cpp_export`, n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose)
}

#' @noRd
#' @keywords internal
.rIndepNormalGammaReg_std_parallel_cpp <- function(n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose) {
  .Call(`_glmbayes_rIndepNormalGammaReg_std_parallel_cpp_export`, n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeCentering_cpp <- function(y, x, mu, P, offset, wt, shape, rate, Gridtype = 2L, verbose = FALSE) {
  .Call(`_glmbayes_EnvelopeCentering_cpp_export`, y, x, mu, P, offset, wt, shape, rate, Gridtype, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeSize_cpp <- function(a, G1, Gridtype, n, n_envopt, use_opencl, verbose) {
  .Call(`_glmbayes_EnvelopeSize_cpp_export`, a, G1, Gridtype, n, n_envopt, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeBuild_cpp <- function(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose) {
  .Call(`_glmbayes_EnvelopeBuild_cpp_export`, bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeBuild_Ind_Normal_Gamma_cpp <- function(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose) {
  .Call(`_glmbayes_EnvelopeBuild_Ind_Normal_Gamma_cpp_export`, bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeEval_cpp <- function(G4, y, x, mu, P, alpha, wt,
                          family, link,
                          use_opencl = FALSE,
                          verbose = FALSE) {
  .Call(`_glmbayes_EnvelopeEval_cpp_export`,
        G4, y, x, mu, P, alpha, wt,
        family, link,
        use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeDispersionBuild_cpp <- function(
    Env,
    Shape,
    Rate,
    P,
    y,
    x,
    alpha,
    n_obs,
    RSS_post,
    RSS_ML,
    mu,
    wt,
    max_disp_perc,
    disp_lower = NULL,
    disp_upper = NULL,
    verbose = FALSE,
    use_parallel = TRUE
) {
  .Call(
    "_glmbayes_EnvelopeDispersionBuild_cpp_export",
    Env,
    Shape,
    Rate,
    P,
    y,
    x,
    alpha,
    n_obs,
    RSS_post,
    RSS_ML,
    mu,
    wt,
    max_disp_perc,
    disp_lower,
    disp_upper,
    verbose,
    use_parallel
  )
}

#' @noRd
#' @keywords internal
.EnvelopeOrchestrator_cpp <- function(bstar2, A, y, x2, mu2, P2, alpha, wt, n, Gridtype, n_envopt, shape, rate, RSS_Post2, RSS_ML, max_disp_perc, disp_lower, disp_upper, use_parallel, use_opencl, verbose) {
  .Call(`_glmbayes_EnvelopeOrchestrator_cpp_export`, bstar2, A, y, x2, mu2, P2, alpha, wt, n, Gridtype, n_envopt, shape, rate, RSS_Post2, RSS_ML, max_disp_perc, disp_lower, disp_upper, use_parallel, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeSet_Grid_cpp <- function(GIndex, cbars, Lint) {
  .Call(`_glmbayes_EnvelopeSet_Grid_cpp_export`, GIndex, cbars, Lint)
}

#' @noRd
#' @keywords internal
.EnvelopeSet_LogP_cpp <- function(logP, NegLL, cbars, G3) {
  .Call(`_glmbayes_EnvelopeSet_LogP_cpp_export`, logP, NegLL, cbars, G3)
}


# =============================================================================
#  Tier 3: Model Utilities
#  Callers: glmb_Standardize_Model
#  User:    Advanced users – model preparation, standardization
# =============================================================================

#' @noRd
#' @keywords internal
.glmb_Standardize_Model_cpp <- function(y, x, P, bstar, A1) {
  .Call(`_glmbayes_glmb_Standardize_Model_cpp_export`, y, x, P, bstar, A1)
}


# =============================================================================
#  Tier 4: OpenCL / GPU
#  Callers: has_opencl, get_opencl_core_count, gpu_names
#  Kernel loading: opencltools (see ?opencltools::load_kernel_source)
# =============================================================================

#' @noRd
#' @keywords internal
.has_opencl_cpp <- function() {
  .Call("_glmbayes_has_opencl_cpp_export")
}

#' @noRd
#' @keywords internal
.get_opencl_core_count_cpp <- function() {
  .Call("_glmbayes_get_opencl_core_count_cpp_export")
}

#' @noRd
#' @keywords internal
.gpu_names_cpp <- function() {
  .Call("_glmbayes_gpu_names_cpp_export")
}


# =============================================================================
#  Phased Out (no R wrappers; C++ exports may still exist for compatibility)
#  - .rss_face_at_disp_cpp, .UB2_cpp
#  - Former RSS/UB2 minimization callbacks; active path uses closed-form C++ bounds
# =============================================================================
