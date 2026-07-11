#include "RcppArmadillo.h"
#include "Envelopefuncs.h"
#include "openclPort.h"
#include "opencl.h"
#include "simfuncs.h"

using namespace openclPort;
using namespace glmbayes::env;
using namespace glmbayes::sim;


// -----------------------------------------------------------------------------
// Wrapper organization mirrors R/rcpp_wrappers.R:
//   Tier 1: Core Simulation   - Main sampling entry points
//   Tier 2: Envelope          - Build/eval; used by rNormalGLM
//   Tier 3: Indep NG std      - Split workflow samplers
//   Tier 4: Model Utilities   - Standardization
//   Tier 5: OpenCL/GPU        - Kernel loading, diagnostics
//   Phased out: rss_face_at_disp, UB2 (no R wrappers)
// -----------------------------------------------------------------------------


// =============================================================================
// Tier 1: Core Simulation
// Callers: rNormal_reg, rNormalGamma_reg, rindepNormalGamma_reg, rGamma_reg,
//          rNormalGLM_reg_block (via .rNormalGLMBlocks_cpp when wired)
// User:    All users - primary paths via rglmb, rlmb, glmb, pfamily
// =============================================================================

// [[Rcpp::export]]
Rcpp::List rNormalGLM_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double dispersion,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::NumericVector& start,
    const std::string& family = "binomial",
    const std::string& link   = "logit",
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
) {
  return rNormalGLM(
    n, y, x, mu, P, offset, wt,
    dispersion,
    f2, f3, start,
    family, link, Gridtype,
    n_envopt, use_parallel, use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rNormalGLMBlocks_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    const Rcpp::NumericVector& dispersion,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::List& P_blocks,
    bool prior_by_block,
    const Rcpp::List& row_blocks,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const std::string& family = "binomial",
    const std::string& link   = "logit",
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
) {
  return rNormalGLMBlocks(
    n, y, x, offset, wt,
    dispersion,
    mu, P_blocks, prior_by_block, row_blocks,
    f2, f3,
    family, link, Gridtype,
    n_envopt, use_parallel, use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rNormalReg_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double dispersion,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::NumericVector& start,
    const std::string& family = "gaussian",
    const std::string& link   = "identity",
    int Gridtype = 2
) {
  return rNormalReg(
    n, y, x, mu, P, offset, wt,
    dispersion, f2, f3, start,
    family, link, Gridtype
  );
}

// [[Rcpp::export]]
Rcpp::List rIndepNormalGammaReg_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<Rcpp::NumericVector> disp_lower,
    Rcpp::Nullable<Rcpp::NumericVector> disp_upper,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    bool progbar
) {
  return rIndepNormalGammaReg(
    n,
    y,
    x,
    mu,
    P,
    offset,
    wt,
    shape,
    rate,
    max_disp_perc,
    disp_lower,
    disp_upper,
    Gridtype,
    n_envopt,
    use_parallel,
    use_opencl,
    verbose,
    progbar
  );
}

// [[Rcpp::export]]
Rcpp::List rNormalGammaReg_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double shape,
    double rate,
    Rcpp::Nullable<double> max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    bool verbose = false
) {
  return glmbayes::sim::rNormalGammaReg(
    n,
    y,
    x,
    mu,
    P,
    offset,
    wt,
    shape,
    rate,
    max_disp_perc,
    disp_lower,
    disp_upper,
    verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rGammaGaussian_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& wt,
    const Rcpp::NumericVector& alpha,
    double shape,
    double rate,
    Rcpp::Nullable<double> disp_lower = R_NilValue,
    Rcpp::Nullable<double> disp_upper = R_NilValue,
    bool verbose = false
) {
  return glmbayes::sim::rGammaGaussian(
    n, y, x, beta, wt, alpha,
    shape, rate,
    disp_lower, disp_upper,
    verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rGammaGamma_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& wt,
    const Rcpp::NumericVector& alpha,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower = R_NilValue,
    Rcpp::Nullable<double> disp_upper = R_NilValue,
    bool verbose = false
) {
  return glmbayes::sim::rGammaGamma(
    n, y, x, beta, wt, alpha,
    shape, rate, max_disp_perc,
    disp_lower, disp_upper,
    verbose
  );
}


// =============================================================================
// Tier 2: Envelope & Standardization
// Callers: EnvelopeSize, EnvelopeBuild, EnvelopeEval, EnvelopeDispersionBuild,
//          EnvelopeOrchestrator, rNormalGLM_std; EnvelopeSet_* are internal
// User:    Advanced users - understanding algorithm, custom envelope workflows
// =============================================================================

// [[Rcpp::export]]
Rcpp::List rNormalGLM_std_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2,
    const Rcpp::List& Envelope,
    const Rcpp::CharacterVector& family,
    const Rcpp::CharacterVector& link,
    int progbar = 1,
    bool verbose = false
) {
  return rNormalGLM_std(
    n, y, x, mu, P, alpha, wt,
    f2, Envelope, family, link,
    progbar, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeSize_cpp_export(
    const arma::vec& a,
    const Rcpp::NumericMatrix& G1,
    int Gridtype,
    int n,
    int n_envopt,
    bool use_opencl,
    bool verbose
) {
  return glmbayes::env::EnvelopeSize(
    a, G1, Gridtype, n, n_envopt, use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeBuild_cpp_export(
    Rcpp::NumericVector bStar,
    Rcpp::NumericMatrix A,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    std::string family,
    std::string link,
    int Gridtype,
    int n,
    int n_envopt,
    bool sortgrid,
    bool use_opencl,
    bool verbose
) {
  return glmbayes::env::EnvelopeBuild(
    bStar, A, y, x, mu, P, alpha, wt,
    family, link, Gridtype, n, n_envopt,
    sortgrid, use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeBuild_Ind_Normal_Gamma_cpp_export(
    const Rcpp::NumericVector& bStar,
    const Rcpp::NumericMatrix& A,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const std::string& family = "binomial",
    const std::string& link   = "logit",
    int Gridtype              = 2,
    int n                     = 1,
    int n_envopt              = -1,
    bool sortgrid             = false,
    bool use_opencl           = false,
    bool verbose              = false
) {
  return EnvelopeBuild_Ind_Normal_Gamma(
    bStar, A, y, x, mu, P, alpha, wt,
    family, link,
    Gridtype, n, n_envopt,
    sortgrid, use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeEval_cpp_export(
    const Rcpp::NumericMatrix& G4,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const std::string& family,
    const std::string& link,
    bool use_opencl = false,
    bool verbose = false
) {
  return EnvelopeEval(
    G4, y, x, mu, P, alpha, wt,
    family, link,
    use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeDispersionBuild_cpp_export(
    const Rcpp::List& Env,
    double Shape,
    double Rate,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    int n_obs,
    double RSS_post,
    double RSS_ML,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericVector& wt,
    double max_disp_perc = 0.99,
    Rcpp::Nullable<double> disp_lower = R_NilValue,
    Rcpp::Nullable<double> disp_upper = R_NilValue,
    bool verbose = false,
    bool use_parallel = true
) {
  return EnvelopeDispersionBuild(
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
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeOrchestrator_cpp_export(
    const Rcpp::NumericVector& bstar2,
    const Rcpp::NumericMatrix& A,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x2,
    const Rcpp::NumericMatrix& mu2,
    const Rcpp::NumericMatrix& P2,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    int n,
    int Gridtype,
    Rcpp::Nullable<int> n_envopt,
    double shape,
    double rate,
    double RSS_Post2,
    double RSS_ML,
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    bool use_parallel,
    bool use_opencl,
    bool verbose
) {
  return EnvelopeOrchestrator(
    bstar2,
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
    use_parallel,
    use_opencl,
    verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeCentering_cpp_export(
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double shape,
    double rate,
    int Gridtype = 2,
    bool verbose = false
) {
  return glmbayes::env::EnvelopeCentering(
    y, x, mu, P, offset, wt,
    shape, rate,
    Gridtype, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeSet_Grid_cpp_export(
    const Rcpp::NumericMatrix& GIndex,
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericMatrix& Lint
) {
  return EnvelopeSet_Grid(
    GIndex,
    cbars,
    Lint
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeSet_LogP_cpp_export(
    const Rcpp::NumericMatrix& logP,
    const Rcpp::NumericVector& NegLL,
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericMatrix& G3
) {
  return EnvelopeSet_LogP(
    logP,
    NegLL,
    cbars,
    G3
  );
}


// =============================================================================
// Tier 3: Standardized Samplers (Indep Normal-Gamma)
// Callers: C++ only (rIndepNormalGammaReg); R wrappers for custom split workflow
// User:    Advanced / developers - after EnvelopeOrchestrator, sample separately
// =============================================================================

// [[Rcpp::export]]
Rcpp::List rIndepNormalGammaReg_std_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2,
    const Rcpp::List& Envelope,
    const Rcpp::List& gamma_list,
    const Rcpp::List& UB_list,
    const Rcpp::CharacterVector& family,
    const Rcpp::CharacterVector& link,
    bool progbar = true,
    bool verbose = false
) {
  return rIndepNormalGammaReg_std(
    n, y, x, mu, P, alpha, wt,
    f2, Envelope, gamma_list, UB_list,
    family, link, progbar, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rIndepNormalGammaReg_std_parallel_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2,
    const Rcpp::List& Envelope,
    const Rcpp::List& gamma_list,
    const Rcpp::List& UB_list,
    const Rcpp::CharacterVector& family,
    const Rcpp::CharacterVector& link,
    bool progbar = true,
    bool verbose = false
) {
  return rIndepNormalGammaReg_std_parallel(
    n, y, x, mu, P, alpha, wt,
    f2, Envelope, gamma_list, UB_list,
    family, link, progbar, verbose
  );
}


// =============================================================================
// Tier 4: Model Utilities
// Callers: glmb_Standardize_Model
// User:    Advanced users - model preparation, standardization
// =============================================================================

// [[Rcpp::export]]
Rcpp::List glmb_Standardize_Model_cpp_export(
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericMatrix& bstar,
    const Rcpp::NumericMatrix& A1
) {
  return glmb_Standardize_Model(
    y, x, P, bstar, A1
  );
}


// =============================================================================
// Tier 5: OpenCL / GPU
// Callers: has_opencl, get_opencl_core_count, gpu_names
// User:    Advanced users - GPU diagnostics for use_opencl
// Kernel loading: opencltools (see ?opencltools::load_kernel_source)
// =============================================================================

// [[Rcpp::export]]
bool has_opencl_cpp_export() {
  return has_opencl();
}

// [[Rcpp::export]]
int get_opencl_core_count_cpp_export() {
  return get_opencl_core_count();
}

// [[Rcpp::export]]
Rcpp::CharacterVector gpu_names_cpp_export() {
  return gpu_names();
}

#ifdef USE_OPENCL
// [[Rcpp::export]]
Rcpp::List debug_likelihood_program_cpp_export(
    std::string family,
    std::string link)
{
  using glmbayes::opencl::load_likelihood_subgradient_program;

  const std::string program =
      load_likelihood_subgradient_program(family, link, "glmbayes");

  return Rcpp::List::create(
      Rcpp::Named("program") = program,
      Rcpp::Named("nchar") = static_cast<int>(program.size()));
}
#endif


// =============================================================================
// Phased Out (no R wrappers; C++ exports commented out)
// - rss_face_at_disp, UB2: former RSS/UB2 minimization; active path uses
//   closed-form C++ bounds.
//
// To fully remove: delete this block, then (1) remove *.o from src/,
// (2) uninstall old glmbayes, (3) Rcpp::compileAttributes(),
// (4) devtools::document(), (5) devtools::install().
// =============================================================================

/*
// [[Rcpp::export]]
double rss_face_at_disp_cpp_export(
    double dispersion,
    const Rcpp::List& cache,
    const Rcpp::NumericVector& cbars_j,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt
) {
  return rss_face_at_disp(
    dispersion,
    cache,
    cbars_j,
    y,
    x,
    alpha,
    wt
  );
}

// [[Rcpp::export]]
double UB2_cpp_export(
    double dispersion,
    const Rcpp::List& cache,
    const Rcpp::NumericVector& cbars_j,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    double rss_min_global
) {
  return UB2(
    dispersion,
    cache,
    cbars_j,
    y,
    x,
    alpha,
    wt,
    rss_min_global
  );
}
*/
