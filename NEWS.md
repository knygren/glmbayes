# glmbayes 0.9.76

## Under development

* Stage 0 migration prep: file-vs-file inventory vs **glmbayesCore** in
  `data-raw/CORE_MIGRATION_DIFF.md`; Stage 1 will Import **glmbayesCore** (>= 0.5.4).

# glmbayes 0.9.75

## Bug fixes

* **Configure (Linux/macOS):** `-DUSE_OPENCL` is set only when a **non-PoCL**
  OpenCL platform exposes at least one **GPU** device. PoCL-only ICD stacks (typical
  on CRAN **debian-gcc** incoming) no longer enable OpenCL at compile time, so
  `has_opencl()` is `FALSE` during check and kernel builds do not touch
  `~/.cache/pocl` (addresses incoming NOTE on `tempfile_*` there). Real GPU installs
  (NVIDIA/AMD/Intel) are unchanged.

* **Vignette Chapter 16:** `diagnose_glmbayes()`, `opencltools::diagnose_glmbayes()`,
  and `example(Cleveland)` chunks are `eval=FALSE` during vignette rebuild (defense
  in depth alongside configure).

# glmbayes 0.9.74

## Bug fixes

* **OpenCL examples:** Cleveland, Boston_centered, and `gpu_diagnostics` kernel-load
  examples run OpenCL code only when `NOT_CRAN="true"` and `has_opencl()` (same
  policy as OpenCL **testthat** skips on CRAN). During CRAN checks they print a
  skip message instead of calling **opencltools** / GPU paths (avoids PoCL cache
  NOTEs on incoming Linux builders when OpenCL is compiled in).

# glmbayes 0.9.73

## Bug fixes and Removal of any external recommendation for Github install

* **Configure:** Removed **`tools/rcpp_include.R`** Rcpp header probing, Function.h
  branch flags, and **`glmbayes_getRegisteredNamespace`** shim on Unix and Windows
  (rely on **`LinkingTo: Rcpp`** and CRAN **Rcpp** instead).

# glmbayes 0.9.72

## Bug fixes

* **Examples:** `Prior_Setup`, `pfamily`, and `simfuncs` examples use
  `use_parallel = FALSE` during CRAN checks.

# glmbayes 0.9.71

## Bug fixes

* **`configure.win`:** Windows builds now set `-DRCPP_PARALLEL_USE_TBB=1` and
  append `RcppParallel::RcppParallelLibs()` to **`PKG_LIBS`**, matching updated
  CRAN/R-Universe **`RcppParallel`** TBB linking requirements.

# glmbayes 0.9.7

## Highlights

* **CRAN dependency on `nmathopencl`:** OpenCL statistical kernels are now
  supplied by the **`nmathopencl`** package (hard **`Imports`** dependency),
  including Windows binaries on CRAN — unblocking GPU/OpenCL builds on Windows
  without vendored nmath in **glmbayes**.
* **Bug fixes** (see below): **`print.lmb()`** call display, truncated
  **`dIndependent_Normal_Gamma`** dispersion prior simulation for
  **`simulate_prior()`**, and related documentation/roxygen cleanup.

## insight integration

* Added S3 methods for **`insight`** accessors on **`glmb`** fits
  (`model_info`, `get_parameters`, `find_parameters`, `find_algorithm`,
  `get_data`, `get_priors`). **`insight`** moves from `Suggests` to
  **`Imports`** and key generics are re-exported. **`get_priors()`** returns
  **`pfamily(model)`** (full prior specification, including complete
  **`Sigma`**) rather than a marginal-only table.

## bayestestR prior-checking integration

* Added **`simulate_prior()`**, **`check_prior()`**, and **`describe_prior()`**
  methods for **`glmb`** fits; **`bayestestR`** moves from `Suggests` to
  **`Imports`**. Each **`pfamily()`** stores a matching **`pfun`** (like
  **`simfun`**); **`simulate_prior.glmb()`** calls it, including the truncated
  inverse-gamma dispersion prior used by **`dIndependent_Normal_Gamma()`**
  fits.

## Bug fixes

* **`print.lmb()`:** Single-response **`lmb()`** fits again show a concise
  **`Call:`** line (formula, `n`, `data`, and related arguments) instead of
  deparsing the full **`pfamily`** object, including embedded **`simfun`**
  source. When multi-response **`lmb()`** / **`mlmb`** support was added in
  0.9.6, **`.mlmb_lmb_display_call()`** omitted **`pfamily`** from the stored
  call for each block fit only; univariate **`lmb()`** still used the internal
  **`.uni_lmb()`** matched call via **`do.call()`**, which inlined the evaluated
  **`pfamily`**. Univariate fits now use the same display-call helper as
  **`mlmb`**.

## Independent Normal-Gamma simulation

* Updates to independent normal gamma simulation to better handle
  non-isotropic priors with highly differentiated implied pweights across
  dimensions.

## CPU nmath phase-out

* Removed vendored CPU R Mathlib sources; CPU routines in C++ now use R's
  libR via `<Rmath.h>`.

## OpenCL kernel loading

* OpenCL nmath comes from **`nmathopencl`**; **glmbayes** likelihood/envelope
  kernels remain under `inst/cl/`, loaded via **opencltools**.

# glmbayes 0.9.6

## Highlights

* **Multi-response `lmb()`:** **`lmb()`** now handles both univariate and
  multivariate responses with a single unified interface, mirroring the behaviour
  of R's **`lm()`**. When the response has a single column the result is an
  **`lmb`** object (unchanged from prior releases). When the formula specifies
  multiple response columns (e.g. `cbind(y1, y2) ~ x`), **`lmb()`** fits a
  separate Bayesian linear model per response column and returns a named list
  with class **`mlmb`**. For the multi-response case, **`pfamily`** must be a
  list of **`pfamily`** objects with exactly one entry per response column;
  passing a single **`pfamily`** object is an error. Summary, print, and
  coefficient methods for **`mlmb`** objects are included.

* **Conjugate GLM priors (Poisson, binomial, Gamma):** New closed-form IID
  sampling paths for intercept-only models with identity links. **`dBeta()`**
  with **`rBeta_reg()`** supports Beta–Binomial(identity) conjugate updates;
  **`dGamma(Inv_Dispersion = FALSE)`** with **`rGamma_Conjugate_reg()`**
  supports Gamma–Poisson(identity) and Gamma–Gamma(identity) rate priors.
  **`Prior_Setup()`** can calibrate conjugate hyperparameters for these
  families (weighted Poisson rate and binomial probability defaults). See
  **`?dBeta`**, **`?dGamma`**, and the Chapter 02 / Chapter 07–11 vignettes.

* **Vignette structure:** Reworked **Chapter 00** as a roadmap across five
  main parts plus technical appendices. **Chapter 02** is now a conceptual
  introduction to single-parameter conjugacy; worked examples move to
  **Chapter 02-S01** through **Chapter 02-S05** (Beta–Binomial, Normal–Normal,
  Gamma–Poisson, exposure-weighted Poisson, and related topics). A **Companion
  textbooks** section in Chapter 00 indexes optional Bayes Rules! and `LearnBayes`
  appendices tied to the main GLM chapters.

* **`opencltools` import:** Core host/runtime OpenCL discovery and diagnostics
  (`detect_*`, PATH helpers, environment checks) now live in the **`opencltools`**
  package (`Imports`, >= 0.8.0). **glmbayes** keeps package-specific entry
  points (`has_opencl()`, `diagnose_glmbayes()`) that report compile-time
  OpenCL status for this build while delegating shared GPU/runtime checks—reducing
  duplicated maintenance in **glmbayes**.

* **Bayes Rules! companion examples:** Optional vignette appendices reproduce
  book datasets and published posterior summaries using **`lmb()`**, **`glmb()`**,
  **`Prior_Setup()`**, and **`dNormal()`** (suggested package **`bayesrules`** for
  data only). Coverage includes **`bikes`** (Ch. 03), **`weather_perth`** (Ch. 08–09),
  **`equality_index`** (Ch. 10), Gamma–Poisson conjugacy (Ch. 02-S04), and a
  scope note for Gamma regression (Ch. 11). Comparison tables use **printed book
  values**, not live **`rstanarm`** fits. See **Chapter 00** § Companion textbooks.

* **`LearnBayes` examples:** **Chapter 02-S04**, Appendix A, maps the
  **`hearttransplants`** example from Albert (2009) / `LearnBayes` (exposure-weighted
  Gamma–Poisson conjugacy) to **`glmb()`** with analytic Albert posteriors for
  verification (suggested package **`LearnBayes`**).

## Other changes

* **Prior-vs-data guard for `dIndependent_Normal_Gamma` sampling:**
  **`rindepNormalGamma_reg()`** now rejects calls where the Gamma (precision)
  part of the prior carries more effective prior observations than the data
  supply: inverting the `Prior_Setup()` calibration
  `shape = (n_prior + 1 + p)/2`, sampling requires
  `n_prior <= n_w = sum(weights)` (equivalently a prior weight
  `pwt <= 0.5`). Rationale: the dispersion envelope caps its log-tilt at
  `n_w/2` - the *data* contribution to the posterior Gamma shape (Remark
  4.1.3 of the ING vignette) - a strengthening of the validity condition
  `lm_log2 < shape2` that presumes a likelihood-dominated regime.
  Prior-dominated calls could previously bind that cap on every envelope
  build (console `UB3A mean slope` warnings) and silently degrade the
  envelope. Note that `n_prior` here is the effective sample size of the
  Gamma component specifically; under the `Prior_Setup()` calibration the
  Gamma and coefficient parts share a common `n_prior`, so the two are not
  fully independent.

* Expanded **testthat** coverage for **`dBeta()`** / binomial(identity) conjugate
  paths and related **`glmb()`** integration.

# glmbayes 0.9.5

* **Tests / CRAN:** All **OpenCL**-specific **testthat** blocks now call
  **`skip_on_cran()`** (in addition to **`skip_if_no_opencl()`**), consistent
  with existing Boston/Cleveland OpenCL tests. OpenCL coverage remains for local
  checks and source builds with OpenCL; CRAN checks avoid parallel/GPU-heavy
  tests that could trigger **CPU time vs elapsed time** NOTES.

# glmbayes 0.9.4

* **Vignettes:** A vignette that previously used the `notangle` engine now
  uses the standard R Markdown vignette machinery (`knitr` /
  `rmarkdown::html_vignette`), so builds align with CRAN expectations and
  vignette index ordering should be consistent with the rest of the package.

* **OpenCL sources (`inst/cl`):** Removed unused or superseded material,
  consolidated kernels and library fragments, and aligned `.cl` layout and
  dependency tagging with the conventions used in 'openclport' and
  'nmathopencl' (prelude, shims, `nmath/` stems, family kernels under
  `src/`). See `inst/cl/README.md` for how the assembled program is stitched.

* **OpenCL program assembly:** Reworked loading so the full OpenCL program is
  built from explicit fragments (global header, `nmath` closure, family/link
  kernels) rather than ad hoc concatenation—clearer ownership of what enters
  GPU compilation and easier parity with CPU paths.

* **Tests:** Added and expanded **testthat** coverage aimed at OpenCL code
  paths (including binomial examples that exercise GPU envelope evaluation),
  complementing existing Cleveland-style checks.

* **Bug fix — binomial OpenCL:** Binomial `f2_f3` OpenCL kernels now evaluate
  the data log-likelihood with the same **proportion × trial-count**
  semantics as **`dbinom_glmb`** on the CPU (`round` successes and trials,
  clamped probability). This fixes envelope / PLSD failures for aggregated
  binomial data (e.g. `cbind(successes, failures)` / `MASS::menarche`) where
  the previous kernels treated **`y`** like a raw success count.

# glmbayes 0.9.3

* Published on CRAN.
* Version bump in response to CRAN resubmission feedback.

# glmbayes 0.9.2

* Version bump in preparation for resubmission incorporating CRAN review feedback.

# glmbayes 0.9.1

* Wrapped OpenCL-dependent examples in `\donttest{}` for CRAN compliance.
* Reduced iteration counts in rlmb Gibbs sampler example to stay within
  CRAN example time limits on slower check machines.

# glmbayes 0.9.0

First CRAN submission. This release is a stable pre-release with a
near-complete feature set relative to earlier development builds.

## Highlights

### Bayesian Generalized Linear (glmb) and Linear (lmb) modeling functions:

  `glmb()` is a Bayesian analog for the classical `glm()` function while
  `lmb()` covers Gaussian models. Calls largely mirror those for the 
  classical functions but leverage pfamilies for prior specifications.
  Method functions largely mirror those for the classical functions. 
  Samples generated by the functions are largely iid samples 
  (no MCMC convergence dignostics are needed).

### Implemented Likelihood families/ link functions:
   
  Most of the families implemented in the `glm()` function are also implemented 
  in the `glmb()` function (the `lmb()` function covers only gaussian() families). 
  Link functions that lead to log-concave likelihood functions are generally 
  implemented.  Specifically, we have the following:
  
  **Supported likelihoods:** gaussian (identity), Poisson / quasi-Poisson
  (log), binomial / quasi-binomial (logit, probit, cloglog), Gamma (log).

### Prior Family functions:

 `pfamily` constructors are used to specify priors and play the same
  kind of role for the prior specifications as `family` constructors 
  and `link` functions play for the likelihoods. Specifically, we
  have the following:

  **Supported Priors:** Normal (all families/links), Normal–Gamma and 
  independent Normal–Gamma (gaussian families), and Gamma-on-precision 
  (gaussian and Gamma families).
  
### Prior_Setup function:
 
  The package comes with a convenient `Prior_Setup()` function that provides 
  default prior input parameters for each of the implemented models. Basic calls
  (without tailoring) mirror traditional calls to the `glmb()` and `lmb()`
  functions respectively and only require the user to provide the model formula
  and (if not the gaussian family) the family/link function. 
  
  The function can also be used to easily adjust prior specifications 
  (see documentation for details).
  
### Extensive Method functions:
  
  The package comes with extensive method functions that mirror those 
  for the classical functions.  These include dedicated `print()`,
  `summary()`, `predict()` and `simulate()` functions.

### Lower Level Modeling functions:

  The package comes with lower level modeling/simulation functions
  that advanced users can use to implement block Gibbs samplers. These
  generally come with less overhead than the `glmb()` and `lmb()` functions 
  and are called internally by the the higher level modeling functions.

### RcppParallel and OpenCL GPU Acceleration Implementations
  
  Some of the simulation functions comes with use_parallel and use_opencl options
  that speed up simulation for higher dimensional models.
  
### Extensive help files, vignettes, examples and demos

  The package also comes with extensive help files for the varios functions 
  that are complemented with a rich set of vignettes. A large number of 
  examples and demos are also availabel (see the READM.md file for a sample).

---

## Earlier development history (0.1.x series)

The notes below summarize major work during the initial development series
before the 0.9.0 pre-release.

### OpenCL and GPU acceleration

- Completed the OpenCL-based grid construction framework for large models.
- Added GPU-aware envelope sizing and improved OpenCL failure handling.
- Introduced diagnostic utilities to assess OpenCL availability and
  performance.
- Improved configure scripts to detect OpenCL and provide informative
  messages.
- Expanded OpenCL documentation and added a dedicated vignette chapter.

### Parallel CPU sampling (RcppParallel)

- Enabled parallel envelope construction and parallel iid sampling.
- Added pilot functions for large-dimension grid estimation.
- Implemented thread-safe parallel sampling for independent normal-gamma
  models.

### Core statistical improvements

- Migrated to an improved independent normal-gamma simulation algorithm.
- Added theoretical derivations for independent normal-gamma regression.
- Improved UB2 and RSS minimization routines, including scaling corrections.
- Enhanced `Prior_Setup()` to support family-specific prior construction.
- Added dedicated envelope evaluation and sizing functions.

### Package infrastructure

- Significant cleanup to remove NOTES and improve CRAN readiness.
- Improved configure and Makevars files for portability.
- Added testthat tests, including OpenCL-specific tests.
- Consolidated envelope-building functions into a cleaner structure.

### Documentation

- Major updates to README and package-level documentation.
- Added multiple new vignettes and expanded existing ones.
- Improved examples for `lmb()`, `rlmb()`, and OpenCL models.

### Bug fixes (0.1.x era)

- Corrected scaling in UB2 minimization.
- Improved error handling for missing OpenCL functionality.
- Fixed various small issues uncovered during parallelization work.
