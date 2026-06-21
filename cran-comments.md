# CRAN submission comments — glmbayes 0.9.6

## Summary

This is an update from glmbayes 0.9.5 (current CRAN release) to 0.9.6.

Main changes:

* **Multi-response `lmb()`:** unified interface for univariate and multivariate
  responses (mirroring `lm()`); multi-column formulas return an `mlmb` object with
  print/summary/coefficient methods.

* **Conjugate GLM priors:** closed-form IID sampling for intercept-only models with
  identity links — `dBeta()` / `rBeta_reg()` (Beta–Binomial) and
  `dGamma(Inv_Dispersion = FALSE)` / `rGamma_Conjugate_reg()` (Gamma–Poisson and
  Gamma–Gamma rate priors); `Prior_Setup()` calibrates conjugate hyperparameters.

* **Vignettes:** reworked Chapter 00 roadmap; Chapter 02 split into conceptual
  introduction plus worked examples (Chapter 02-S01–S05); optional Bayes Rules!
  and LearnBayes companion appendices.

* **New import:** `opencltools (>= 0.8.1)` for shared OpenCL host/runtime
  discovery and diagnostics. Package-specific entry points (`has_opencl()`,
  `diagnose_glmbayes()`) remain in glmbayes and report compile-time OpenCL status
  for this build.

* **Bug fix / guard:** `rindepNormalGamma_reg()` now rejects prior-dominated
  calls where the Gamma component carries more effective prior observations than
  the data supply (`n_prior > sum(weights)`), preventing silent envelope degradation.

There is **no** dependency on unreleased packages (`nmathopencl`, `glmbayesCore`).
OpenCL statistical kernels remain bundled in this package under `inst/cl/`.

OpenCL support is **optional** (`SystemRequirements: Optional OpenCL support`).
CRAN binary builds are expected to compile without OpenCL; GPU acceleration is
available only when the package is built from source on a system with OpenCL
headers and a working runtime. All OpenCL-specific **testthat** tests call
`skip_on_cran()` (in addition to `skip_if_no_opencl()`), consistent with 0.9.5.

Suggested packages (`bayesrules`, `LearnBayes`, `bayestestR`, etc.) are used only
in vignettes/examples and remain in `Suggests`.

## Test environments

* local Windows 10, `R CMD check --as-cran`: 0 errors | 0 warnings | 0 notes

* win-builder (CRAN): 0 errors | 0 warnings | 0 notes on **r-devel**, **r-release**,
  and **r-oldrel**

* R-universe (https://knygren.r-universe.dev/glmbayes): checks passed cleanly across
  the build matrix

## Check notes addressed in this release

* OpenCL / parallel tests skipped on CRAN to avoid CPU-time vs elapsed-time
  NOTEs on check farms without GPU.
* Spelling check uses `skip_on_cran = TRUE` in `tests/spelling.R`.
* OpenCL-dependent examples are wrapped in `\donttest{}` where appropriate.

---
_This file is listed in `.Rbuildignore` and is not included in the built source
tarball. When submitting, paste the **Summary** and **Test environments**
sections above into the “Optional comments” field on the CRAN submission form at_
https://cran.r-project.org/submit.html
