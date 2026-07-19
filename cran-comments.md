# CRAN submission comments — glmbayes 0.9.7

## Summary

This is an update from glmbayes 0.9.6 (current CRAN release) to 0.9.7.

Main changes:

* **Hard dependency on `nmathopencl`:** OpenCL statistical kernels are now
  supplied by the **`nmathopencl`** package (on CRAN, including Windows
  binaries). **glmbayes**-specific envelope/likelihood kernels remain under
  `inst/cl/`. CPU nmath uses R's libR via `<Rmath.h>`; vendored CPU `nmath`
  sources were removed in this development cycle.

* **`insight` and `bayestestR` integration:** New S3 methods for
  `model_info()`, `get_parameters()`, `get_priors()`, and related **`insight`**
  accessors on **`glmb`** fits, plus **`simulate_prior()`**, **`check_prior()`**,
  and **`describe_prior()`** for prior checking. Both packages move from
  `Suggests` to `Imports`.

* **Prior simulation (`pfun`):** Each **`pfamily()`** constructor stores a
  matching prior simulation function; **`simulate_prior.glmb()`** delegates to
  it (mirroring how **`simfun`** drives posterior sampling).

* **Bug fixes:**
  - **`print.lmb()`** again shows a concise **`Call:`** line instead of
    deparsing the full **`pfamily`** object (regression from 0.9.6 multi-response
    work).
  - **`simulate_prior.glmb()`** for **`dIndependent_Normal_Gamma()`** now draws
    dispersion from the same truncated inverse-gamma the sampler uses (via
    bounds written back into **`pfamily`** after fitting), not the untruncated
    prior.
  - Documentation/roxygen cleanup (e.g. stray content in **`rglmb.Rd`**).

OpenCL support remains **optional** (`SystemRequirements: Optional OpenCL support`).
CRAN binary builds compile without OpenCL; GPU acceleration requires a source
install with OpenCL available at compile time. All OpenCL-specific **testthat**
tests call `skip_on_cran()` (in addition to `skip_if_no_opencl()`).

## Test environments

* local Windows 10, `R CMD check --as-cran`: 0 errors | 0 warnings | 0 notes

* win-builder (CRAN): 0 errors | 0 warnings | 0 notes on **r-devel**, **r-release**,
  and **r-oldrel** (to be confirmed on resubmission)

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
