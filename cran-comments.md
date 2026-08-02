# CRAN submission comments — glmbayes 0.9.75

## Summary

**Resubmission** after archival. **glmbayes 0.9.75** follows **0.9.74** (PoCL NOTE
still on **debian-gcc** incoming), **0.9.73** (auto-rejected **2026-08-02**), and
**0.9.72**, archived on **2026-07-31** for a configure policy violation (see
correspondence below).

**Configure policy (0.9.73):** No GitHub install recommendations and no custom Rcpp
version checks. **`configure` / `configure.win`** rely on standard **`LinkingTo:
Rcpp`** and CRAN **Rcpp** only (Windows still passes `-DRCPP_PARALLEL_USE_TBB=1`
and `RcppParallel::RcppParallelLibs()` for TBB linking).

**PoCL-aware OpenCL (0.9.75, Linux/macOS `configure`):** `-DUSE_OPENCL` is defined
only when a **non-PoCL** platform exposes at least one **GPU** device. PoCL-only
stacks on CRAN **debian-gcc** no longer compile OpenCL into **glmbayes**, so
`has_opencl()` is `FALSE` during check and this package should not invoke
`clBuildProgram` / PoCL cache files under `~/.cache/pocl/uncached/tempfile_*`.
**0.9.74** example `NOT_CRAN` gating alone did not clear that NOTE (compile-time
OpenCL was still enabled). Windows **`configure.win`** unchanged.

**Vignette (0.9.75):** Chapter 16 diagnostic / `example(Cleveland)` chunks are
`eval=FALSE` during rebuild (defense in depth).

No other functional or API changes.

OpenCL support remains **optional**.

### Correspondence with CRAN (2026-07-31)

CRAN (Prof Ripley) reported the following **`configure` output** from **0.9.72**
(R 4.7.0):

```
configure: R version: 4.7.0 (svn: 90279)
configure: Rcpp Repository: CRAN
configure: Rcpp Function.h: line with R_getVarEx + R_UnboundValue present = FALSE
configure: WARNING: Rcpp looks like a CRAN install (no GitHub Remote* fields).
configure: WARNING: On R-devel / R >= 4.5, stale CRAN headers can be incompatible
...
```

**0.9.73** implements the fix offered in that reply: all configure logic that
emitted those warnings is removed; **`configure` / `configure.win`** no longer
mention GitHub, `install_github`, Remote* fields, or custom Rcpp version checks.

### Auto-rejection (2026-08-02, **0.9.73** / **0.9.74**)

Incoming pretest reported **2 NOTEs** on Windows and **debian-gcc** (in addition to
incoming feasibility): Windows *checking compiled code* (`'cc' not on the path` /
`ccE`); debian-gcc *for new files in some other directories*
(`~/.cache/pocl/uncached/tempfile_*`). **0.9.75** addresses the PoCL path via
configure (above); the Windows `ccE` NOTE appears environmental (same on
**nmathopencl** / **opencltools** r-devel win-builder).

### CRAN incoming feasibility NOTE

Local **`R CMD check --as-cran`** reports:

```
checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Kjell Nygren <kjell.a.nygren@gmail.com>'

  New submission

  Package was archived on CRAN

  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2026-07-31 for policy violation.

    Recommends packages from GitHub.
```

This NOTE is expected for a resubmission after archival. The archival reason is
addressed in **0.9.73** / **0.9.75** as described above.

### Win-builder (2026-08-02)

**glmbayes 0.9.73** checked on win-builder (reference for **0.9.75** resubmission):

- **R 4.6.1 (release, ucrt):** `Status: 1 NOTE` (install ~295s, check ~247s).
  The NOTE is *checking CRAN incoming feasibility* only: new submission, package
  was archived on CRAN, and the repository `X-CRAN-Comment` for the prior archival
  (GitHub Rcpp policy violation). This is the **resubmission context** addressed
  in **0.9.73** (configure cleanup above); install and build succeed.

- **R Under development (2026-07-30 r90327, ucrt):** `Status: 2 NOTEs`
  (install ~312s, check ~244s):

  1. *checking CRAN incoming feasibility* — as on release, plus a standard
     “possibly misspelled words in DESCRIPTION” line (`Nygren`, `R's`,
     `subgradients`; maintainer name and technical terms).

  2. *checking compiled code*:

```
Error in ccE(...): 'cc' is not on the path
```

  Install and check complete with **2 NOTEs** on r-devel. R code uses **Rcpp /
  `.Call()`** only (no `.C()` / `.Fortran()`). The *checking compiled code* NOTE
  comes from **`tools:::ccE`** during R's DLL API scan when `cc` is missing from
  PATH—we do not believe this reflects a defect in **glmbayes** source (same
  class of NOTE on **nmathopencl 0.8.4** win-builder r-devel, 2026-08-02).

## Test environments

* **win-builder (2026-08-02):** release — 1 NOTE (incoming feasibility /
  resubmission after archival only); r-devel — 2 NOTEs (incoming feasibility +
  compiled code `cc` / `ccE` as above).

* **local Windows 10, `R CMD check --as-cran`:** 1 NOTE (incoming feasibility /
  resubmission after archival only).

* **R-universe (2026-08-02, glmbayes 0.9.73 / prior):** passes cleanly (0 errors,
  warnings, or NOTEs) on all built platforms except **wasm** (not supported;
  package requires native compiled C++/OpenCL stack).

---
_This file is listed in `.Rbuildignore` and is not included in the built source
tarball. When submitting, paste the content above into the “Optional comments”
field on the CRAN submission form at_
https://cran.r-project.org/submit.html
