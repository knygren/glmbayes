# CRAN submission comments — glmbayes 0.9.73

## Summary

**Resubmission** after archival. **glmbayes 0.9.73** follows **0.9.72**, archived on
**2026-07-31** for a configure policy violation (see correspondence below).

**This submission removes the flagged behaviour entirely.** There are no GitHub
install recommendations and no custom Rcpp version checks.

* **`configure` / `configure.win`:** Removed custom Rcpp header probing
  (`tools/rcpp_include.R`, Function.h patches, `glmbayes_getRegisteredNamespace`
  shim). Builds rely on standard **`LinkingTo: Rcpp`** and CRAN **Rcpp** only.
  Windows still passes `-DRCPP_PARALLEL_USE_TBB=1` and
  `RcppParallel::RcppParallelLibs()` for TBB linking.

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
configure: WARNING: with R (e.g. R_NamespaceRegistry). Consider
configure: WARNING: remotes::install_github("RcppCore/Rcpp") or ensure
configure: WARNING: install_github actually replaced the library.
```

CRAN noted that recommending installation from GitHub violates CRAN policy (CRAN
versions of packages must work with current CRAN/Bioconductor dependencies and
must not recommend development versions from other repositories). Repository note:

> X-CRAN-Comment: Archived on 2026-07-31 for policy violation. Recommends packages from GitHub.

Maintainer reply (2026-07-31): the GitHub suggestion was **not** intended as
ongoing user guidance. It was temporary build-time logic for an R 4.6.0
pre-release / Rcpp compatibility window (R 4.6.0, SVN < 89746; see
[RcppCore/Rcpp#1473](https://github.com/RcppCore/Rcpp/issues/1473)). By the
final **R 4.6** release and current CRAN **Rcpp**, GitHub installs were no
longer required; the configure messages should have been removed before
submission. Apologies for the policy violation.

**0.9.73** implements the fix offered in that reply: all configure logic that
emitted those warnings is removed; **`configure` / `configure.win`** no longer
mention GitHub, `install_github`, Remote* fields, or custom Rcpp version checks.

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
addressed in **0.9.73** as described above.

### Win-builder (2026-08-02)

**glmbayes 0.9.73** (this submission tarball) checked on win-builder:

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

* **R-universe (2026-08-02, glmbayes 0.9.73):** passes cleanly (0 errors,
  warnings, or NOTEs) on all built platforms except **wasm** (not supported;
  package requires native compiled C++/OpenCL stack).

---
_This file is listed in `.Rbuildignore` and is not included in the built source
tarball. When submitting, paste the content above into the “Optional comments”
field on the CRAN submission form at_
https://cran.r-project.org/submit.html
