# CRAN submission comments — glmbayes 0.9.73

## Summary

This is a patch release from glmbayes 0.9.72 (current CRAN release) to 0.9.73.

0.9.72 was flagged for a policy violation as a configure warning recommended github install
of Rcpp on certain dev versions of R (present because of issue with Rcpp during the transition to R 4.6).
This issues was described on the Rcpp github page. https://github.com/RcppCore/Rcpp/issues/1473

The changes in this submission is aimed at removing this policy violatin as per the below.

Main changes:

* **`configure.win`:** Windows builds now pass `-DRCPP_PARALLEL_USE_TBB=1` and
  append `RcppParallel::RcppParallelLibs()` to **`PKG_LIBS`**, so linking matches
  updated CRAN/R-Universe **`RcppParallel`** TBB requirements on Windows.
  No longer invokes **`tools/rcpp_include.R`** (Rcpp via standard **`LinkingTo`**).

* **`tools/rcpp_include.R`:** Removed configure warnings that recommended
  installing **Rcpp** from GitHub for certain .

No other functional or API changes.

OpenCL support remains **optional**. OpenCL-specific **testthat** tests call
`skip_on_cran()` (in addition to `skip_if_no_opencl()`).

## Test environments

* local Windows 10, `R CMD check --as-cran`: 

* win-builder (CRAN): clean on r-release and r-devel (0.9.73)

---
_This file is listed in `.Rbuildignore` and is not included in the built source
tarball. When submitting, paste the **Summary** section above into the “Optional
comments” field on the CRAN submission form at_
https://cran.r-project.org/submit.html
