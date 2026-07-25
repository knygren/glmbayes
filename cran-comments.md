# CRAN submission comments — glmbayes 0.9.71

## Summary

This is a patch release from glmbayes 0.9.7 (current CRAN release) to 0.9.71.

Main change:

* **`configure.win`:** Windows builds now pass `-DRCPP_PARALLEL_USE_TBB=1` and
  append `RcppParallel::RcppParallelLibs()` to **`PKG_LIBS`**, so linking matches
  updated CRAN/R-Universe **`RcppParallel`** TBB requirements on Windows.

No other functional or API changes.

OpenCL support remains **optional**. OpenCL-specific **testthat** tests call
`skip_on_cran()` (in addition to `skip_if_no_opencl()`).

## Test environments

* local Windows 10, `R CMD check --as-cran`: (to be run before submit)

* win-builder (CRAN): (to be confirmed on resubmission)

---
_This file is listed in `.Rbuildignore` and is not included in the built source
tarball. When submitting, paste the **Summary** section above into the “Optional
comments” field on the CRAN submission form at_
https://cran.r-project.org/submit.html
