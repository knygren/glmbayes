# Phase 2 engine snapshot (glmbayes → glmbayesCore backend)

Frozen copies of **glmbayes** engine sources before Phase 2 re-exports.
Not loaded by the installed package (`legacy_code/` is excluded from `R CMD build`).

- **Date:** 2026-09-29
- **Pin:** `glmbayesCore (>= 0.5.5)`
- **Active API:** re-exports in `R/reexports-core-engine.R` and `R/reexports.R` (`Prior_Setup`)

## Contents

- `R/` — former `glmb.R`, `lmb.R`, `pfamily.R`, sim/envelope pipeline, C++ R bindings, etc.
- `src/` — former glmbayes native sampler (OpenCL + Rcpp)
- `configure`, `configure.win` — former package configure scripts

## Recovery

Copy from here or use git history. Phase 3 may restore formula-layer code into `R/glmb.R` / `R/lmb.R`
while keeping sampling on Core.
