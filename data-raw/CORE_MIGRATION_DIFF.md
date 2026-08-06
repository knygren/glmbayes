# Stage 0 — glmbayes vs glmbayesCore file inventory

**Date:** 2026-08-06  
**glmbayes:** 0.9.76 (development)  
**glmbayesCore (baseline compared):** 0.5.3  
**Stage 1 pin (after Stage 0 Core updates):** `Imports: glmbayesCore (>= 0.5.4)`

Inventory of phase-out candidates under `R/`, `src/`, `inst/cl/`, plus top-level
`configure` / `configure.win` / `DESCRIPTION` / `NAMESPACE`. Build artifacts
(`*.o`, `*.dll`) excluded. Supporting path lists: `data-raw/_stage0_lists/`.

**Counts:** 131 matching paths → **70 identical**, **61 differ**; **26 glmbayes-only**; **12 Core-only**.

---

## Decision legend

| Decision | Meaning |
|----------|---------|
| **same** | Byte-identical; no action |
| **Keep Core** | Consume Core as-is; glmbayes delta is packaging, docs, or older style |
| **Update Core** | Port glmbayes improvement into Core before relying on that path |
| **keep-glmbayes** | Not a phase-out candidate (formula / S3 / ecosystem) |
| **Core-only** | Stays in Core; no glmbayes counterpart |

Expected packaging noise (do **not** merge): package name, `_glmbayes_*` vs
`_glmbayesCore_*`, `GLMBAYES_R_NS` / `package_ns.h`, Rdpack `{glmbayes}` vs
`{glmbayesCore}` cites.

---

## Update Core (required before later Stages)

| File | Why | Before Stage |
|------|-----|--------------|
| `configure` | glmbayes has **non-PoCL GPU** OpenCL probe (CRAN PoCL NOTE fix); Core still enables OpenCL on any platform. Also remove Core `tools/rcpp_include.R` probing (CRAN policy; already removed from glmbayes). Keep Core `-include glmbayes_getRegisteredNamespace.h`. | **0 / 1** (landed in Core 0.5.4 as part of Stage 0) |
| `configure.win` | glmbayes removed `rcpp_include` / GitHub-style Rcpp probing; Core still has it. Keep TBB flags + Core namespace shim include. | **0 / 1** (with configure) |
| `R/pfamily.R` (+ `R/prior_simfunction.R`) | glmbayes embeds **`pfun`** for bayestestR `simulate_prior`; Core constructors lack `pfun`. | **3** (before re-exporting pfamily) |
| `tools/rcpp_include.R`, `tools/patch_rcpp_function_h.R` | Delete from Core once configure no longer calls them (policy cleanup). | **0 / 1** |

---

## Keep Core (Core is source of truth / ahead)

| File | Notes |
|------|-------|
| `R/gpu_diagnostics.R` | Core: object + `print` (CRAN-clean). glmbayes: ungated `cat()`. |
| `R/prior.R` | Core: `message()` for status. glmbayes: `print()`. |
| `R/envelopeorchestrator.R` | Core: `message()`. glmbayes: `cat()`. Docs/cites only otherwise. |
| `R/rglmb.R`, `R/rlmb.R` | Mostly ns/docs; Core namespace-safe patterns. |
| `R/simfunction.R`, `R/simulationpipeline.R`, `R/summary.rglmb.R`, `R/formula.summary.rglmb.R` | Large diffs; Core multi-response / cleanup ahead. No glmbayes-only algorithm flagged for port except via pfamily `pfun` (above). |
| `R/rcpp_wrappers.R`, `R/RcppExports.R` | Dynlib symbol names; Core correct for Core. glmbayes extras = Blocks export only. |
| `R/zzz.R`, `R/globals.R`, CT helpers (`normal_ct`, `gamma_ct`, `invgamma_ct`) | Package-name / trivial. |
| `R/data-*.R` | Example path / cite packaging; datasets not migration-critical. |
| `src/R_interface.h` | **Keep Core** — registered-namespace callbacks (`package_ns`). |
| `src/Envelopefuncs.h`, `src/EnvelopeDispersionBuild.cpp` | **Keep Core** — `check_disp_bounds_or_stop`, ub2 diagnostics. |
| `src/rIndepNormalGammaReg.cpp`, related Envelope/sim `.cpp` | **Keep Core** — Core has more guard/diagnostic logic. |
| `src/rNormalGLM.cpp` | **Keep Core** — Core suppresses noisy `kappa_H` NOTE (mixed-model paths). |
| `src/progress_utils.*`, `src/rng_utils.cpp`, `src/famfuncs_*.cpp`, `src/kernel_*.cpp`, `src/opencl*.h`, `src/OpenCL_helper.cpp` | Packaging / comment / Core OpenCL layout. |
| `src/export_wrappers.cpp`, `src/simfuncs.h` | Core without Blocks; Blocks stays glmbayes-only until ported. |
| `src/Makevars`, `src/Makevars.win` | Generated / local; ignore for merge. |
| `inst/cl/README.md`, `inst/cl/cpp/*` | Package name / minor; entry `.cl` kernels mostly **same**. |
| `DESCRIPTION`, `NAMESPACE` | Package identity — not merge targets. |

---

## Identical (`same`) — no action

70 paths. Highlights: most `inst/cl/nmath|R_*|src/f2_f3_*.cl`, plus
`R/compute_gaussian_prior.R`, `R/fitter_functions.R`, `R/internal_rcppparallel.R`,
`R/summary.rgamma_reg.R`, and several `src/` files
(`EnvelopeEval.cpp`, `EnvelopeSort.cpp`, `famfuncs.h`, `famfuncs_poisson.cpp`,
`kernel_runners.cpp`, `rNormalReg.cpp`, `Set_Grid.cpp`, `Set_LogP.cpp`,
`configure_OpenCL.cpp`, `opencl_detect.cpp`, CT cpp, etc.).

Full list: [`_stage0_lists/same.txt`](./_stage0_lists/same.txt).

---

## glmbayes-only — keep in glmbayes (or defer)

| File | Fate |
|------|------|
| `R/glmb.R`, `R/lmb.R` | **keep-glmbayes** — formula API (Stages 1+) |
| `R/*glmb*.R` S3, insight/bayestestR methods, `reexports.R`, `glmbayes-package.R` | **keep-glmbayes** |
| `R/directional_tail.R`, `R/extractDIC.R`, `R/prior_simfunction.R`, `R/get_opencl_core_count.R` | **keep-glmbayes** until optionally moved; `prior_simfunction.R` / `pfun` → Core before Stage 3 |
| `src/rNormalGLMBlocks.cpp` (+ wrapper exports) | **keep-glmbayes** for now; **not called** from non-wrapper R code. Port to Core or drop before Stage 6 if unused. |
| `src/backup/load_likelihood_subgradient_program_v1.cpp` | backup; ignore |

---

## Core-only — stay in Core

| File | Notes |
|------|-------|
| `R/multi_*.R`, `R/summary.mrglmb.R`, `R/multi_prior_setup.R` | Multi-response API |
| `R/dic_info.R`, `R/ing_prior_guard.R`, `R/residuals.rglmb.R` | Core helpers |
| `R/glmbayesCore-package.R` | Package docs |
| `src/package_ns.h`, `src/glmbayes_getRegisteredNamespace.*` | Namespace-safe C++→R |
| `src/backup/kernel_loader_fat_pre_opencltools.cpp` | backup |

---

## Parity checklist (for Stages 1+)

Re-run after each Stage that changes call paths:

1. Gaussian LM — conjugate (`lmb` / `rlmb` + `dNormal_Gamma` / ING as applicable)
2. Binomial logit — envelope (`glmb` / `rglmb` + `dNormal`)
3. Poisson — envelope
4. Gamma — as supported
5. `Prior_Setup()` defaults on a small `glm`/`lm`
6. OpenCL on/off when a non-PoCL GPU is present (`use_opencl = TRUE` smoke)
7. `diagnose_glmbayes()` returns a printable object (Core behavior)

---

## Stage 0 actions taken

1. This inventory and decisions recorded.
2. **glmbayesCore 0.5.4:** port PoCL GPU probe; remove `rcpp_include` configure path (policy); keep Core namespace shim in Makevars flags; drop unused `tools/rcpp_include.R` / patch helper once unused.
3. Stage 1 will use `glmbayesCore (>= 0.5.4)`.
