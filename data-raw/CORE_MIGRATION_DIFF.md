# glmbayes → glmbayesCore backend migration inventory

Generated for the outside-in relocation (see `Package_Marketing/README.md` and
`glmbayesCore/README.md` § Function overview). Pin: **glmbayesCore (>= 0.5.4)**.

## Summary counts (R + src + inst/cl; excludes build artifacts in diff list)

| Class | Count | Action |
|-------|------:|--------|
| Byte-identical paths | 73 | Phase out glmbayes copy when wiring complete |
| Differing paths | 90 | Mostly namespace/dynlib/build noise; source twins → **Keep Core** |
| glmbayes-only | 28 | Keep in glmbayes (formula/S3/vignettes/blocks) |
| Core-only | 9 | Core source of truth (`multi_*`, `residuals.rglmb.R`, …) |

## glmbayes-only (keep)

- **User / S3 / vignettes:** `R/glmb.R`, `R/lmb.R` (formula + mlmb helpers), `R/summary.glmb.R`, `R/predict.glmb.R`, `R/residuals.glmb.R` (glmb/lmb methods), `R/prior_simfunction.R`, `R/directional_tail.R`, `R/dic_info.R` (until DIC fully delegated), insight/bayestestR methods, etc.
- **Blocks sampler (not in Core):** `src/rNormalGLMBlocks.cpp`, `inst/DESIGN_RGLM_BLOCKS.md` — **removed with local engine**; port to Core in Stage 7 if block-Gibbs benchmarks need it (no formula-path caller today).

## Core-only (consume via Imports / S3)

- `R/multi_*.R`, `R/summary.mrglmb.R`, `R/residuals.rglmb.R` (S3 for rglmb/rlmb)
- `src/package_ns.h` (Core dynlib identity)

## Differing phase-out files — decision

| Area | Decision |
|------|----------|
| `R/prior.R`, `R/pfamily.R`, `R/simfunction.R`, `R/simulationpipeline.R` | **Keep Core** — glmbayes forwards/re-exports |
| `R/rglmb.R`, `R/rlmb.R` | **Keep Core** — thin wrappers keep glmbayes documentation |
| `R/rcpp_wrappers.R`, `R/RcppExports.R` | **Delete** with `src/` (Phase 6) |
| `src/*.cpp`, `configure*`, `inst/cl/` | **Delete** from glmbayes (Phase 6); OpenCL via Core + nmathopencl + opencltools |
| `R/gpu_diagnostics.R` | **Keep Core** — `has_opencl()` aliases `glmbayesCore_has_opencl()` |
| `R/get_opencl_core_count.R` | **opencltools** (same as Core compile-time policy) |
| Shared datasets (`R/data-*.R`) | **Keep glmbayes** copies (LazyData for vignettes) |

## Rule of thumb applied

- Algorithm / sampler / envelope code → **glmbayesCore**
- `glmb()` / `lmb()` + S3 + docs → **glmbayes**
- GPU prelude/nmath → **nmathopencl**; host loaders → **opencltools** (Dependencies of Core, not duplicated in glmbayes after Phase 6)

## OpenCL / CRAN binaries (delegation note)

After wiring `glmb()` / `lmb()` to **`glmbayesCore::rglmb()` / `rlmb()`**, OpenCL compile-time support comes from **glmbayesCore** (and **opencltools**, **nmathopencl**). **CRAN/R-Universe binaries are CPU-only**; `use_opencl = TRUE` requires **source installs** of that stack. Parity and CI should default to **`use_opencl = FALSE`**; GPU smoke tests are optional on source-built environments only.

## Three-type map (see plan § R function taxonomy)

- **Type (1) — Core-only:** envelope, simfuncs, `*_ct`, `simfunction`, … — **no** `glmbayes::` export; migrate **first**.
- **Type (2) — Core + glmbayes re-export:** `Prior_Setup`, `pfamily`, `d*`, `rglmb`, `rlmb`, `diagnose_glmbayes`, …
- **Type (3) — glmbayes-only:** `glmb`, `lmb`, S3 methods, `prior_simfunction`, insight/bayestestR, …

## Implementation order

1. **Type (1):** remove glmbayes duplicates/exports; examples use **`glmbayesCore::`**; **`glmb()`** may stay on local engine temporarily.
2. **Type (2):** re-export from Core.
3. **`glmb()` / `lmb()`** → Core last; drop local **`src/`** when no `_glmbayes_*` path remains.
