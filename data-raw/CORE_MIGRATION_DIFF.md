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

---

## Checklist — `glmb()` / `lmb()` via Core backend (three phases)

**Goal:** One sampler stack and one OpenCL dynlib (**`glmbayesCore`**) before splitting
formula-layer code back into **glmbayes**.

**Golden rule for Phase 3:** Move **formula / object shape / glmb-only S3** back into
glmbayes one function at a time. **Never** move **`rglmb` / `rlmb` / simfuncs / envelope /
`src/`** back into glmbayes. Engine calls stay **`glmbayesCore::…`** (or re-exports
that assign from Core).

**Pin during this work:** `Imports: glmbayesCore (>= 0.5.5)` (raise when Core gains staged
fitters).

### Phase 1 — Stage formula fitters in glmbayesCore; validate there

Use internal export names in Core if you want to avoid committing to permanent
`glmbayesCore::glmb` / **`lmb`** as permanent Core API — Phase 1 uses the **same
names** temporarily in Core, then **glmbayes** re-exports in Phase 2 and owns formula
code again in Phase 3.

- [x] **Core: port post-sample assembly** (inlined in **`R/glmb.R`** and **`.uni_lmb`** in
  **`R/lmb.R`**: Prior, **`DIC_Info`**, dispersion branches, `outlist`, classes).
- [x] **Core: port `glmb()`** in **`R/glmb.R`** (`glm.fit` → **`rglmb()`**).
- [x] **Core: port `lmb()`** in **`R/lmb.R`** (single + **`mlmb`** via **`.uni_lmb`**).
- [x] **`multi_prior_setup()`** — already in Core; not duplicated in **`lmb.R`** extract.
- [x] **Core: COPYRIGHTS** — **skipped** for Phase 1 temporary staging (notices stay on **glmbayes**).
- [x] **Core: do not** add insight / bayestestR Imports for staged fitters.
- [x] **Core tests:** **`tests/testthat/test-glmb.R`**, **`test-lmb.R`** (CPU); **`test-opencl-glmb.R`**, **`test-opencl-lmb.R`**
  when **`glmbayesCore_has_opencl()`** is **`TRUE`** (smoke + structure checks; skipped on CRAN).
- [ ] **Core tests (optional CI job):** one smoke fit with **`use_opencl = TRUE`** on
  source-built Core + opencltools + nmathopencl.
- [x] **Assert simfun namespace:** tests check **`environment(pfamily$simfun)`** is Core.

### Phase 2 — Re-export from glmbayes; disable local copies; single backend

- [ ] **glmbayes: Type (2) re-export** (before or with fitters): **`pfamily`**, **`d*`**, prior
  **`r*_prior`**, **`rglmb`**, **`rlmb`**, **`simfunction`**, envelope exports, **`Prior_Setup`**
  (done), **`diagnose_glmbayes`**, etc. — `@inherit` + assignment pattern in **`R/reexports.R`**
  (or dedicated doc stubs where S3-only).
- [ ] **glmbayes: re-export staged fitters**, e.g. `glmb <- glmbayesCore::glmb` and
  `lmb <- glmbayesCore::lmb` with roxygen `@inherit` (same pattern as **`Prior_Setup`**).
- [ ] **glmbayes: comment out** (archive) local implementations in **`R/glmb.R`**, **`R/lmb.R`**
  (and duplicate **`multi_prior_setup`** in **`lmb.R`** if Core is canonical) so roxygen does
  not register duplicate exports.
- [ ] **glmbayes: remove duplicate S3** from NAMESPACE where Core registers via Imports
  (**`pfamily.default`**, **`simfunction.default`**, summary/residuals for rglmb, etc.).
- [ ] **glmbayes: drop `src/`** when no R code path calls **`useDynLib(glmbayes)`** symbols
  (remove **`R/RcppExports.R`**, **`R/rcpp_wrappers.R`**, configure OpenCL in glmbayes, etc.).
- [ ] **glmbayes tests:** full **`devtools::test()`** / examples with **`library(glmbayes)`**
  only (no attached Core).
- [ ] **Confirm OpenCL path:** **`diagnose_glmbayes()`** / **`has_opencl()`** reflect Core +
  opencltools; main fits do not load glmbayes duplicate C++.
- [ ] **Re-run parity checklist** (§ below) on glmbayes after re-exports.

### Phase 3 — Move formula layer back to glmbayes one function at a time

Each item: restore **glmbayes** source + docs; **keep** sampling and priors on Core.

- [x] **`.uni_lmb()` / single-response `lmb()`** — local **`model.frame` / `lm.fit`** in
  glmbayes; sampling line **`glmbayesCore::rlmb(...)`**; assembly local; **`mlmb`**
  helpers in **`R/lmb.R`**.
- [x] **`glmb()`** — local glm preamble; **`glmbayesCore::rglmb(...)`**; assembly in **`R/glmb.R`**.
- [ ] **Multi-response `lmb()` / `mlmb`** — local orchestration; per-column
  **`glmbayesCore::rlmb`** (or one **`multi_rlmb`** when aligned with **`.mlmb_assemble`**).
- [ ] **`multi_prior_setup`** — glmbayes docs + **`glmbayesCore::multi_prior_setup`** or
  thin wrapper only.
- [x] **Remove temporary Core exports** — **`glmb()`** / **`lmb()`** removed from Core 0.5.5 (keep
  **`rglmb` / `rlmb`** in Core permanently).
- [ ] After **each** Phase 3 function: parity tests + note in NEWS; bump Core version pin if
  assembly API changed.

### Parity gates (run after Phase 1, 2, and each Phase 3 step)

Same as § **Parity checklist** below; minimum:

- [ ] Gaussian — conjugate **`lmb`** / **`rlmb`** + **`dNormal_Gamma`** (or ING as used in tests).
- [ ] Binomial logit — **`glmb`** / **`rglmb`** + **`dNormal`** envelope path.
- [ ] Poisson — envelope path; Poisson + **`dGamma`** coercion if covered in tests.
- [ ] **`Prior_Setup()`** → **`glmb` / `lmb`** on small model.
- [ ] **`simulate_prior.glmb`** still finds **`pfamily$pfun`** (Core constructors with **`pfun`**).
- [ ] **`summary.glmb`**, **`directional_tail`**, **`extractDIC`** on object from new path.

### Anti-patterns (do not do in Phase 3)

- [ ] Do **not** copy **`rglmb` / `rlmb` / `R/simfunction.R` / `simulationpipeline.R`** back into
  glmbayes.
- [ ] Do **not** re-enable glmbayes **`src/`** for the main formula path.
- [ ] Do **not** build **`pfamily`** with local **`dNormal()`** while calling Core **`rglmb`**
  unless **`dNormal`** is a re-export from Core (otherwise **`simfun`** binds to glmbayes simfuncs).

---

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
