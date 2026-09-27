# C++ → R function call inventory (glmbayes & glmbayesCore)

Canonical reference for which **R functions** are invoked from **`src/*.cpp`**.  
Regenerate checks with:

```text
Rscript data-raw/cpp_r_callback_inventory.R
Rscript data-raw/cpp_r_callback_inventory.R C:/Rpackages/glmbayesCore
```

The script matches `(Rcpp::)Function optionalName("symbol")`, `namespace_env("pkg")["symbol"]`,
and `pkg["symbol"]` / `glmbayes_ns["symbol"]` (line comments stripped). Core files using
`glmbayes_R::` accessors are reported but not whitelisted as literals.

**Out of scope for this document:** direct **libR / Rmath** use (`R::qnorm`, `R::rgamma`, `Rf_pgamma`, `R::runif`, …). Those are C API calls to Mathlib, not R function lookup.

**Mechanisms:**

| Mechanism | Description |
|-----------|-------------|
| `Rcpp::Function("name")` | Resolve `name` on the R search path |
| `namespace_env("pkg")["symbol"]` | Explicit package namespace |
| `glmbayes_R::r_*()` | Core: cached accessor → `pkg_env()["symbol"]` |
| `Rcpp::Function` parameters | Closures passed from R (e.g. `f2`, `f3`) |

---

## 1. Package-local R (engine)

Type 1–style symbols implemented in R in each package (`simulationpipeline.R`, `gamma_ct.R`, `fitter_functions.R`, …).

| R symbol | glmbayes `src/` | glmbayesCore `src/` | Notes |
|----------|-----------------|---------------------|--------|
| **EnvelopeOpt** | `EnvelopeSize.cpp`, `EnvelopeBuild.cpp`, `EnvelopeBuild_Ind_Normal_Gamma.cpp` — `Rcpp::Function("EnvelopeOpt")` | Same files — `glmbayes_R::r_envelope_opt()` → `pkg_env()["EnvelopeOpt"]` | Grid size optimization |
| **EnvelopeSort** | `EnvelopeBuild.cpp`, `EnvelopeBuild_Ind_Normal_Gamma.cpp` — `"EnvelopeSort"`; `EnvelopeOrchestrator.cpp` — `namespace_env("glmbayes")["EnvelopeSort"]` | `EnvelopeOrchestrator.cpp`, EnvelopeBuild* — `r_envelope_sort()` | Sort envelope grid |
| **glmbfamfunc** | `rNormalGammaReg.cpp` — `"glmbfamfunc"`; `rIndepNormalGammaReg.cpp` — `glmbayes_ns["glmbfamfunc"]` | `rNormalGammaReg.cpp`, `rIndepNormalGammaReg.cpp` — `r_glmbfamfunc()` | Supplies **f2** / **f3** list elements |
| **rNormal_reg.wfit** | `rNormalGammaReg.cpp` — `"rNormal_reg.wfit"` | `rNormalGammaReg.cpp` — `r_rNormal_reg_wfit()` | Normal–Gamma WLS |
| **rgamma_ct** | `rGammaGaussian.cpp`, `rGammaGamma.cpp` — `"rgamma_ct"` | Same — `r_rgamma_ct()` → `pkg_env()["rgamma_ct"]` | Truncated Gamma on precision scale; **only CT helper called as R from C++** |

glmbayes declares matching accessors in `src/R_interface.h` but most `.cpp` files still use string / global lookup; Core uses `package_ns.h` + `pkg_env()`.

**Not R callbacks:** truncated Normal / inverse-Gamma in accept–reject paths use C++ (`rnorm_ct.cpp`, `rinvgamma_ct_safe` in `rng_utils.cpp`).

---

## 2. Other R packages

| R symbol | Package | Files |
|----------|---------|--------|
| **get_opencl_core_count** | **opencltools** | glmbayes: `kernel_loader.cpp`; glmbayesCore: `OpenCL_helper.cpp` via `namespace_env("opencltools")` |

Legacy (Core `src/backup/` only, not active compile): `system.file` in `kernel_loader_fat_pre_opencltools.cpp`.

---

## 3. Base R / stats

Same usage patterns in both packages unless noted.

| R symbol | Typical use | `.cpp` files |
|----------|-------------|--------------|
| **optim** | Posterior mode (ING / GLM setup) | `rIndepNormalGammaReg.cpp`, `rNormalGLM.cpp` |
| **try** | Wrap **optim** | `rNormalGLM.cpp` |
| **gaussian** | Argument to **glmbfamfunc()** | `rNormalGammaReg.cpp`, `rIndepNormalGammaReg.cpp` |
| **lm.wfit** | **EnvelopeCentering** | `EnvelopeCentering.cpp` |
| **lm.fit** | **rNormalReg** | `rNormalReg.cpp` |
| **as.matrix**, **as.vector** | Coercion | `rNormalGLM.cpp`, `rNormalReg.cpp` |
| **expand.grid** | Envelope grids | `EnvelopeBuild.cpp`, `EnvelopeBuild_Ind_Normal_Gamma.cpp` |
| **qgamma**, **runif** | **rGammaGamma** accept–reject | `rGammaGamma.cpp` |
| **interactive**, **readline** | Long-run prompts | `EnvelopeEval.cpp`, `rNormalGLM.cpp`, `rIndepNormalGammaReg.cpp` |
| **format**, **Sys.time**, **as.numeric** | Timing (some commented) | `rIndepNormalGammaReg.cpp`, `EnvelopeDispersionBuild.cpp` (`run_ub2_pilot_block`) |

---

## 4. R closures passed from R (no name lookup in C++)

| Parameter | Origin | C++ use |
|-----------|--------|---------|
| **f2**, **f3** | `glmbfamfunc(family)` → `R/simfunction.R` → `.Call(..., f2, f3, ...)` | **optim** `fn` / `gr` in `rNormalGLM.cpp`, `rIndepNormalGammaReg.cpp`; signatures for `rNormalGLM_std*`. Parallel/std accept–reject often uses C++ `f2_*` in `famfuncs_*.cpp` instead of the R closure. |
| **ub2_parallel_fn** | Would be R | `run_ub2_pilot_block` in `EnvelopeDispersionBuild.cpp` calls it but **no caller** in either repo (reserved / dead). |

`export_wrappers.cpp` / `RcppExports.cpp` only marshal **f2** / **f3** SEXP; they do not invoke R.

---

## 5. glmbayes-only

| File | R calls |
|------|---------|
| `rNormalGLMBlocks.cpp` | Uses `Rcpp::Function` type only; no `Rcpp::Function("...")` lookups (blocks path uses C++ famfuncs). |

glmbayesCore has no active `rNormalGLMBlocks.cpp`.

---

## 6. glmbayes vs glmbayesCore (lookup style)

| Concern | glmbayes | glmbayesCore |
|---------|----------|--------------|
| Package-local symbols | `Rcpp::Function("Symbol")` or `namespace_env("glmbayes")` | `glmbayes_R::` + `pkg_env()` in `R_interface.h` |
| **glmbfamfunc** | Mixed | `r_glmbfamfunc()` |
| **rgamma_ct** | Global `"rgamma_ct"` | Namespace-bound **rgamma_ct** |
| OpenCL scaling | `kernel_loader.cpp` | `OpenCL_helper.cpp` |

The **set of R symbols** touched from C++ is the same; Core binds the package namespace explicitly.

---

## 7. Reference for later migration (no code changes implied here)

When phasing Type 1 symbols out of glmbayes exports, C++ (or Core) must still provide **EnvelopeOpt**, **EnvelopeSort**, **glmbfamfunc**, **rNormal_reg.wfit**, **rgamma_ct**, and **f2**/**f3** until the engine no longer needs them. **opencltools::get_opencl_core_count** remains a cross-package dependency. Removing **rgamma_ct** R sources requires either keeping R registration for C++ callbacks or a future C++ port (not part of this inventory).

When adding a new `Rcpp::Function("...")` literal in `src/`, update this file and the whitelist in `data-raw/cpp_r_callback_inventory.R`, then re-run the scan script.
