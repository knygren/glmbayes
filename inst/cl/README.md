# OpenCL sources (`inst/cl`)

This directory is installed as **`cl/`** in the built package. Kernel paths passed to `system.file("cl", …, package = "glmbayes")` resolve here.

It supplies the **OpenCL C** sources used for GPU evaluation of **standard-form `f2` / `f3`** (negative log-posterior fragment and its gradient in \(\beta\)) during **likelihood-subgradient envelope** sampling (Nygren & Nygren, 2006). The CPU statistical path is unchanged; OpenCL is optional (`USE_OPENCL`).

---

## Layout

| Path | Role |
|------|------|
| **`OPENCL.cl`** | Prelude: extensions (e.g. `cl_khr_fp64`), feature macros (`HAVE_EXPM1`, `HAVE_WORKING_ISFINITE`, …), `ML_NAN` / infinities, `INLINE`, `R_UNUSED`. Stitched **first** into every assembled program. |
| **`libR_shims/`** | Minimal stand-ins for symbols expected by ported code. |
| **`R_ext_types/`**, **`R_shims/`**, **`R_ext_runtime/`**, **`R_ext_internals/`**, **`System/`** | Headers/runtime fragments adapted for OpenCL C so R/nmath-style sources compile on-device. Loaded as whole libraries (dependency-sorted `.cl` files per directory). |
| **`nmath/`** | Ported **nmath-related** routines (densities, helpers, `Rmath.cl`, etc.). **Not** always loaded in full: see below. |
| **`nmath/kernel_dependency_index.tsv`** | Stem load order for selective inclusion of `nmath/*.cl` files. |
| **`src/`** | Entry kernels `f2_f3_*.cl` (one **\_\_kernel** per family/link). Each file lists **`@all_depends_nmath`** so only the needed `nmath` stems are concatenated. |

Legacy layouts (`rmath/`, `dpq/` as separate trees for old concatenation) are **not** used by the current default GPU path.

---

## Building one executable OpenCL program

At runtime, **`glmbayes::opencl::load_likelihood_subgradient_program(family, link, app_package, nmath_package)`** (implemented in `src/kernel_loader.cpp`) returns a **single character string**. Prelude, shims, and selective **`nmath/`** are read from **`nmathopencl`**; the entry kernel from **`glmbayes`** (`inst/cl/src/f2_f3_*.cl`).

Assembly order:

1. **`OPENCL.cl`** (from **nmathopencl**)
2. **`libR_shims`**, **`R_ext_types`**, **`R_shims`**, **`R_ext_runtime`**, **`R_ext_internals`**, **`System`** (from **nmathopencl**)
3. **`nmath`** subset — stems from `@all_depends_nmath` on the chosen entry kernel; files and TSV from **nmathopencl**
4. **Entry kernel** — e.g. `src/f2_f3_binomial_logit.cl` (from **glmbayes**)

Family/link strings match R’s GLM conventions (`"binomial"` / `"logit"`, `"poisson"`, `"Gamma"`, `"gaussian"`, …) and map internally to the correct `src/` path.

That string is passed to **`clCreateProgramWithSource`** (see `glmbayes::opencl::f2_f3_kernel_runner` in `src/kernel_runners.cpp`), built for the device, and the matching **`\_\_kernel`** (`f2_f3_binomial_logit`, etc.) is launched from **`f2_f3_opencl`** (`src/kernel_wrappers.cpp`).

---

## Relation to R helpers

- **`opencltools::load_kernel_source()`** / **`opencltools::load_kernel_library()`** (with `package = "glmbayes"` for entry kernels, or `"nmathopencl"` for prelude/nmath) are **generic** exploration helpers.
- **`load_likelihood_subgradient_program()`** is **application-specific**: the exact C++ recipe for envelope GPU evaluation (delegates loading to **opencltools** via C API). The pre-nmathopencl loader is archived under **`src/backup/`**.

---

## Editing conventions

- **`nmath/`** `.cl` files and **`src/f2_f3_*.cl`** use comment tags (`@provides`, `@depends`, `@all_depends_nmath`, …) consumed by the loader. If you add or rename stems, update **`kernel_dependency_index.tsv`** so transitive ordering stays consistent.
- Keep **`OPENCL.cl`** compatible with both the prelude expectations of **`nmath/`** ports and the entry kernels (double precision, etc.).

---

## References

- Nygren, K. N., & Nygren, L. M. (2006). Likelihood subgradient densities. *Journal of the American Statistical Association*, 101(475), 1144–1156. https://doi.org/10.1198/016214506000000357  
- Package vignettes: OpenCL chapter and appendix on kernel assembly (`Chapter-16`, `Chapter-A10`).
