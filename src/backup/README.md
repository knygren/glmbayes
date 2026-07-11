# Archived OpenCL program loaders (not compiled)

Reference copies of superseded GPU program assembly code. These files are
**not** linked into the glmbayes shared library.

## `load_likelihood_subgradient_program_v1.cpp`

Pre-nmathopencl loader: prelude, shims, and selective `nmath/` were all read
from the **glmbayes** package (`inst/cl/`). Replaced by the production loader
in `src/kernel_loader.cpp`, which loads prelude and nmath from **nmathopencl**
via `program_preload_manifest.tsv` and `load_library_for_kernel_cross_package`.

Archived: 2026-06-21.
