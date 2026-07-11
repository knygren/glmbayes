# Dump production OpenCL program strings via opencltools R loaders.
suppressPackageStartupMessages({
  if (!requireNamespace("opencltools", quietly = TRUE)) {
    stop("opencltools required")
  }
  if (!requireNamespace("nmathopencl", quietly = TRUE)) {
    stop("nmathopencl required")
  }
  if (!requireNamespace("glmbayes", quietly = TRUE)) {
    stop("glmbayes required")
  }
})

app_pkg <- "glmbayes"
nmath_pkg <- "nmathopencl"

assemble_program <- function(kernel_rel) {
  preload <- opencltools::load_program_preload(source_package = nmath_pkg)
  nmath_src <- opencltools::load_library_for_kernel_cross_package(
    kernel_relative_path = kernel_rel,
    kernel_package = app_pkg,
    library_subdir = "nmath",
    library_package = nmath_pkg,
    depends_tag = "all_depends_nmath"
  )
  ksrc <- opencltools::load_kernel_source(kernel_rel, app_pkg)
  paste(preload, nmath_src, ksrc, sep = "\n")
}

out_dir <- file.path("data-raw", "compare_v1_v2_cpp_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

kernels <- list(
  poisson = "src/f2_f3_poisson.cl",
  gaussian = "src/f2_f3_gaussian.cl",
  gamma = "src/f2_f3_gamma.cl",
  binomial_logit = "src/f2_f3_binomial_logit.cl",
  binomial_probit = "src/f2_f3_binomial_probit.cl",
  binomial_cloglog = "src/f2_f3_binomial_cloglog.cl"
)

for (name in names(kernels)) {
  program <- assemble_program(kernels[[name]])
  path <- file.path(out_dir, paste0(name, ".cl"))
  writeLines(program, path, useBytes = TRUE)
  message(sprintf("%s: nchar=%d -> %s", name, nchar(program, type = "bytes"), path))
}
