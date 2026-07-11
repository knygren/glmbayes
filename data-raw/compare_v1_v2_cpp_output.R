# Dump production OpenCL program strings (R-level assembly via opencltools).
# Mirrors load_likelihood_subgradient_program() without a package C++ export.
suppressPackageStartupMessages({
  if (!requireNamespace("opencltools", quietly = TRUE)) {
    stop("opencltools required")
  }
  if (!requireNamespace("nmathopencl", quietly = TRUE)) {
    stop("nmathopencl required")
  }
})

app_pkg <- "glmbayes"
nmath_pkg <- "nmathopencl"

assemble_preload <- function(pkg = nmath_pkg) {
  manifest_path <- system.file("cl", "program_preload_manifest.tsv", package = pkg)
  if (!nzchar(manifest_path)) {
    stop("program_preload_manifest.tsv not found in ", pkg)
  }
  manifest <- read.delim(manifest_path, stringsAsFactors = FALSE)
  pieces <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, ]
    pieces[[i]] <- if (row$kind == "file") {
      opencltools::load_kernel_source(row$rel_path, pkg)
    } else if (row$kind == "library") {
      opencltools::load_kernel_library(row$rel_path, pkg)
    } else {
      stop("Unknown kind: ", row$kind)
    }
  }
  paste(unlist(pieces), collapse = "\n")
}

assemble_program <- function(kernel_rel, app = app_pkg, nmath = nmath_pkg) {
  preload <- assemble_preload(nmath)
  kernel_path <- system.file("cl", kernel_rel, package = app)
  nmath_dir <- system.file("cl/nmath", package = nmath)
  nmath_src <- opencltools::load_library_for_kernel(
    kernel_path, nmath_dir, depends_tag = "all_depends_nmath"
  )
  ksrc <- opencltools::load_kernel_source(kernel_rel, app)
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
