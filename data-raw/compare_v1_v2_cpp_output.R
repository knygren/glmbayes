# Dump C++ assembled OpenCL program strings (production loader).
# Former v1 vs v2 comparison removed; v1 archived in src/backup/.
suppressPackageStartupMessages({
  pkgload::load_all(quiet = TRUE)
})

if (!has_opencl()) {
  stop("OpenCL not enabled in this glmbayes build.")
}

out_dir <- file.path("data-raw", "compare_v1_v2_cpp_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

kernels <- list(
  poisson = list(family = "poisson", link = "logit"),
  gaussian = list(family = "gaussian", link = "logit"),
  gamma = list(family = "Gamma", link = "logit"),
  binomial_logit = list(family = "binomial", link = "logit"),
  binomial_probit = list(family = "binomial", link = "probit"),
  binomial_cloglog = list(family = "binomial", link = "cloglog")
)

for (name in names(kernels)) {
  spec <- kernels[[name]]
  out <- glmbayes:::debug_likelihood_program_cpp_export(spec$family, spec$link)
  path <- file.path(out_dir, paste0(name, ".cl"))
  writeLines(out$program, path, useBytes = TRUE)
  message(sprintf("%s: nchar=%d -> %s", name, out$nchar, path))
}
