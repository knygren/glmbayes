# Compare OpenCL program assembly: v1 (glmbayes-local) vs v2 (nmathopencl preload/nmath)
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

cat("=== Package paths ===\n")
cat("glmbayes OPENCL.cl:", system.file("cl", "OPENCL.cl", package = app_pkg), "\n")
cat("nmathopencl manifest:", system.file("cl", "program_preload_manifest.tsv", package = nmath_pkg), "\n")
cat("glmbayes nmath dir:", system.file("cl/nmath", package = app_pkg), "\n")
cat("nmathopencl nmath dir:", system.file("cl/nmath", package = nmath_pkg), "\n")
cat("nmathopencl version:", as.character(utils::packageVersion(nmath_pkg)), "\n\n")

assemble_preload_v1 <- function(pkg = app_pkg) {
  pieces <- c(
    opencltools::load_kernel_source("OPENCL.cl", pkg),
    opencltools::load_kernel_library("libR_shims", pkg),
    opencltools::load_kernel_library("R_ext_types", pkg),
    opencltools::load_kernel_library("R_shims", pkg),
    opencltools::load_kernel_library("R_ext_runtime", pkg),
    opencltools::load_kernel_library("R_ext_internals", pkg),
    opencltools::load_kernel_library("System", pkg)
  )
  paste(pieces, collapse = "\n")
}

assemble_preload_v2 <- function(pkg = nmath_pkg) {
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

assemble_nmath_v1 <- function(kernel_rel, pkg = app_pkg) {
  kernel_path <- system.file("cl", kernel_rel, package = pkg)
  nmath_dir <- system.file("cl/nmath", package = pkg)
  opencltools::load_library_for_kernel(
    kernel_path, nmath_dir, depends_tag = "all_depends_nmath"
  )
}

assemble_nmath_v2 <- function(kernel_rel, app = app_pkg, nmath = nmath_pkg) {
  kernel_path <- system.file("cl", kernel_rel, package = app)
  nmath_dir <- system.file("cl/nmath", package = nmath)
  opencltools::load_library_for_kernel(
    kernel_path, nmath_dir, depends_tag = "all_depends_nmath"
  )
}

assemble_program_v1 <- function(kernel_rel) {
  paste(
    assemble_preload_v1(),
    assemble_nmath_v1(kernel_rel),
    opencltools::load_kernel_source(kernel_rel, app_pkg),
    sep = "\n"
  )
}

assemble_program_v2 <- function(kernel_rel) {
  paste(
    assemble_preload_v2(),
    assemble_nmath_v2(kernel_rel),
    opencltools::load_kernel_source(kernel_rel, app_pkg),
    sep = "\n"
  )
}

md5_str <- function(x) {
  tmp <- tempfile(fileext = ".cl")
  on.exit(unlink(tmp))
  writeLines(x, tmp, useBytes = TRUE)
  unname(tools::md5sum(tmp))
}

first_diff_info <- function(a, b) {
  if (identical(a, b)) {
    return(list(same = TRUE))
  }
  len_a <- nchar(a, type = "bytes")
  len_b <- nchar(b, type = "bytes")
  min_len <- min(len_a, len_b)
  # compare by lines for readability
  la <- strsplit(a, "\n", fixed = TRUE)[[1]]
  lb <- strsplit(b, "\n", fixed = TRUE)[[1]]
  n <- min(length(la), length(lb))
  line_diff <- which(la[seq_len(n)] != lb[seq_len(n)])[1]
  list(
    same = FALSE,
    nchar_v1 = len_a,
    nchar_v2 = len_b,
    nchar_delta = len_b - len_a,
    nlines_v1 = length(la),
    nlines_v2 = length(lb),
    first_line_diff = line_diff,
    v1_line = if (!is.na(line_diff)) la[line_diff] else NA_character_,
    v2_line = if (!is.na(line_diff)) lb[line_diff] else NA_character_
  )
}

report_section <- function(label, v1, v2) {
  cat("\n--- ", label, " ---\n", sep = "")
  cat("v1 nchar:", nchar(v1, type = "bytes"), " md5:", md5_str(v1), "\n")
  cat("v2 nchar:", nchar(v2, type = "bytes"), " md5:", md5_str(v2), "\n")
  info <- first_diff_info(v1, v2)
  if (info$same) {
    cat("IDENTICAL\n")
  } else {
    cat("DIFFER\n")
    cat("  nchar delta (v2 - v1):", info$nchar_delta, "\n")
    cat("  nlines v1/v2:", info$nlines_v1, "/", info$nlines_v2, "\n")
    if (!is.na(info$first_line_diff)) {
      cat("  first differing line:", info$first_line_diff, "\n")
      cat("  v1:", substr(info$v1_line, 1, 120), "\n")
      cat("  v2:", substr(info$v2_line, 1, 120), "\n")
    } else if (info$nlines_v1 != info$nlines_v2) {
      cat("  lines differ in count only (content equal up to min length)\n")
    }
  }
}

kernels <- list(
  poisson = "src/f2_f3_poisson.cl",
  gaussian = "src/f2_f3_gaussian.cl",
  gamma = "src/f2_f3_gamma.cl",
  binomial_logit = "src/f2_f3_binomial_logit.cl",
  binomial_probit = "src/f2_f3_binomial_probit.cl",
  binomial_cloglog = "src/f2_f3_binomial_cloglog.cl"
)

cat("=== Preload comparison (all kernels share this) ===\n")
preload_v1 <- assemble_preload_v1()
preload_v2 <- assemble_preload_v2()
report_section("preload", preload_v1, preload_v2)

cat("\n=== Per-kernel comparison ===\n")
results <- list()
for (name in names(kernels)) {
  kernel_rel <- kernels[[name]]
  cat("\n========== ", name, " (", kernel_rel, ") ==========\n", sep = "")

  nmath_v1 <- assemble_nmath_v1(kernel_rel)
  nmath_v2 <- assemble_nmath_v2(kernel_rel)
  report_section("nmath", nmath_v1, nmath_v2)

  ksrc <- opencltools::load_kernel_source(kernel_rel, app_pkg)
  cat("\n--- kernel (same in v1/v2) ---\n")
  cat("nchar:", nchar(ksrc, type = "bytes"), " md5:", md5_str(ksrc), "\n")

  prog_v1 <- assemble_program_v1(kernel_rel)
  prog_v2 <- assemble_program_v2(kernel_rel)
  report_section("full program", prog_v1, prog_v2)

  results[[name]] <- list(
    preload_same = identical(preload_v1, preload_v2),
    nmath_same = identical(nmath_v1, nmath_v2),
    program_same = identical(prog_v1, prog_v2),
    program_nchar_v1 = nchar(prog_v1, type = "bytes"),
    program_nchar_v2 = nchar(prog_v2, type = "bytes")
  )
}

cat("\n=== Summary table ===\n")
summary_df <- do.call(rbind, lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(
    kernel = nm,
    preload_same = r$preload_same,
    nmath_same = r$nmath_same,
    program_same = r$program_same,
    nchar_v1 = r$program_nchar_v1,
    nchar_v2 = r$program_nchar_v2,
    nchar_delta = r$program_nchar_v2 - r$program_nchar_v1,
    row.names = NULL
  )
}))
print(summary_df, row.names = FALSE)
