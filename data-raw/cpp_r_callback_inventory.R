#!/usr/bin/env Rscript
## Scan src/*.cpp for C++ -> R lookups (Rcpp::Function literals, namespace_env).
## Whitelist must match data-raw/CPP_R_CALLBACK_INVENTORY.md (sections 1-3).
##
## Usage:
##   Rscript data-raw/cpp_r_callback_inventory.R
##   Rscript data-raw/cpp_r_callback_inventory.R C:/Rpackages/glmbayesCore
##
## Exit 1 if an unknown symbol literal appears.

args <- commandArgs(trailingOnly = TRUE)

resolve_pkg_root <- function() {
  if (length(args) >= 1L && nzchar(args[[1L]])) {
    p <- normalizePath(args[[1L]], mustWork = TRUE)
    if (!file.exists(file.path(p, "DESCRIPTION"))) {
      stop("Not an R package root (no DESCRIPTION): ", p)
    }
    return(p)
  }
  if (file.exists("DESCRIPTION")) {
    return(normalizePath("."))
  }
  if (file.exists("../DESCRIPTION")) {
    return(normalizePath(".."))
  }
  stop("Pass package root as first argument or run from package / data-raw directory.")
}

pkg_root <- resolve_pkg_root()
src_dir <- file.path(pkg_root, "src")
if (!dir.exists(src_dir)) {
  stop("No src/ directory under ", pkg_root)
}

cpp_files <- list.files(src_dir, pattern = "\\.cpp$", full.names = TRUE, recursive = TRUE)
cpp_files <- cpp_files[!grepl("[/\\\\]backup[/\\\\]", cpp_files, ignore.case = TRUE)]

strip_cpp_line_comments <- function(lines) {
  vapply(lines, function(ln) {
    if (grepl("//", ln, fixed = TRUE)) {
      sub("//.*$", "", ln)
    } else {
      ln
    }
  }, character(1))
}

read_text <- function(path) {
  lines <- readLines(path, warn = FALSE)
  paste(strip_cpp_line_comments(lines), collapse = "\n")
}

fn_literal_re <- '(?:Rcpp::)?Function(?:\\s+\\w+)?\\s*\\(\\s*"([^"]+)"\\s*\\)'
ns_env_re <- 'namespace_env\\s*\\(\\s*"[^"]+"\\s*\\)\\s*\\[\\s*"([^"]+)"\\s*\\]'
pkg_bracket_re <- '(?:pkg|glmbayes_ns)\\s*\\[\\s*"([^"]+)"\\s*\\]'

known_base_r <- c(
  "optim", "try", "gaussian", "lm.wfit", "lm.fit",
  "as.matrix", "as.vector", "as.numeric",
  "expand.grid", "qgamma", "runif",
  "interactive", "readline",
  "format", "Sys.time", "system.file"
)

known_pkg_local <- c(
  "EnvelopeOpt", "EnvelopeSort", "glmbfamfunc",
  "rNormal_reg.wfit", "rgamma_ct"
)

known_other_pkg <- c("get_opencl_core_count")

known_all <- c(known_base_r, known_pkg_local, known_other_pkg)

norm_path <- function(f) {
  gsub("\\\\", "/", f)
}

rel_path <- function(f) {
  root <- gsub("\\\\", "/", pkg_root)
  sub(paste0("^", root, "/?"), "", norm_path(f))
}

hits <- list()
for (f in cpp_files) {
  txt <- read_text(f)
  rel <- rel_path(f)

  m <- gregexpr(fn_literal_re, txt, perl = TRUE)[[1]]
  if (m[1] != -1) {
    for (i in seq_along(m)) {
      start <- m[i]
      len <- attr(m, "match.length")[i]
      chunk <- substr(txt, start, start + len - 1)
      nm <- sub(fn_literal_re, "\\1", chunk, perl = TRUE)
      hits[[length(hits) + 1L]] <- list(
        file = rel, symbol = nm, kind = "Rcpp::Function literal"
      )
    }
  }

  m_pkg <- gregexpr(pkg_bracket_re, txt, perl = TRUE)[[1]]
  if (m_pkg[1] != -1) {
    for (i in seq_along(m_pkg)) {
      start <- m_pkg[i]
      len <- attr(m_pkg, "match.length")[i]
      chunk <- substr(txt, start, start + len - 1)
      nm <- sub(pkg_bracket_re, "\\1", chunk, perl = TRUE)
      hits[[length(hits) + 1L]] <- list(
        file = rel, symbol = nm, kind = "pkg_env[]"
      )
    }
  }

  m2 <- gregexpr(ns_env_re, txt, perl = TRUE)[[1]]
  if (m2[1] != -1) {
    for (i in seq_along(m2)) {
      start <- m2[i]
      len <- attr(m2, "match.length")[i]
      chunk <- substr(txt, start, start + len - 1)
      nm <- sub(ns_env_re, "\\1", chunk, perl = TRUE)
      hits[[length(hits) + 1L]] <- list(
        file = rel, symbol = nm, kind = "namespace_env[]"
      )
    }
  }

  if (grepl("glmbayes_R::", txt, fixed = TRUE)) {
    hits[[length(hits) + 1L]] <- list(
      file = rel,
      symbol = "(glmbayes_R accessor)",
      kind = "glmbayes_R:: (see R_interface.h)"
    )
  }
}

cat("=== C++ -> R callback scan:", pkg_root, "===\n")
cat("Scanned", length(cpp_files), "cpp file(s) under src/\n\n")

if (!length(hits)) {
  cat("No Rcpp::Function literals or namespace_env[] lookups found.\n")
  quit(status = 0)
}

sym_df <- unique(do.call(rbind, lapply(hits, function(h) {
  data.frame(
    file = h$file, symbol = h$symbol, kind = h$kind,
    stringsAsFactors = FALSE
  )
})))

print(sym_df, row.names = FALSE)

unknown <- sym_df[
  !sym_df$symbol %in% known_all &
    sym_df$symbol != "(glmbayes_R accessor)",
]

if (nrow(unknown)) {
  cat("\n*** UNKNOWN symbols (update CPP_R_CALLBACK_INVENTORY.md and whitelist): ***\n")
  print(unique(unknown[, c("symbol", "file")]), row.names = FALSE)
  quit(status = 1)
}

cat("\nAll literal lookups match the canonical inventory.\n")
quit(status = 0)
