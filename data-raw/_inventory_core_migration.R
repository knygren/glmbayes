gb <- "C:/Rpackages/glmbayes"
gc <- "C:/Rpackages/glmbayesCore"
dirs <- c("R", "src", "inst/cl")
collect <- function(root, sub) {
  p <- file.path(root, sub)
  if (!dir.exists(p)) return(character())
  list.files(p, recursive = TRUE, full.names = TRUE)
}
rel <- function(root, f) sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", gsub("\\\\", "/", f))
gb_files <- unlist(lapply(dirs, function(d) collect(gb, d)))
gc_files <- unlist(lapply(dirs, function(d) collect(gc, d)))
gb_rel <- rel(gb, gb_files)
gc_rel <- rel(gc, gc_files)
all_rel <- sort(unique(c(gb_rel, gc_rel)))
same <- character()
diff <- character()
gb_only <- character()
gc_only <- character()
for (r in all_rel) {
  p1 <- file.path(gb, r)
  p2 <- file.path(gc, r)
  in1 <- file.exists(p1)
  in2 <- file.exists(p2)
  if (in1 && !in2) { gb_only <- c(gb_only, r); next }
  if (!in1 && in2) { gc_only <- c(gc_only, r); next }
  if (identical(readBin(p1, raw(), file.info(p1)$size),
                readBin(p2, raw(), file.info(p2)$size))) {
    same <- c(same, r)
  } else {
    diff <- c(diff, r)
  }
}
cat("SAME", length(same), "DIFF", length(diff), "GB_ONLY", length(gb_only), "GC_ONLY", length(gc_only), "\n")
writeLines(diff, "C:/Rpackages/glmbayes/data-raw/_diff_paths.txt")
writeLines(gb_only, "C:/Rpackages/glmbayes/data-raw/_gb_only_paths.txt")
writeLines(gc_only, "C:/Rpackages/glmbayes/data-raw/_gc_only_paths.txt")
