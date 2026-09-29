#' GPU and OpenCL diagnostics for \pkg{glmbayes}
#'
#' @description
#' Compile-time OpenCL probing for \pkg{glmbayes}, plus \code{diagnose_glmbayes()}
#' --- a readable report that combines \pkg{opencltools} host/runtime checks with
#' this package's OpenCL build status (\code{\link{has_opencl}}).
#'
#' Workstation probes (GPU vendor detection, drivers, ICD/PATH, and related helpers)
#' live in \pkg{opencltools}; call \code{opencltools::…} or see \code{?opencltools}.
#'
#' @section Diagnostics exported from \pkg{glmbayes}:
#' \itemize{
#'   \item \code{\link{diagnose_glmbayes}()} --- full report including compile-time
#'     OpenCL status for this package.
#'   \item \code{\link{has_opencl}()} --- \code{TRUE} if this build was compiled
#'     with OpenCL support.
#' }
#'
#' @section Host / runtime checks (\pkg{opencltools}):
#' \itemize{
#'   \item \code{\link[opencltools:gpu_diagnostics]{detect_environment_and_gpus}()}
#'   \item \code{\link[opencltools:gpu_diagnostics]{detect_compute_runtimes}()}
#'   \item \code{\link[opencltools:gpu_diagnostics]{verify_opencl_runtime}()}
#'   \item \code{\link[opencltools:gpu_diagnostics]{check_runtime_env}()}
#'   \item \code{\link[opencltools:load_kernel_source]{load_kernel_source}()},
#'     \code{\link[opencltools:load_kernel_source]{load_kernel_library}()}
#'     (pass \code{package = "glmbayes"} for entry kernels under \code{inst/cl/src/};
#'     prelude and \code{nmath/} come from \pkg{nmathopencl} in the production GPU path)
#'   \item \code{\link[opencltools:add_to_path]{add_to_path_windows}()} and related PATH helpers
#' }
#'
#' @details
#' GPU acceleration speeds up envelope construction and grid evaluation when you pass
#' \code{use_opencl = TRUE} in \code{\link{glmb}} and related functions. CPU-only
#' builds remain fully usable for standard modelling.
#'
#' Start with \code{\link{diagnose_glmbayes}()} for a single readable report;
#' use \code{\link{has_opencl}()} for a quick boolean when scripting. Full install
#' notes: \code{vignette("Chapter-16", package = "glmbayes")}
#' (\insertCite{glmbayesChapter12}{glmbayes}).
#'
#' @return
#' \code{diagnose_glmbayes()} returns an object of class
#' \code{"diagnose_glmbayes"} (also inheriting from \code{"list"}) with components
#' \code{environment_info}, \code{driver_status}, \code{runtime_status},
#' \code{env_diag}, \code{opencl_runtime_probe} (logical or \code{NA}), and
#' \code{opencl_enabled} (logical). A human-readable report is produced by
#' \code{\link{print.diagnose_glmbayes}} when the result is printed (for example
#' at the top level in an interactive session). Assigned results stay quiet until
#' \code{print()} is called.
#'
#' \code{has_opencl()} returns a length-1 logical: \code{TRUE}
#' if this build was compiled with OpenCL support.
#'
#' @seealso
#' \code{\link{diagnose_glmbayes}}, \code{\link{has_opencl}}, \pkg{opencltools},
#' \code{\link{glmb}}, \code{\link{rglmb}}.
#'
#' @references
#' \insertAllCited{}
#' @example inst/examples/Ex_load_kernel_source.R
#' @importFrom Rdpack reprompt
#' @keywords diagnostics gpu opencl environment
#' @name gpu_diagnostics
NULL


#' @export
#' @rdname gpu_diagnostics
#' @order 1
diagnose_glmbayes <- function() {
  info     <- opencltools::detect_environment_and_gpus()
  drivers  <- opencltools::detect_or_install_gpu_drivers(info)
  runtimes <- opencltools::detect_compute_runtimes(info)
  env_diag <- opencltools::check_runtime_env(runtimes)

  gpu_vendor <- if (info$nvidia$present) "nvidia"
  else if (info$amd$present) "amd"
  else if (info$intel$present) "intel"
  else NULL

  diag <- NULL
  runtime_ok <- NA

  if (!is.null(gpu_vendor)) {
    drv  <- drivers$drivers[[gpu_vendor]]
    rt   <- runtimes$runtimes[[gpu_vendor]]
    diag <- env_diag$diagnostics[[gpu_vendor]]

    paths_ok <- (length(diag$opencl$missing_path_dirs) == 0 &&
                   length(diag$opencl$missing_lib_dirs) == 0)

    if (paths_ok && tolower(info$environment) %in% c("linux", "wsl")) {
      runtime_ok <- opencltools::verify_opencl_runtime(rt$opencl$lib_dirs)
    }
  }

  opencl_enabled <- has_opencl()

  missing_items <- !is.null(diag) &&
    (length(diag$opencl$missing_path_dirs) > 0 ||
       length(diag$opencl$missing_lib_dirs) > 0)

  if (missing_items && !isTRUE(opencl_enabled) && interactive()) {
    message("Missing PATH/lib entries detected and OpenCL is not enabled.")

    if (length(diag$opencl$missing_path_dirs) > 0) {
      message("Missing PATH entries:")
      message(paste(" -", diag$opencl$missing_path_dirs, collapse = "\n"))
      ans <- readline("Would you like to permanently add missing PATH dirs? [y/N]: ")
      if (tolower(ans) == "y") {
        if (tolower(info$environment) == "windows") {
          opencltools::add_to_path_windows(diag$opencl$missing_path_dirs)
        } else {
          opencltools::add_to_path_linux(diag$opencl$missing_path_dirs)
        }
      }
    }

    if (length(diag$opencl$missing_lib_dirs) > 0 &&
        tolower(info$environment) %in% c("linux", "wsl")) {
      message("Missing library dirs:")
      message(paste(" -", diag$opencl$missing_lib_dirs, collapse = "\n"))
      ans_lib <- readline("Would you like to permanently add missing library dirs to LD_LIBRARY_PATH? [y/N]: ")
      if (tolower(ans_lib) == "y") {
        opencltools::add_to_libpath_linux(diag$opencl$missing_lib_dirs)
      }
    }
  }

  out <- list(
    environment_info      = info,
    driver_status         = drivers,
    runtime_status        = runtimes,
    env_diag              = env_diag,
    opencl_runtime_probe  = runtime_ok,
    opencl_enabled        = opencl_enabled
  )
  class(out) <- c("diagnose_glmbayes", "list")
  out
}


## print.* S3 method provided by glmbayesCore (Imports); local registration commented out.
# #' @export
# #' @rdname gpu_diagnostics
# #' @order 2
# #' @param x An object of class \code{"diagnose_glmbayes"}.
# #' @param ... Unused; for S3 method compatibility.
# print.diagnose_glmbayes <- function(x, ...) {
#   info <- x$environment_info
#   drivers <- x$driver_status
#   runtimes <- x$runtime_status
#   env_diag <- x$env_diag
#   runtime_ok <- x$opencl_runtime_probe
#   opencl_enabled <- x$opencl_enabled
# 
#   cat("=== glmbayes OpenCL Diagnostic Report ===\n")
#   cat("Environment:", info$environment, "\n\n")
# 
#   gpu_vendor <- if (isTRUE(info$nvidia$present)) "nvidia"
#   else if (isTRUE(info$amd$present)) "amd"
#   else if (isTRUE(info$intel$present)) "intel"
#   else NULL
# 
#   if (!is.null(gpu_vendor)) {
#     cat("GPU:", toupper(gpu_vendor), "\n")
#     drv  <- drivers$drivers[[gpu_vendor]]
#     rt   <- runtimes$runtimes[[gpu_vendor]]
#     diag <- env_diag$diagnostics[[gpu_vendor]]
# 
#     if (isTRUE(drv$installed)) {
#       cat("  [OK] Driver installed\n")
#     } else {
#       cat("  [FAIL] Driver not installed\n")
#       if (length(drv$issues) > 0)
#         cat("    Issues:", paste(drv$issues, collapse = ", "), "\n")
#     }
# 
#     hdr <- isTRUE(rt$opencl$headers_present)
#     rtm <- isTRUE(rt$opencl$runtime_present)
#     inst <- isTRUE(rt$opencl$installed)
# 
#     if (hdr) {
#       cat("  [OK] OpenCL headers found (CL/cl.h)\n")
#     } else {
#       cat("  [FAIL] OpenCL headers not found (CL/cl.h missing)\n")
#     }
# 
#     if (rtm) {
#       cat("  [OK] OpenCL runtime found (OpenCL.dll / ICD)\n")
#     } else {
#       cat("  [FAIL] OpenCL runtime not found\n")
#     }
# 
#     if (inst) {
#       cat("  [OK] OpenCL fully available (headers + runtime)\n")
#     } else {
#       cat("  [FAIL] OpenCL incomplete (missing headers or runtime)\n")
#     }
# 
#     paths_ok <- (length(diag$opencl$missing_path_dirs) == 0 &&
#                    length(diag$opencl$missing_lib_dirs) == 0)
# 
#     if (paths_ok) {
#       cat("  [OK] Required PATH and library dirs present\n")
#     } else {
#       if (length(diag$opencl$missing_path_dirs) > 0)
#         cat("  [WARN] Missing PATH entries:",
#             paste(diag$opencl$missing_path_dirs, collapse = ", "), "\n")
#       if (length(diag$opencl$missing_lib_dirs) > 0)
#         cat("  [WARN] Missing library dirs:",
#             paste(diag$opencl$missing_lib_dirs, collapse = ", "), "\n")
#     }
# 
#     if (paths_ok && tolower(info$environment) %in% c("linux", "wsl")) {
#       if (isTRUE(runtime_ok)) {
#         cat("  [OK] OpenCL runtime probe succeeded (platform available)\n")
#       } else if (isFALSE(runtime_ok)) {
#         cat("  [FAIL] OpenCL runtime probe failed (no usable platform)\n")
#       } else {
#         cat("  [SKIP] Runtime probe result unavailable\n")
#       }
#     } else if (!paths_ok) {
#       cat("  [SKIP] Runtime probe skipped (missing PATH/lib dirs)\n")
#     } else {
#       cat("  [SKIP] Runtime probe skipped on Windows\n")
#     }
# 
#   } else {
#     cat("[FAIL] No supported GPU detected. glmbayes will run in CPU-only mode.\n")
#   }
# 
#   if (isTRUE(opencl_enabled)) {
#     cat("\n[OK] glmbayes was compiled with OpenCL support.\n")
#   } else {
#     cat("\n[FAIL] glmbayes was compiled without OpenCL support.\n")
#   }
# 
#   cat("\n=== End of Diagnostic Report ===\n")
#   invisible(x)
# }


#' @export
#' @rdname gpu_diagnostics
#' @order 3
has_opencl <- function() {
  .has_opencl_cpp()
}
