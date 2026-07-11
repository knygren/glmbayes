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
  cat("=== glmbayes OpenCL Diagnostic Report ===\n")

  info     <- opencltools::detect_environment_and_gpus()
  drivers  <- opencltools::detect_or_install_gpu_drivers(info)
  runtimes <- opencltools::detect_compute_runtimes(info)
  env_diag <- opencltools::check_runtime_env(runtimes)

  cat("Environment:", info$environment, "\n\n")

  gpu_vendor <- if (info$nvidia$present) "nvidia"
  else if (info$amd$present) "amd"
  else if (info$intel$present) "intel"
  else NULL

  diag <- NULL
  runtime_ok <- NA

  if (!is.null(gpu_vendor)) {
    cat("GPU:", toupper(gpu_vendor), "\n")
    drv  <- drivers$drivers[[gpu_vendor]]
    rt   <- runtimes$runtimes[[gpu_vendor]]
    diag <- env_diag$diagnostics[[gpu_vendor]]

    if (drv$installed) {
      cat("  [OK] Driver installed\n")
    } else {
      cat("  [FAIL] Driver not installed\n")
      if (length(drv$issues) > 0)
        cat("    Issues:", paste(drv$issues, collapse=", "), "\n")
    }

    hdr <- rt$opencl$headers_present
    rtm <- rt$opencl$runtime_present
    inst <- rt$opencl$installed

    if (hdr) {
      cat("  [OK] OpenCL headers found (CL/cl.h)\n")
    } else {
      cat("  [FAIL] OpenCL headers not found (CL/cl.h missing)\n")
    }

    if (rtm) {
      cat("  [OK] OpenCL runtime found (OpenCL.dll / ICD)\n")
    } else {
      cat("  [FAIL] OpenCL runtime not found\n")
    }

    if (inst) {
      cat("  [OK] OpenCL fully available (headers + runtime)\n")
    } else {
      cat("  [FAIL] OpenCL incomplete (missing headers or runtime)\n")
    }

    paths_ok <- (length(diag$opencl$missing_path_dirs) == 0 &&
                   length(diag$opencl$missing_lib_dirs) == 0)

    if (paths_ok) {
      cat("  [OK] Required PATH and library dirs present\n")
    } else {
      if (length(diag$opencl$missing_path_dirs) > 0)
        cat("  [WARN] Missing PATH entries:",
            paste(diag$opencl$missing_path_dirs, collapse=", "), "\n")
      if (length(diag$opencl$missing_lib_dirs) > 0)
        cat("  [WARN] Missing library dirs:",
            paste(diag$opencl$missing_lib_dirs, collapse=", "), "\n")
    }

    if (paths_ok && tolower(info$environment) %in% c("linux", "wsl")) {
      runtime_ok <- opencltools::verify_opencl_runtime(rt$opencl$lib_dirs)
      if (runtime_ok) {
        cat("  [OK] OpenCL runtime probe succeeded (platform available)\n")
      } else {
        cat("  [FAIL] OpenCL runtime probe failed (no usable platform)\n")
      }
    } else if (!paths_ok) {
      cat("  [SKIP] Runtime probe skipped (missing PATH/lib dirs)\n")
    } else {
      cat("  [SKIP] Runtime probe skipped on Windows\n")
    }

  } else {
    cat("[FAIL] No supported GPU detected. glmbayes will run in CPU-only mode.\n")
  }

  opencl_enabled <- has_opencl()
  if (opencl_enabled) {
    cat("\n[OK] glmbayes was compiled with OpenCL support.\n")
  } else {
    cat("\n[FAIL] glmbayes was compiled without OpenCL support.\n")
  }

  missing_items <- !is.null(diag) &&
    (length(diag$opencl$missing_path_dirs) > 0 ||
       length(diag$opencl$missing_lib_dirs) > 0)

  if (missing_items && !isTRUE(opencl_enabled)) {
    cat("\n[INFO] Missing PATH/lib entries detected and OpenCL is not enabled.\n")

    if (length(diag$opencl$missing_path_dirs) > 0) {
      cat("  Missing PATH entries:\n")
      cat("   -", paste(diag$opencl$missing_path_dirs, collapse="\n   - "), "\n")
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
      cat("  Missing library dirs:\n")
      cat("   -", paste(diag$opencl$missing_lib_dirs, collapse="\n   - "), "\n")
      ans_lib <- readline("Would you like to permanently add missing library dirs to LD_LIBRARY_PATH? [y/N]: ")
      if (tolower(ans_lib) == "y") {
        opencltools::add_to_libpath_linux(diag$opencl$missing_lib_dirs)
      }
    }
  }

  cat("\n=== End of Diagnostic Report ===\n")

  invisible(list(
    environment_info      = info,
    driver_status         = drivers,
    runtime_status        = runtimes,
    env_diag              = env_diag,
    opencl_runtime_probe  = runtime_ok,
    opencl_enabled        = opencl_enabled
  ))
}


#' @export
#' @rdname gpu_diagnostics
#' @order 2
has_opencl <- function() {
  .has_opencl_cpp()
}
