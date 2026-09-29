# Package attach hook: optional OpenCL advisory for end-user workflows (CPU-first).

.opencl_startup_quiet <- function() {
  isTRUE(getOption("glmbayes.quiet_opencl_startup", FALSE)) || !interactive()
}

.opencl_runtime_sniff <- function() {
  tryCatch(
    {
      info <- opencltools::detect_environment_and_gpus()
      rt_list <- opencltools::detect_compute_runtimes(info)
      gpu <- isTRUE(info$nvidia$present) ||
        isTRUE(info$amd$present) ||
        isTRUE(info$intel$present)
      stack_ok <- FALSE
      for (vendor in c("nvidia", "amd", "intel")) {
        ocl <- rt_list$runtimes[[vendor]]$opencl
        if (isTRUE(ocl$installed) ||
            (isTRUE(ocl$headers_present) && isTRUE(ocl$runtime_present))) {
          stack_ok <- TRUE
          break
        }
      }
      list(gpu = gpu, stack_ok = stack_ok)
    },
    error = function(e) {
      list(gpu = FALSE, stack_ok = FALSE)
    }
  )
}

.opencl_startup_message <- function() {
  if (.opencl_startup_quiet()) {
    return(invisible())
  }
  if (isTRUE(getOption("glmbayes.opencl_startup_checked", FALSE))) {
    return(invisible())
  }
  options(glmbayes.opencl_startup_checked = TRUE)

  if (isTRUE(tryCatch(has_opencl(), error = function(e) FALSE))) {
    return(invisible())
  }

  sniff <- .opencl_runtime_sniff()
  if (!isTRUE(sniff$stack_ok) && !isTRUE(sniff$gpu)) {
    return(invisible())
  }

  packageStartupMessage(
    "Note: glmbayes provides full CPU capability in this session ",
    "(e.g. glmb(), lmb(), Prior_Setup()). GPU acceleration is recommended ",
    "for bigger models and appears available. Reinstall glmbayesCore from ",
    "source with OpenCL at compile time to enable it; see vignette(\"Chapter-16\", ",
    "\"glmbayes\") for install instructions."
  )
  invisible()
}

.onAttach <- function(libname, pkgname) {
  .opencl_startup_message()
}
