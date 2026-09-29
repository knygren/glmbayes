#' GPU and OpenCL diagnostics for \pkg{glmbayes}
#'
#' @description
#' Compile-time OpenCL status comes from **glmbayesCore** (sampler backend).
#' \code{\link{diagnose_glmbayes}()} is re-exported from Core; \code{\link{has_opencl}()}
#' reports whether that backend was built with OpenCL support.
#'
#' Host/runtime probes (GPU vendor, drivers, ICD/PATH) live in \pkg{opencltools}.
#'
#' @return
#' \code{has_opencl()} returns a length-1 logical.
#'
#' @seealso \code{\link{diagnose_glmbayes}}, \pkg{opencltools}, \pkg{glmbayesCore},
#'   \code{\link{rglmb}}, \code{\link{rlmb}}.
#' @name gpu_diagnostics
NULL

#' @inherit glmbayesCore::diagnose_glmbayes return title description details format note references examples
#' @export
#' @rdname gpu_diagnostics
#' @order 1
diagnose_glmbayes <- function() {
  glmbayesCore::diagnose_glmbayes()
}

#' @export
#' @rdname gpu_diagnostics
#' @order 2
has_opencl <- function() {
  glmbayesCore::glmbayesCore_has_opencl()
}
