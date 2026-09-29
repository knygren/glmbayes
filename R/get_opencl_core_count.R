#' Get the number of available OpenCL compute units
#'
#' Delegates to \pkg{opencltools} (same policy as **glmbayesCore**).
#'
#' @return Integer count of compute units.
#' @export
get_opencl_core_count <- function() {
  opencltools::get_opencl_core_count()
}
