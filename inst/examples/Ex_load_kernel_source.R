############################### Start of load_kernel_source example ####################

if (has_opencl()) {
  # Family entry kernels live in glmbayes; prelude/nmath in nmathopencl (production v2 path)
  ksrc <- opencltools::load_kernel_source(
    "src/f2_f3_binomial_logit.cl",
    package = "glmbayes"
  )
  nmath_stem <- opencltools::load_kernel_source(
    "nmath/bd0.cl",
    package = "nmathopencl"
  )
  nmath_lib <- opencltools::load_kernel_library(
    "nmath",
    package = "nmathopencl"
  )
  cat("Entry kernel length:", nchar(ksrc), "\n")
  cat("nmath stem length:", nchar(nmath_stem), "\n")
  cat("nmath library length:", nchar(nmath_lib), "\n")
} else {
  message("OpenCL not available in this build; skipping kernel load example.")
}

## End of load_kernel_source example
