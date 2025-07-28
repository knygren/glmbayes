///*
// * OPENCL.CL – Global OpenCL Extensions & Shared Utilities
// * Stitch this at the top of every combined .cl module.
// */

#ifndef OPENCL_CL
#define OPENCL_CL

// -----------------------------------------------------------------------------
// 1. Enable OpenCL extensions
//    - cl_khr_fp64   → double-precision math
//    - cl_khr_printf → on-device printf debugging
// -----------------------------------------------------------------------------
#pragma OPENCL EXTENSION cl_khr_fp64   : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable

// -----------------------------------------------------------------------------
// 2. Standard headers
// -----------------------------------------------------------------------------
//#include <math.h>

// -----------------------------------------------------------------------------
// 3. IEEE constants
//    Mirrors R's ML_NAN / INF definitions for use throughout.
// -----------------------------------------------------------------------------
#ifndef ML_NAN
  #define ML_NAN      nan("")      /* quiet NaN */
  #define ML_POSINF   INFINITY     /* +∞ */
  #define ML_NEGINF  -INFINITY     /* –∞ */
#endif

// -----------------------------------------------------------------------------
// 4. Basic helpers
//
//    INLINE    → use for device‐side inlining
//    R_UNUSED → silence unused‐parameter warnings
// -----------------------------------------------------------------------------
#ifndef INLINE
  #define INLINE static inline
#endif

#ifndef R_UNUSED
  #define R_UNUSED(x) (void)(x)
#endif

// -----------------------------------------------------------------------------
// 5. Placeholder for other global utilities (e.g., work-item macros, typedefs)
// -----------------------------------------------------------------------------

#endif // OPENCL_CL