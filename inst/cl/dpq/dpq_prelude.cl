// @provides ML_NEGINF, ML_POSINF, ML_NAN, give_log, R_FINITE, R_nonint, safe_log, safe_log1p, safe_expm1

// OpenCL-safe definitions for constants - No longer needed as set on OPENCL.CL
//#define ML_NEGINF   (-INFINITY)
//#define ML_POSINF   (INFINITY)
//#define ML_NAN      (NAN)

// Shim for give_log logic (set externally)
#define give_log    log_p

// Compatibility fallback for R_FINITE
#define R_FINITE(x) isfinite(x)

// Check for non-integer values (float-safe)
inline bool R_nonint(float x) {
  return floor(x) != x;
}

// Shim error macros (replace R-side macros)
#define ML_WARN_return_NAN return NAN;

// Guard values against domain errors or log violations
inline float safe_log(float x) {
  return (x > 0.0f) ? log(x) : ML_NEGINF;
}

inline float safe_log1p(float x) {
  return ((1.0f + x) > 0.0f) ? log1p(x) : ML_NEGINF;
}

inline float safe_expm1(float x) {
  return exp(x) - 1.0f;
}