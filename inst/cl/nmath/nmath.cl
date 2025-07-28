// Stub out R’s error/warning hooks for OpenCL
#ifndef MATHLIB_ERROR
  #define MATHLIB_ERROR(fmt, ...)     /* no-op or printf(fmt, __VA_ARGS__) */
#endif

#ifndef MATHLIB_WARNING
  #define MATHLIB_WARNING(fmt, ...)   /* no-op or printf(fmt, __VA_ARGS__) */
#endif

#ifdef DEBUG_MODE
  #undef MATHLIB_WARNING
  #define MATHLIB_WARNING(fmt, ...) printf("Rmath warning: " fmt, __VA_ARGS__)
#endif


// OpenCL-safe R_forceint equivalent
#define R_forceint(x) round(x)



#ifndef ML_ERROR
  #define ML_ERROR(code, fname)                                       \
    do {                                                               \
      const char *_msg = "";                                           \
      switch(code) {                                                   \
        case ME_DOMAIN:    _msg = "argument out of domain in %s\n";    break; \
        case ME_RANGE:     _msg = "value out of range in %s\n";        break; \
        case ME_NOCONV:    _msg = "convergence failed in %s\n";        break; \
        case ME_PRECISION: _msg = "precision lost in %s\n";            break; \
        case ME_UNDERFLOW: _msg = "underflow occurred in %s\n";        break; \
        default:          _msg = "math error %d in %s\n";              break; \
      }                                                                \
      MATHLIB_WARNING(_msg, fname);                                    \
    } while(0)
#endif


// nmath.cl - OpenCL math constants, macros & remaps for GPU kernels
//@provides: ML_POSINF,ML_NEGINF,ML_NAN,ME_NONE,ME_DOMAIN,ME_RANGE,ME_NOCONV,ME_PRECISION,ME_UNDERFLOW,ISNAN,R_FINITE,ML_VALID,ML_ERR_return_NAN,ML_ERROR,WILCOX_MAX,_

// ──────────────────────────────
// 📦 Numerical Constants (Quiet Override)
// ──────────────────────────────

#ifdef ML_POSINF
#undef ML_POSINF
#endif
#define ML_POSINF   (1.0 / 0.0)

#ifdef ML_NEGINF
#undef ML_NEGINF
#endif
#define ML_NEGINF   (-1.0 / 0.0)

#ifdef ML_NAN
#undef ML_NAN
#endif
#define ML_NAN      (0.0 / 0.0)


// Error Codes
#define ME_NONE        0
#define ME_DOMAIN      1
#define ME_RANGE       2
#define ME_NOCONV      4
#define ME_PRECISION   8
#define ME_UNDERFLOW  16

// Validation Macros
#ifdef ISNAN
#undef ISNAN
#endif
#define ISNAN(x) (isnan(x))

#ifdef R_FINITE
#undef R_FINITE
#endif
#define R_FINITE(x) (isfinite(x))

#ifdef ML_VALID
#undef ML_VALID
#endif
#define ML_VALID(x) (R_FINITE(x) && !ISNAN(x))

// Error-handling Macros

#undef ML_ERR_return_NAN
#define ML_ERR_return_NAN                          \
  do {                                             \
    ML_ERROR(ME_DOMAIN, fname);                    \
    return ML_NAN;                                 \
  } while(0)



#ifndef ML_WARN_return_NAN
  #define ML_WARN_return_NAN                            \
    do {                                               \
      ML_WARNING(ME_DOMAIN, fname);                    \
      return ML_NAN;                                   \
    } while(0)
#endif


// Wilcoxon Distribution Constant
#define WILCOX_MAX 50

// gettext stub
#define _(String) (String)