//@provides: ML_POSINF,ML_NEGINF,ML_NAN,ME_NONE,ME_DOMAIN,ME_RANGE,ME_NOCONV,ME_PRECISION,ME_UNDERFLOW,ISNAN,R_FINITE,ML_VALID,ML_ERR_return_NAN,ML_ERROR,WILCOX_MAX,_

// Numerical Constants
#define ML_POSINF   (1.0 / 0.0)
#define ML_NEGINF   (-1.0 / 0.0)
#define ML_NAN      (0.0 / 0.0)

// Error Codes
#define ME_NONE        0
#define ME_DOMAIN      1
#define ME_RANGE       2
#define ME_NOCONV      4
#define ME_PRECISION   8
#define ME_UNDERFLOW  16

// Validation Macros
#define ISNAN(x)       (isnan(x))
#define R_FINITE(x)    (isfinite(x))
#define ML_VALID(x)    (!ISNAN(x))

// Error-handling Macros
#define ML_ERR_return_NAN return ML_NAN
#define ML_ERROR(code, msg) /* stubbed macro for symbolic error messaging */

// Wilcoxon Distribution Constant
#define WILCOX_MAX 50

// gettext stub
#define _(String) (String)