// @depends dpq_prelude
// @provides R_D__0,R_D__1,R_DT_0,R_DT_1,R_D_half,R_D_Lval,R_D_Cval,R_D_val,R_D_qIv,R_D_exp,R_D_log,R_D_Clog,R_Log1_Exp,R_D_LExp,R_DT_val,R_DT_Cval,R_DT_qIv,R_DT_CIv,R_DT_exp,R_DT_Cexp,R_DT_log,R_DT_Clog,R_DT_Log,R_Q_P01_check,R_Q_P01_boundaries,R_P_bounds_01,R_P_bounds_Inf_01,R_D_fexp,R_D_negInonint,R_D_nonint_check

// --- DPQ Macros for Distribution and Tail Handling ---

#define R_D__0      (log_p ? ML_NEGINF : 0.0f)
#define R_D__1      (log_p ? 0.0f : 1.0f)
#define R_DT_0      (lower_tail ? R_D__0 : R_D__1)
#define R_DT_1      (lower_tail ? R_D__1 : R_D__0)
#define R_D_half    (log_p ? -M_LN2 : 0.5f)

#define R_D_Lval(p)     (lower_tail ? (p) : (0.5f - (p) + 0.5f))
#define R_D_Cval(p)     (lower_tail ? (0.5f - (p) + 0.5f) : (p))

#define R_D_val(x)      (log_p ? safe_log(x) : (x))
#define R_D_qIv(p)      (log_p ? exp(p) : (p))
#define R_D_exp(x)      (log_p ? (x) : exp(x))
#define R_D_log(p)      (log_p ? (p) : safe_log(p))
#define R_D_Clog(p)     (log_p ? safe_log1p(-(p)) : (0.5f - (p) + 0.5f))

#define R_Log1_Exp(x)   ((x) > -M_LN2 ? safe_log(-safe_expm1(x)) : safe_log1p(-exp(x)))
#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : safe_log1p(-x))

#define R_DT_val(x)     (lower_tail ? R_D_val(x) : R_D_Clog(x))
#define R_DT_Cval(x)    (lower_tail ? R_D_Clog(x) : R_D_val(x))

#define R_DT_qIv(p)     (log_p ? (lower_tail ? exp(p) : -safe_expm1(p)) : R_D_Lval(p))
#define R_DT_CIv(p)     (log_p ? (lower_tail ? -safe_expm1(p) : exp(p)) : R_D_Cval(p))

#define R_DT_exp(x)     R_D_exp(R_D_Lval(x))
#define R_DT_Cexp(x)    R_D_exp(R_D_Cval(x))

#define R_DT_log(p)     (lower_tail ? R_D_log(p) : R_D_LExp(p))
#define R_DT_Clog(p)    (lower_tail ? R_D_LExp(p) : R_D_log(p))
#define R_DT_Log(p)     (lower_tail ? (p) : R_Log1_Exp(p))

// --- Input Boundary Macros ---
// Note, put everything in a singl line so that the kernel interprets correctly

#define R_Q_P01_check(p)    if ((log_p && p > 0.0f) ||(!log_p && (p < 0.0f || p > 1.0f))) ML_WARN_return_NAN

#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  if (log_p) {if (p > 0.0f) ML_WARN_return_NAN;if (p == 0.0f) return lower_tail ? _RIGHT_ : _LEFT_; if (p == ML_NEGINF) return lower_tail ? _LEFT_ : _RIGHT_;} else {if (p < 0.0f || p > 1.0f) ML_WARN_return_NAN; if (p == 0.0f) return lower_tail ? _LEFT_ : _RIGHT_; if (p == 1.0f) return lower_tail ? _RIGHT_ : _LEFT_; }

#define R_P_bounds_01(x, x_min, x_max)  if (x <= x_min) return R_DT_0;if (x >= x_max) return R_DT_1

#define R_P_bounds_Inf_01(x) if (!R_FINITE(x)) {if (x > 0.0f) return R_DT_1;return R_DT_0;}

// --- Density-Specific Helpers ---

#define R_D_fexp(f, x)      (give_log ? -0.5f * safe_log(f) + (x) : exp(x) / sqrt(f))
#define R_D_rtxp(rf, x)     (give_log ? -safe_log(rf) + (x) : exp(x) / (rf))

#define R_D_negInonint(x)   ((x) < 0.0f || R_nonint(x))

#define R_D_nonint_check(x) if (R_nonint(x)) {  MATHLIB_WARNING("non-integer x detected");return R_D__0;  }
    
