// Rmath.cl — OpenCL math constants, macros & remaps for GPU kernels
// SPDX-License-Identifier: LGPL-2.1-or-later
// @provides: Rboolean, M_E, M_LOG2E, M_LOG10E, M_LN2, M_LN10, M_PI, M_2PI, M_PI_2, M_PI_4, M_1_PI, M_2_PI, M_2_SQRTPI, M_SQRT2, M_SQRT1_2, M_SQRT_3, M_SQRT_32, M_LOG10_2, M_SQRT_PI, M_1_SQRT_2PI, M_SQRT_2dPI, M_LN_2PI, M_LN_SQRT_PI, M_LN_SQRT_2PI, M_LN_SQRT_PId2, ISNAN, R_FINITE, ML_VALID,  bessel_i, bessel_j, bessel_k, bessel_y, bessel_i_ex, bessel_j_ex, bessel_k_ex, bessel_y_ex, beta, choose, dbeta, dbinom, dbinom_raw, dcauchy, dchisq, dexp, df, dgamma, dgeom, dhyper, digamma, dlnorm, dlogis, dnbeta, dnbinom, dnbinom_mu, dnchisq, dnf, dnorm4, dnt, dpois_raw, dpois, dpsifn, dsignrank, dt, dtukey, dunif, dweibull, dwilcox, fmax2, fmin2, fprec, fround, ftrunc, fsign, gammafn, imax2, imin2, lbeta, lchoose, lgammafn, lgammafn_sign, lgamma1p, log1pexp, log1pmx, logspace_add, logspace_sub, logspace_sum, pbeta, pbeta_raw, pbinom, pcauchy, pchisq, pentagamma, pexp, pf, pgamma, pgeom, phyper, plnorm, plogis, pnbeta, pnbinom, pnbinom_mu, pnchisq, pnf, pnorm5, pnorm_both, pnt, ppois, psignrank, psigamma, pt, ptukey, punif, pythag, pweibull, pwilcox, qbeta, qbinom, qcauchy, qchisq, qchisq_appr, qexp, qf, qgamma, qgeom, qhyper, qlnorm, qlogis, qnbeta, qnbinom, qnbinom_mu, qnchisq, qnf, qnorm5, qnt, qpois, qsignrank, qt, qtukey, qunif, qweibull, qwilcox, rbeta, rbinom, rcauchy, rchisq, rexp, rf, rgamma, rgeom, rhyper, rlnorm, rlogis, rmultinom, rnbeta, rnbinom, rnbinom_mu, rnchisq, rnf, rnorm, rnt, rpois, rsignrank, rt, rtukey, runif, rweibull, rwilcox, sign, tetragamma, trigamma

//#pragma OPENCL EXTENSION cl_khr_fp64 : enable

//─────────────────────────────────────────────────────────────────────────────
// 1) R‐style Boolean
//─────────────────────────────────────────────────────────────────────────────
typedef enum { FALSE = 0, TRUE = 1 } Rboolean;

//─────────────────────────────────────────────────────────────────────────────
// 2) 30‐decimal‐place mathematical constants (SVID, X/Open & R-specific)
//─────────────────────────────────────────────────────────────────────────────
// ───────────────────────────────────────────────────────────────
// Clean override of 30-decimal-place mathematical constants
// ───────────────────────────────────────────────────────────────

#ifdef M_E
#undef M_E
#endif
#define M_E              2.718281828459045235360287471353

#ifdef M_LOG2E
#undef M_LOG2E
#endif
#define M_LOG2E          1.442695040888963407359924681002

#ifdef M_LOG10E
#undef M_LOG10E
#endif
#define M_LOG10E         0.434294481903251827651128918917

#ifdef M_LN2
#undef M_LN2
#endif
#define M_LN2            0.693147180559945309417232121458

#ifdef M_LN10
#undef M_LN10
#endif
#define M_LN10           2.302585092994045684017991454684

#ifdef M_PI
#undef M_PI
#endif
#define M_PI             3.141592653589793238462643383280

#ifdef M_2PI
#undef M_2PI
#endif
#define M_2PI            6.283185307179586476925286766559

#ifdef M_PI_2
#undef M_PI_2
#endif
#define M_PI_2           1.570796326794896619231321691640

#ifdef M_PI_4
#undef M_PI_4
#endif
#define M_PI_4           0.785398163397448309615660845820

#ifdef M_1_PI
#undef M_1_PI
#endif
#define M_1_PI           0.318309886183790671537767526745

#ifdef M_2_PI
#undef M_2_PI
#endif
#define M_2_PI           0.636619772367581343075535053490

#ifdef M_2_SQRTPI
#undef M_2_SQRTPI
#endif
#define M_2_SQRTPI       1.128379167095512573896158903122

#ifdef M_SQRT2
#undef M_SQRT2
#endif
#define M_SQRT2          1.414213562373095048801688724210

#ifdef M_SQRT1_2
#undef M_SQRT1_2
#endif
#define M_SQRT1_2        0.707106781186547524400844362105

// ───────────────────────────────────────────────────────────────
// R-specific constants
// ───────────────────────────────────────────────────────────────

#ifdef M_SQRT_3
#undef M_SQRT_3
#endif
#define M_SQRT_3         1.732050807568877293527446341506

#ifdef M_SQRT_32
#undef M_SQRT_32
#endif
#define M_SQRT_32        5.656854249492380195206754896838

#ifdef M_LOG10_2
#undef M_LOG10_2
#endif
#define M_LOG10_2        0.301029995663981195213738894724

#ifdef M_SQRT_PI
#undef M_SQRT_PI
#endif
#define M_SQRT_PI        1.772453850905516027298167483341

#ifdef M_1_SQRT_2PI
#undef M_1_SQRT_2PI
#endif
#define M_1_SQRT_2PI     0.398942280401432677939946059934

#ifdef M_SQRT_2dPI
#undef M_SQRT_2dPI
#endif
#define M_SQRT_2dPI      0.797884560802865355879892119869

#ifdef M_LN_2PI
#undef M_LN_2PI
#endif
#define M_LN_2PI         1.837877066409345483560659472811

#ifdef M_LN_SQRT_PI
#undef M_LN_SQRT_PI
#endif
#define M_LN_SQRT_PI     0.572364942924700087071713675677

#ifdef M_LN_SQRT_2PI
#undef M_LN_SQRT_2PI
#endif
#define M_LN_SQRT_2PI    0.918938533204672741780329736406

#ifdef M_LN_SQRT_PId2
#undef M_LN_SQRT_PId2
#endif
#define M_LN_SQRT_PId2   0.225791352644727432363097614947
//─────────────────────────────────────────────────────────────────────────────
// 3) Validation macros (match R’s ISNAN, R_FINITE, etc.)
//─────────────────────────────────────────────────────────────────────────────
#define ISNAN(x)     isnan(x)
#define R_FINITE(x)  isfinite(x)
#define ML_VALID(x)  (R_FINITE(x) && !ISNAN(x))



//─────────────────────────────────────────────────────────────────────────────
// 5) Remapping “friendly” names → Rf_… backend
//    (identical to Rmath.h’s #defines; left intact for on‐host linking via
//     clCreateProgramWithSource when building against standalone Rmath)
//─────────────────────────────────────────────────────────────────────────────
#if !defined(MATHLIB_STANDALONE) && !defined(R_NO_REMAP_RMATH)

#define bessel_i        Rf_bessel_i
#define bessel_j        Rf_bessel_j
#define bessel_k        Rf_bessel_k
#define bessel_y        Rf_bessel_y
#define bessel_i_ex     Rf_bessel_i_ex
#define bessel_j_ex     Rf_bessel_j_ex
#define bessel_k_ex     Rf_bessel_k_ex
#define bessel_y_ex     Rf_bessel_y_ex

#define beta            Rf_beta
#define choose          Rf_choose

#define dbeta           Rf_dbeta
#define dbinom          Rf_dbinom
#define dbinom_raw      Rf_dbinom_raw
#define dcauchy         Rf_dcauchy
#define dchisq          Rf_dchisq
#define dexp            Rf_dexp
#define df              Rf_df
#define dgamma          Rf_dgamma
#define dgeom           Rf_dgeom
#define dhyper          Rf_dhyper
#define digamma         Rf_digamma
#define dlnorm          Rf_dlnorm
#define dlogis          Rf_dlogis
#define dnbeta          Rf_dnbeta
#define dnbinom         Rf_dnbinom
#define dnbinom_mu      Rf_dnbinom_mu
#define dnchisq         Rf_dnchisq
#define dnf             Rf_dnf
#define dnorm4          Rf_dnorm4
#define dnt             Rf_dnt
#define dpois_raw       Rf_dpois_raw
#define dpois           Rf_dpois
#define dpsifn          Rf_dpsifn
#define dsignrank       Rf_dsignrank
#define dt              Rf_dt
#define dtukey          Rf_dtukey
#define dunif           Rf_dunif
#define dweibull        Rf_dweibull
#define dwilcox         Rf_dwilcox

#define fmax2           Rf_fmax2
#define fmin2           Rf_fmin2
#define fprec           Rf_fprec
#define fround          Rf_fround
#define ftrunc          Rf_ftrunc
#define fsign           Rf_fsign

#define gammafn         Rf_gammafn

#define imax2           Rf_imax2
#define imin2           Rf_imin2

#define lbeta           Rf_lbeta
#define lchoose         Rf_lchoose
#define lgammafn        Rf_lgammafn
#define lgammafn_sign   Rf_lgammafn_sign
#define lgamma1p        Rf_lgamma1p

#define log1pexp        Rf_log1pexp
#define log1pmx         Rf_log1pmx
#define logspace_add    Rf_logspace_add
#define logspace_sub    Rf_logspace_sub
#define logspace_sum    Rf_logspace_sum

#define pbeta           Rf_pbeta
#define pbeta_raw       Rf_pbeta_raw
#define pbinom          Rf_pbinom
#define pcauchy         Rf_pcauchy
#define pchisq          Rf_pchisq
#define pentagamma      Rf_pentagamma
#define pexp            Rf_pexp
#define pf              Rf_pf
#define pgamma          Rf_pgamma
#define pgeom           Rf_pgeom
#define phyper          Rf_phyper
#define plnorm          Rf_plnorm
#define plogis          Rf_plogis
#define pnbeta          Rf_pnbeta
#define pnbinom         Rf_pnbinom
#define pnbinom_mu      Rf_pnbinom_mu
#define pnchisq         Rf_pnchisq
#define pnf             Rf_pnf
#define pnorm5          Rf_pnorm5
#define pnorm_both      Rf_pnorm_both
#define pnt             Rf_pnt
#define ppois           Rf_ppois
#define psignrank       Rf_psignrank
#define psigamma        Rf_psigamma
#define pt              Rf_pt
#define ptukey          Rf_ptukey
#define punif           Rf_punif
#define pythag          Rf_pythag
#define pweibull        Rf_pweibull
#define pwilcox         Rf_pwilcox

#define qbeta           Rf_qbeta
#define qbinom          Rf_qbinom
#define qcauchy         Rf_qcauchy
#define qchisq          Rf_qchisq
#define qchisq_appr     Rf_qchisq_appr
#define qexp            Rf_qexp
#define qf              Rf_qf
#define qgamma          Rf_qgamma
#define qgeom           Rf_qgeom
#define qhyper          Rf_qhyper
#define qlnorm          Rf_qlnorm
#define qlogis          Rf_qlogis
#define qnbeta          Rf_qnbeta
#define qnbinom         Rf_qnbinom
#define qnbinom_mu      Rf_qnbinom_mu
#define qnchisq         Rf_qnchisq
#define qnf             Rf_qnf
#define qnorm5          Rf_qnorm5
#define qnt             Rf_qnt
#define qpois           Rf_qpois
#define qsignrank       Rf_qsignrank
#define qt              Rf_qt
#define qtukey          Rf_qtukey
#define qunif           Rf_qunif
#define qweibull        Rf_qweibull
#define qwilcox         Rf_qwilcox

#define rbeta           Rf_rbeta
#define rbinom          Rf_rbinom
#define rcauchy         Rf_rcauchy
#define rchisq          Rf_rchisq
#define rexp            Rf_rexp
#define rf              Rf_rf
#define rgamma          Rf_rgamma
#define rgeom           Rf_rgeom
#define rhyper          Rf_rhyper
#define rlnorm          Rf_rlnorm
#define rlogis          Rf_rlogis
#define rmultinom       Rf_rmultinom
#define rnbeta          Rf_rnbeta
#define rnbinom         Rf_rnbinom
#define rnbinom_mu      Rf_rnbinom_mu
#define rnchisq         Rf_rnchisq
#define rnf             Rf_rnf
#define rnorm           Rf_rnorm
#define rnt             Rf_rnt
#define rpois           Rf_rpois
#define rsignrank       Rf_rsignrank
#define rt              Rf_rt
#define rtukey          Rf_rtukey
#define runif           Rf_runif
#define rweibull        Rf_rweibull
#define rwilcox         Rf_rwilcox

#define sign            Rf_sign
#define tetragamma      Rf_tetragamma
#define trigamma        Rf_trigamma

#endif  // MATHLIB_STANDALONE || R_NO_REMAP_RMATH

//─────────────────────────────────────────────────────────────────────────────
// 6) DROPPED items (cannot run on‐device)
//─────────────────────────────────────────────────────────────────────────────
//  * External function prototypes (e.g. R_pow, R_pow_di, distribution routines)
//  * Random number generators & set_seed/get_seed
//  * C++ extern "C" linkage blocks
//  * Inline wrappers for pnorm/dnorm etc.
//
// If you need any of these on‐device, you must port or inline their bodies
// directly into your OpenCL kernel source.
//
// End of Rmath.cl