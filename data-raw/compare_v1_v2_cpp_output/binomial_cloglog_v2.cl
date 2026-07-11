/*
 * OPENCL.cl - global OpenCL configuration header
 *
 * This file is intended to be stitched at the top of an assembled OpenCL
 * program. Packages can copy and customize this header for their own kernels.
 */

#ifndef OPENCLPORT_OPENCL_CL
#define OPENCLPORT_OPENCL_CL

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#if defined(__OPENCL_C_VERSION__) && (__OPENCL_C_VERSION__ >= 120)
  #define HAVE_EXPM1 1
  #define HAVE_LOG1P 1
  #define HAVE_WORKING_ISFINITE 1
#else
  #define HAVE_EXPM1 0
  #define HAVE_LOG1P 0
  #define HAVE_WORKING_ISFINITE 0
#endif

/*
 * OpenCL C has no portable long double math surface. Keep HAVE_LONG_DOUBLE
 * undefined so nmath sources stay on double-precision branches (log/log1p/exp/fabs)
 * instead of long-double branches (logl/log1pl/expl/fabsl).
 */
#ifdef HAVE_LONG_DOUBLE
  #undef HAVE_LONG_DOUBLE
#endif

#ifndef ML_NAN
  #define ML_NAN (0.0/0.0)
  #define ML_POSINF INFINITY
  #define ML_NEGINF -INFINITY
#endif

#ifndef INLINE
  #define INLINE static inline
#endif

#ifndef R_UNUSED
  #define R_UNUSED(x) (void)(x)
#endif

#endif

// @source_type: shim
// @source_origin: libR
// @provides: R_pow, R_pow_di, R_CheckStack

/*
 * Minimal device-side libR runtime shim for OpenCL.
 * Keep this focused to symbols required by active kernel call paths.
 */

#ifndef OPENCLPORT_LIBR_SHIMS_LIBR_CL
#define OPENCLPORT_LIBR_SHIMS_LIBR_CL

INLINE double R_pow(double x, double y) {
    return pow(x, y);
}

INLINE double R_pow_di(double x, int n) {
    double p = 1.0;

    if (isnan(x)) return x;
    if (n != 0) {
        if (!isfinite(x)) return R_pow(x, (double)n);
        if (n < 0) { n = -n; x = 1.0 / x; }
        for (;;) {
            if (n & 1) p *= x;
            n >>= 1;
            if (!n) break;
            x *= x;
        }
    }
    return p;
}

INLINE void R_CheckStack(void) {
    /*
     * No-op on device: host stack checks are not meaningful for OpenCL kernels.
     * This satisfies translated R math paths that call R_CheckStack().
     */
}

#endif


// @source_type: h
// @source_origin: Boolean.h
// @includes: R_ext.h, stdbool.h, Rconfig.h
// @depends:
// @provides: R_EXT_BOOLEAN_H_, Rboolean, int, FALSE, TRUE
// @used: Rboolean, int, FALSE, TRUE
// @used_includes: stdbool.h, Rconfig.h
// @to_shim: Rboolean, int, FALSE, TRUE
// @to_shim_deterministic: FALSE, TRUE
// @to_shim_reason: Rboolean, int
// @to_shim_kind: Rboolean=reason, int=reason, FALSE=literal_constant, TRUE=literal_constant
// @all_depends_count: 0
// @load_order: 5

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000, 2026 The R Core Team.
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Included by R.h: API */

#ifndef R_EXT_BOOLEAN_H_
#define R_EXT_BOOLEAN_H_
// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
#if !defined(R_INCLUDE_BOOLEAN_H) || R_INCLUDE_BOOLEAN_H

// NB: there is a version of this in Rmath.h0[.in]

#undef FALSE
#undef TRUE

/* Ensuure a 'bool' type is available.  We could use
   __bool_true_false_are_defined, 
   but that was declared obsolescent in C23.
*/
#if defined __STDC_VERSION__ && __STDC_VERSION__ > 202000L
// C23 so bool is a keyword
#elif defined __cplusplus
// part of C++ >= 11, which is all R supports.
#else
// openclport-disabled-include: # include <stdbool.h>
// stdbool.h is C99, so available everywhere.
#endif

// openclport-disabled-include: #include <Rconfig.h> /* for HAVE_ENUM_BASE_TYPE */
/*
  Setting the underlying aka base type is supported in C23, C++11 
  and some C compilers based on clang and some versions of GCC.
  What matters here is the C compiler used to build R.
 */
#ifdef  __cplusplus
extern "C" {
#endif
#ifdef HAVE_ENUM_BASE_TYPE
// Apple clang 17 warns even in C23 mode: gcc warns about #pragma clang
// LLVM clang no longer warns.
// Apple clang 21 no longer has that warning.
# if defined  __APPLE__ && defined __clang__ && __clang_major__ < 21
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wfixed-enum-extension"
# endif

  typedef enum :int { FALSE = 0, TRUE } Rboolean;  // so NOT NA

# if defined  __APPLE__ && defined __clang__ && __clang_major__ < 21
#  pragma clang diagnostic pop
# endif
#else
    typedef enum { FALSE = 0, TRUE } Rboolean;  // so NOT NA
#endif
#ifdef  __cplusplus
}
#endif

#else
/* The Rbolean type is used in too many R headers to condition them
 * all.  However, people defining R_INCLUDE_BOOLEAN_H=0 should not be
 * using it in their own code, and its base type is expected to be int
 * (and guaranteed to be on most platforms as from R 4.5.0). */

    typedef Rboolean int;
#endif /* R_INCLUDE_BOOLEAN_H = 0 */
#endif /* R_EXT_BOOLEAN_H_ */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: Complex.h
// @includes: R_ext.h
// @depends:
// @provides: R_COMPLEX_H
// @all_depends_count: 0
// @load_order: 6

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2023   The R Core Team
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Included by R.h: Part of the API. */

#ifndef R_COMPLEX_H
#define R_COMPLEX_H

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
#ifdef  __cplusplus
extern "C" {
#endif

# ifdef R_LEGACY_RCOMPLEX

/* This definition does not work with optimizing compilers which take
advantage of strict aliasing rules.  It is not safe to use with Fortran
COMPLEX*16 (PR#18430) or in arguments to library calls expecting C99
_Complex double.  This definition should not be used, but if it were still
necessary, one should at least disable LTO.
*/

typedef struct {
 	double r;
 	double i;
 } Rcomplex;

# else

/* This definition uses an anonymous structure, which is defined in C11 (but
not C99).  It is, however, supported at least by GCC, clang and icc.  The
private_data_c member should never be used in code, but tells the compiler
about type punning when accessing the .r and .i elements, so is safer to use
when interfacing with Fortran COMPLEX*16 or directly C99 _Complex double
(PR#18430).

This form of static initialization works with both definitions:
Rcomplex z = { .r = 1, .i = 2 };

Anonymous structures and C99 _Complex have not been incorporated into C++
standard.  While they are usually supported as compiler extensions, warnings
are typically issued (-pedantic) by a C++ compiler.
*/

#ifdef __cplusplus
// Look for clang first as it defines __GNUC__ and reacts to #pragma GCC
# if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wgnu-anonymous-struct"
#  pragma clang diagnostic ignored "-Wc99-extensions"
# elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wpedantic"
# endif
#endif

typedef union {
    struct {
	double r;
	double i;
    };
    double _Complex private_data_c;
} Rcomplex;

#ifdef __cplusplus
# if defined(__clang__)
#  pragma clang diagnostic pop
# elif defined(__GNUC__)
#  pragma GCC diagnostic pop
# endif
#endif

# endif 

#ifdef  __cplusplus
}
#endif

#endif /* R_COMPLEX_H */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: Constants.h
// @includes: R_ext.h, cfloat, float.h
// @depends:
// @provides: R_EXT_CONSTANTS_H_, M_PI
// @used: M_PI
// @used_includes: float.h
// @builtin: M_PI
// @all_depends_count: 0
// @load_order: 7

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998-2024   The R Core Team.
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Included by R.h: Part of the API. */

#ifndef R_EXT_CONSTANTS_H_
#define R_EXT_CONSTANTS_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
/* usually in math.h, but not with strict C99/C++11 compliance.
   Also in Rmath.h
 */
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375
#endif

/*
  S-compatibility defines.
 */
#ifdef __cplusplus
// openclport-disabled-include: # include <cfloat>   /* Defines the RHSs, C++11 and later */
#else
// openclport-disabled-include: # include <float.h>  /* Defines the RHSs, C99 and later */
#endif

/* #ifndef STRICT_R_HEADERS
# define PI             M_PI
#endif
*/

/* The DOUBLE_* defines were deprecated in R 4.2.0 and removed in 4.3.0.
#define DOUBLE_DIGITS  DBL_MANT_DIG
#define DOUBLE_EPS     DBL_EPSILON
#define DOUBLE_XMAX    DBL_MAX
#define DOUBLE_XMIN    DBL_MIN
*/

#endif /* R_EXT_CONSTANTS_H_ */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: Visibility.h
// @includes: R_ext.h, Rconfig.h
// @depends:
// @provides: R_EXT_VISIBILITY_H_, attribute_visible, attribute_hidden
// @used: attribute_hidden
// @used_includes: Rconfig.h
// @to_shim: attribute_hidden
// @to_shim_reason: attribute_hidden
// @to_shim_kind: attribute_hidden=reason
// @all_depends_count: 0
// @load_order: 4

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2008    the R Core Team
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/*
  Definitions controlling visibility on some platforms.

  Part of the API.
*/

#ifndef R_EXT_VISIBILITY_H_
#define R_EXT_VISIBILITY_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
// openclport-disabled-include: #include <Rconfig.h>

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_visible __attribute__ ((visibility ("default")))
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_visible
# define attribute_hidden
#endif

#endif /* R_EXT_VISIBILITY_H_ */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: libextern.h
// @includes: R_ext.h
// @depends:
// @provides: LibImport, LibExport, LibExtern, extern
// @used: extern
// @to_shim: extern
// @to_shim_reason: extern
// @to_shim_kind: extern=reason
// @all_depends_count: 0
// @load_order: 8

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001, 2022  The R Core Team.
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Included by R.h: API on Windows */

/* don't disallow including this one more than once */

/* This is intended to be called from other header files, so not callable
   from C++ */

#undef LibExtern
#undef LibImport
#undef LibExport

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
#ifdef _WIN32 /* _WIN32 as does not depend on config.h */
#define LibImport __declspec(dllimport)
/* exporting is now done via .def file in R */
#define LibExport /* __declspec(dllexport) */
#else
#define LibImport
#define LibExport
#endif

#ifdef __MAIN__
#define LibExtern LibExport
#define extern
#elif defined(R_DLL_BUILD)
#define LibExtern extern
#else
#define LibExtern extern LibImport
#endif

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----


// @source_type: h
// @source_origin: Rconfig.h
// @provides: R_CONFIG_H, SIZEOF_SIZE_T, R_INLINE, IEEE_754

#ifndef R_CONFIG_H
#define R_CONFIG_H

/*
 * Minimal shim config for OpenCL portability scaffolding.
 * Keep this intentionally tiny: just enough for dependent typedef logic.
 */
#ifndef SIZEOF_SIZE_T
#define SIZEOF_SIZE_T 8
#endif

#ifndef R_INLINE
#define R_INLINE inline
#endif

#ifndef IEEE_754
#define IEEE_754 1
#endif

#endif /* R_CONFIG_H */


// @source_type: h
// @source_origin: Rinternals.h
// @includes: Rconfig.h, Boolean.h, Complex.h
// @depends: Rconfig
// @provides: R_INTERNALS_H_, Rbyte, R_len_t, R_xlen_t, SEXP, PROTECT_INDEX

#ifndef R_INTERNALS_H_
#define R_INTERNALS_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <Rconfig.h>
// openclport-disabled-include: #include <R_ext/Boolean.h>
// openclport-disabled-include: #include <R_ext/Complex.h>

/*
 * Minimal R internals shim surface for parsing selected R_ext headers.
 * Opaque SEXP plus length/type aliases only.
 */
typedef unsigned char Rbyte;
typedef int R_len_t;

#if (SIZEOF_SIZE_T > 4)
typedef long long R_xlen_t;
#else
typedef int R_xlen_t;
#endif

typedef struct SEXPREC *SEXP;
typedef int PROTECT_INDEX;

#endif /* R_INTERNALS_H_ */

// @source_type: h
// @source_origin: Rdefines.h
// @includes: Rinternals.h
// @depends: Rinternals
// @provides: R_DEFINES_H_, NA_STRING

#ifndef R_DEFINES_H_
#define R_DEFINES_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <Rinternals.h>

/* Minimal compatibility aliases used by some downstream headers. */
#ifndef NA_STRING
#define NA_STRING ((SEXP)0)
#endif

#endif /* R_DEFINES_H_ */


// @source_type: h
// @source_origin: Arith.h
// @includes: R_ext.h, libextern.h, math.h
// @depends:
// @provides: R_ARITH_H_, NA_LOGICAL, NA_INTEGER, NA_REAL, ISNA, ISNAN, R_FINITE, R_IsNA		, R_IsNaN		, R_finite		, R_isnancpp 
// @used: NA_INTEGER, NA_REAL, ISNAN, R_FINITE, R_finite, R_isnancpp
// @used_includes: libextern.h, math.h, R_ext/libextern.h
// @to_shim: NA_INTEGER, NA_REAL, ISNAN, R_FINITE, R_finite, R_isnancpp
// @to_shim_deterministic: NA_INTEGER, NA_REAL, ISNAN, R_FINITE
// @to_shim_reason: R_finite, R_isnancpp
// @to_shim_kind: NA_INTEGER=define_object_identifier, NA_REAL=define_object_identifier, ISNAN=define_function_macro, R_FINITE=define_function_macro, R_finite=reason, R_isnancpp=reason
// @all_depends_count: 0
// @load_order: 16

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998--2016  The R Core Team.
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Included by R.h: Part of the API. */

#ifndef R_ARITH_H_
#define R_ARITH_H_

/* 
   This used to define _BSD_SOURCE to make declarations of isfinite
   and isnan visible in glibc.  But that was deprecated in glibc 2.20,
   and --std=c99 suffices nowadays.
*/

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
// openclport-disabled-include: #include <R_ext/libextern.h>
#ifdef  __cplusplus
extern "C" {
#else
/* needed for isnan and isfinite, neither of which are used under C++ */
// openclport-disabled-include: # include <math.h>
#endif

/* implementation of these : ../../main/arithmetic.c */
LibExtern double R_NaN;		/* IEEE NaN */
LibExtern double R_PosInf;	/* IEEE Inf */
LibExtern double R_NegInf;	/* IEEE -Inf */
LibExtern double R_NaReal;	/* NA_REAL: IEEE */
LibExtern int	 R_NaInt;	/* NA_INTEGER:= INT_MIN currently */
#ifdef __MAIN__
#undef extern
#undef LibExtern
#endif

#define NA_LOGICAL	R_NaInt
#define NA_INTEGER	R_NaInt
/* #define NA_FACTOR	R_NaInt  unused */
#define NA_REAL		R_NaReal
/* NA_STRING is a SEXP, so defined in Rinternals.h */

int R_IsNA(double);		/* True for R's NA only */
int R_IsNaN(double);		/* True for special NaN, *not* for NA */
int R_finite(double);		/* True if none of NA, NaN, +/-Inf */
#define ISNA(x)	       R_IsNA(x)

/* ISNAN(): True for *both* NA and NaN.
   NOTE: some systems do not return 1 for TRUE.
   Also note that C++ math headers specifically undefine
   isnan if it is a macro (it is on macOS and in C99),
   hence the workaround.  This code also appears in Rmath.h
*/
#ifdef __cplusplus
  int R_isnancpp(double); /* in arithmetic.c */
#  define ISNAN(x)     R_isnancpp(x)
#else
#  define ISNAN(x)     (isnan(x)!=0)
#endif

/* Configure-time feature gate:
   host R builds may define this via config headers; OpenCL/device builds
   should provide it from prelude/shim configuration. */
#ifdef HAVE_WORKING_ISFINITE
/* isfinite is defined in <math.h> according to C99 */
# define R_FINITE(x)    isfinite(x)
#else
# define R_FINITE(x)    R_finite(x)
#endif

#ifdef  __cplusplus
}
#endif

#endif /* R_ARITH_H_ */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: MathThreads.h
// @includes: R_ext.h, libextern.h
// @depends:
// @provides: R_EXT_MATHTHREADS_H_
// @used_includes: libextern.h, R_ext/libextern.h
// @all_depends_count: 0
// @load_order: 9

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000-2026 The R Core Team.
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/*
  Experimental: included by src/library/stats/src/distance.c

  This is not used currently on Windows.
*/

#ifndef R_EXT_MATHTHREADS_H_
#define R_EXT_MATHTHREADS_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
// openclport-disabled-include: #include <R_ext/libextern.h>

#ifdef  __cplusplus
extern "C" {
#endif

#ifdef USE_MATH_THREADS
LibExtern int R_num_math_threads;
LibExtern int R_max_num_math_threads;
#endif

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_MATHTHREADS_H_ */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: Memory.h
// @includes: R_ext.h, cstddef, stddef.h
// @depends:
// @provides: R_EXT_MEMORY_H_, R_SIZE_T, vmaxget , vmaxset , R_gc, R_gc_running, R_alloc, S_alloc, S_realloc, R_malloc_gc, R_calloc_gc, R_realloc_gc
// @used: vmaxget, vmaxset, R_alloc
// @used_includes: stddef.h
// @to_shim: vmaxget, vmaxset, R_alloc
// @to_shim_reason: vmaxget, vmaxset, R_alloc
// @to_shim_kind: vmaxget=reason, vmaxset=reason, R_alloc=reason
// @all_depends_count: 0
// @load_order: 10

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2024  The R Core Team
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 * Memory Allocation (garbage collected) --- INCLUDING S compatibility ---
 */

/* Included by R.h: Part of the API. */

#ifndef R_EXT_MEMORY_H_
#define R_EXT_MEMORY_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
#if defined(__cplusplus) && !defined(DO_NOT_USE_CXX_HEADERS)
// openclport-disabled-include: # include <cstddef>
# define R_SIZE_T std::size_t
#else
// openclport-disabled-include: # include <stddef.h> /* for size_t */
# define R_SIZE_T size_t
#endif

#ifdef  __cplusplus
extern "C" {
#endif

void*	vmaxget(void); // not remapped
void	vmaxset(const void *); // not re-mapped

void	R_gc(void);
#ifdef USE_BASE_R_SUPPORT
int	R_gc_running(void);
#endif

char*	R_alloc(R_SIZE_T, int);
long double *R_allocLD(R_SIZE_T nelem);
char*	S_alloc(long, int);
char*	S_realloc(char *, long, long, int);

void *  R_malloc_gc(size_t);
void *  R_calloc_gc(size_t, size_t);
void *  R_realloc_gc(void *, size_t);

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_MEMORY_H_ */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
void* vmaxget(void) {
    return (void*)0;
}

void vmaxset(const void* v) {
    (void)v;
}

char* R_alloc(R_SIZE_T n, int size) {
    (void)n;
    (void)size;
    return (char*)0;
}
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: Print.h
// @includes: R_ext.h, cstdarg, stdarg.h
// @depends:
// @provides: R_EXT_PRINT_H_, R_USE_C99_IN_CXX, R_VA_LIST, R_PRINTF_FORMAT, Rprintf, REprintf, Rvprintf, REvprintf
// @used: REprintf
// @used_includes: cstdarg, stdarg.h
// @to_shim: REprintf
// @to_shim_reason: REprintf
// @to_shim_kind: REprintf=reason
// @all_depends_count: 0
// @load_order: 12

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2024    The R Core Team
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Included by R.h: Part of the API. */

#ifndef R_EXT_PRINT_H_
#define R_EXT_PRINT_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
#ifdef  __cplusplus
/* If the vprintf interface is defined at all in C++ it may only be
   defined in namespace std.  It is part of the C++11 standard. */
# if __cplusplus >= 201103L && !defined(R_USE_C99_IN_CXX)
#  define R_USE_C99_IN_CXX
# endif
# ifdef R_USE_C99_IN_CXX
// openclport-disabled-include: #  include <cstdarg>
#  define R_VA_LIST std::va_list
# endif
extern "C" {
#else
// openclport-disabled-include: # include <stdarg.h>
# define R_VA_LIST va_list
#endif

#ifdef __GNUC__
# ifdef _WIN32
#  if defined(_UCRT) || ((__MSVCRT_VERSION__ >= 0x1400) || \
                        (__MSVCRT_VERSION__ >= 0xE00 && __MSVCRT_VERSION__ < 0x1000))
#   if defined(__clang__)
#    define R_PRINTF_FORMAT(M,N) __attribute__ ((format (printf, M, N)))    
#   else
#    define R_PRINTF_FORMAT(M,N) __attribute__ ((format (gnu_printf, M, N)))    
#   endif
#  else
#   define R_PRINTF_FORMAT(M,N)
#  endif
# else
#  define R_PRINTF_FORMAT(M,N) __attribute__ ((format (printf, M, N)))
# endif
#else
# define R_PRINTF_FORMAT(M,N)
#endif

#ifdef __OPENCL_VERSION__
/* OpenCL C disallows variadic function declarations/prototypes. */
# undef Rprintf
# undef REprintf
# define Rprintf(...) ((void)0)
# define REprintf(...) ((void)0)

# if !defined(__cplusplus) || defined R_USE_C99_IN_CXX
#  undef Rvprintf
#  undef REvprintf
#  define Rvprintf(...) ((void)0)
#  define REvprintf(...) ((void)0)
# endif
#else
void Rprintf(const char *, ...) R_PRINTF_FORMAT(1, 2);
void REprintf(const char *, ...) R_PRINTF_FORMAT(1, 2);

# if !defined(__cplusplus) || defined R_USE_C99_IN_CXX

void Rvprintf(const char *, R_VA_LIST) R_PRINTF_FORMAT(1, 0);
void REvprintf(const char *, R_VA_LIST) R_PRINTF_FORMAT(1, 0);

# endif
#endif

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_PRINT_H_ */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: RS.h
// @includes: R_ext.h, cstring, cstddef, string.h, stddef.h, Rconfig.h
// @depends:
// @provides: R_RS_H, R_SIZE_T, R_Calloc, R_Realloc, R_Free, Memcpy, Memzero, CallocCharBuf, F77_CALL, F77_NAME, F77_SUB, R_chk_free
// @used: F77_SUB
// @used_includes: string.h, stddef.h, Rconfig.h
// @to_shim: F77_SUB
// @to_shim_deterministic: F77_SUB
// @to_shim_kind: F77_SUB=define_function_macro
// @all_depends_count: 0
// @load_order: 3

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999-2025 The R Core Team.
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 * 
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Included by R.h: nowadays almost all API */

#ifndef R_RS_H
#define R_RS_H

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
#if defined(__cplusplus) && !defined(DO_NOT_USE_CXX_HEADERS)
// openclport-disabled-include: # include <cstring>
// openclport-disabled-include: # include <cstddef>
# define R_SIZE_T std::size_t
#else
// openclport-disabled-include: # include <string.h>		/* for memcpy, memset */
// openclport-disabled-include: # include <stddef.h> /* for size_t */
# define R_SIZE_T size_t
#endif

// openclport-disabled-include: #include <Rconfig.h>		/* for HAVE_F77_UNDERSCORE */

#ifdef  __cplusplus
extern "C" {
#endif

/* S Like Memory Management */

/* not of themselves API */
extern void *R_chk_calloc(R_SIZE_T, R_SIZE_T);
extern void *R_chk_realloc(void *, R_SIZE_T);
extern void R_chk_free(void *);
extern void *R_chk_memcpy(void *, const void *, R_SIZE_T);
extern void *R_chk_memset(void *, int, R_SIZE_T);

/* API */
#define R_Calloc(n, t)   (t *) R_chk_calloc( (R_SIZE_T) (n), sizeof(t) )
#define R_Realloc(p,n,t) (t *) R_chk_realloc( (void *)(p), (R_SIZE_T)((n) * sizeof(t)) )
#define R_Free(p)      (R_chk_free( (void *)(p) ), (p) = NULL)

/* Nowadays API: undocumented until 4.1.2: widely used. */
#define Memcpy(p,q,n)  R_chk_memcpy( p, q, (R_SIZE_T)(n) * sizeof(*p) )

/* Nowadays API: added for 3.0.0 but undocumented until 4.1.2. */
#define Memzero(p,n)  R_chk_memset(p, 0, (R_SIZE_T)(n) * sizeof(*p))

/* API: Added in R 2.6.0 */
#define CallocCharBuf(n) (char *) R_chk_calloc(((R_SIZE_T)(n))+1, sizeof(char))

/* S Like Fortran Interface */
/* These may not be adequate everywhere. Convex had _ prepending common
   blocks, and some compilers may need to specify Fortran linkage.

   HP-UX did not add a trailing underscore.  (It still existed in
   2024, but R ports had not been seen for many years.)

   Note that this is an F77 interface, intended only for valid F77
   names of <= 6 ASCII characters (and no underscores) and there is an
   implicit assumption that the Fortran compiler maps names to
   lower-case (and 'x' is lower-case when called).

   The configure code has

   HAVE_F77_EXTRA_UNDERSCORE
   Define if your Fortran compiler appends an extra_underscore to
   external names containing an underscore.

   but that is not used here (and none of gfortran, flang-new nor
   x86_64 ifx do so: earlier Intel x86 compilere might have).  It is
   used in Rdynload.c to support .Fortran.

   These macros have always been the same in R.  Their documented uses are

   F77_SUB to define a function in C to be called from Fortran 
   F77_NAME to declare a Fortran routine in C before use 
   F77_CALL to call a Fortran routine from C
 */

#ifdef HAVE_F77_UNDERSCORE
# define F77_CALL(x)	x ## _
#else
# define F77_CALL(x)	x
#endif
#define F77_NAME(x)    F77_CALL(x)
#define F77_SUB(x)     F77_CALL(x)
/* Last two were historical from S, not used in R, deprecated in 4.4.2, removed in 4.5.0
#define F77_COM(x)     F77_CALL(x)
#define F77_COMDECL(x) F77_CALL(x)
*/
 
/* call_R was deprecated in R 2.15.0, removed in R 4.2.0 */

#ifdef  __cplusplus
}
#endif

#endif /* R_RS_H */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: Random.h
// @includes: R_ext.h, Boolean.h
// @depends:
// @provides: R_RANDOM_H, Int32, WICHMANN_HILL, MARSAGLIA_MULTICARRY, SUPER_DUPER, MERSENNE_TWISTER, KNUTH_TAOCP, USER_UNIF, KNUTH_TAOCP2, LECUYER_CMRG, BUGGY_KINDERMAN_RAMAGE, AHRENS_DIETER, BOX_MULLER, USER_NORM, INVERSION, KINDERMAN_RAMAGE, ROUNDING, REJECTION, R_sample_kind, GetRNGstate, PutRNGstate, unif_rand, R_unif_index, norm_rand, exp_rand, user_unif_rand, user_unif_init, user_unif_nseed, user_unif_seedloc, user_norm_rand
// @used: BUGGY_KINDERMAN_RAMAGE, AHRENS_DIETER, BOX_MULLER, USER_NORM, INVERSION, KINDERMAN_RAMAGE, unif_rand, R_unif_index, norm_rand, exp_rand
// @used_includes: Boolean.h, R_ext/Boolean.h
// @to_shim: BUGGY_KINDERMAN_RAMAGE, AHRENS_DIETER, BOX_MULLER, USER_NORM, INVERSION, KINDERMAN_RAMAGE, unif_rand, R_unif_index, norm_rand, exp_rand
// @to_shim_reason: BUGGY_KINDERMAN_RAMAGE, AHRENS_DIETER, BOX_MULLER, USER_NORM, INVERSION, KINDERMAN_RAMAGE, unif_rand, R_unif_index, norm_rand, exp_rand
// @to_shim_kind: BUGGY_KINDERMAN_RAMAGE=reason, AHRENS_DIETER=reason, BOX_MULLER=reason, USER_NORM=reason, INVERSION=reason, KINDERMAN_RAMAGE=reason, unif_rand=reason, R_unif_index=reason, norm_rand=reason, exp_rand=reason
// @all_depends_count: 0
// @load_order: 13

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2022    The R Core Team
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Included by R.h: Part of the API. */

#ifndef R_RANDOM_H
#define R_RANDOM_H

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
// openclport-disabled-include: #include <R_ext/Boolean.h>

#ifdef  __cplusplus
extern "C" {
#endif

typedef enum {
    WICHMANN_HILL,
    MARSAGLIA_MULTICARRY,
    SUPER_DUPER,
    MERSENNE_TWISTER,
    KNUTH_TAOCP,
    USER_UNIF,
    KNUTH_TAOCP2,
    LECUYER_CMRG
} RNGtype;

/* Different kinds of "N(0,1)" generators :*/
typedef enum {
    BUGGY_KINDERMAN_RAMAGE,
    AHRENS_DIETER,
    BOX_MULLER,
    USER_NORM,
    INVERSION,
    KINDERMAN_RAMAGE
} N01type;

/* Different ways to generate discrete uniform samples */
typedef enum {
    ROUNDING,
    REJECTION
} Sampletype;
Sampletype R_sample_kind(void);

void GetRNGstate(void);
void PutRNGstate(void);

double unif_rand(void);
double R_unif_index(double);
/* These are also defined in Rmath.h */
double norm_rand(void);
double exp_rand(void);

typedef unsigned int Int32;
double * user_unif_rand(void);
void user_unif_init(Int32);
int * user_unif_nseed(void);
int * user_unif_seedloc(void);

double * user_norm_rand(void);

#ifdef  __cplusplus
}
#endif

#endif /* R_RANDOM_H */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: Riconv.h
// @includes: R_ext.h
// @depends:
// @provides: R_ICONV_H, Riconv_open, Riconv_close
// @all_depends_count: 0
// @load_order: 2

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2005     the R Core Team
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/*
  Interface to R's platform-independent implementation of iconv.

  Part of the API.
*/

#ifndef R_ICONV_H
#define R_ICONV_H

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
#ifdef  __cplusplus
extern "C" {
#endif

/* from sysutils.c */
#undef Riconv_open
#undef Riconv
#undef Riconv_close
void * Riconv_open (const char* tocode, const char* fromcode);
size_t Riconv (void * cd, const char **inbuf, size_t *inbytesleft,
	       char  **outbuf, size_t *outbytesleft);
int Riconv_close (void * cd);

#ifdef  __cplusplus
}
#endif

#endif /* R_ICONV_H */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: Utils.h
// @includes: R_ext.h, Boolean.h, Complex.h, cstddef, stddef.h
// @depends:
// @provides: R_EXT_UTILS_H_, R_SIZE_T, revsort, iPsort, rPsort, cPsort, IndexWidth, StringFalse, StringTrue, isBlankString, R_isort, R_rsort, R_csort, rsort_with_index , R_qsort, R_qsort_I, R_qsort_int, R_qsort_int_I, F77_NAME, F77_SUB, StringFalse , StringTrue , isBlankString , R_atof, R_strtod, R_free_tmpnam, R_CheckUserInterrupt, R_CheckStack, R_CheckStack2, R_max_col
// @used: F77_SUB, R_CheckUserInterrupt, R_CheckStack
// @used_includes: Boolean.h, stddef.h, R_ext/Boolean.h
// @to_shim: F77_SUB, R_CheckUserInterrupt, R_CheckStack
// @to_shim_reason: F77_SUB, R_CheckUserInterrupt, R_CheckStack
// @to_shim_kind: F77_SUB=reason, R_CheckUserInterrupt=reason, R_CheckStack=reason
// @all_depends_count: 0
// @load_order: 15

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2026    The R Core Team
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.

 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 * Generally useful  UTILITIES  *NOT* relying on R internals (from Defn.h)
 */

/* Included by R.h: some are API (documented in R-exts), 
   others are noted below. */

#ifndef R_EXT_UTILS_H_
#define R_EXT_UTILS_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
// openclport-disabled-include: #include <R_ext/Boolean.h>
// openclport-disabled-include: #include <R_ext/Complex.h>

#if defined(__cplusplus) && !defined(DO_NOT_USE_CXX_HEADERS)
// openclport-disabled-include: # include <cstddef>
# define R_SIZE_T std::size_t
#else
// openclport-disabled-include: # include <stddef.h>
# define R_SIZE_T size_t
#endif

#define revsort       Rf_revsort
#define iPsort        Rf_iPsort
#define rPsort        Rf_rPsort
#define cPsort        Rf_cPsort
#define IndexWidth    Rf_IndexWidth
//#define setIVector    Rf_setIVector
//#define setRVector    Rf_setRVector
#define StringFalse   Rf_StringFalse
#define StringTrue    Rf_StringTrue
#define isBlankString Rf_isBlankString

#ifdef  __cplusplus
extern "C" {
#endif

/* ../../main/sort.c : */
void	R_isort(int*, int);
void	R_rsort(double*, int);
void	R_csort(Rcomplex*, int);
void    rsort_with_index(double *, int *, int); // not remapped.
void	revsort(double*, int*, int);/* reverse; sort i[] alongside */
void	iPsort(int*,    int, int);
void	rPsort(double*, int, int);
void	cPsort(Rcomplex*, int, int);

/* ../../main/qsort.c : */
/* dummy renamed to II to avoid problems with g++ on Solaris */
void R_qsort    (double *v,         R_SIZE_T i, R_SIZE_T j);
void R_qsort_I  (double *v, int *II, int i, int j);
void R_qsort_int  (int *iv,         R_SIZE_T i, R_SIZE_T j);
void R_qsort_int_I(int *iv, int *II, int i, int j);
#ifdef R_RS_H
void F77_NAME(qsort4)(double *v, int *indx, int *ii, int *jj);
void F77_NAME(qsort3)(double *v,            int *ii, int *jj);
#endif

// listed as callable from C in WRE
#ifdef R_RS_H
int F77_SUB(i1mach)(int *i);
double F77_SUB(d1mach)(int *i);
#endif
    
/* ../../main/util.c  and others : */
const char *R_ExpandFileName(const char *);
#ifdef Win32
// not API
const char *R_ExpandFileNameUTF8(const char *);
#endif
/*  attribute_hidden and no longer used.
void	setIVector(int*, int, int);
void	setRVector(double*, int, double);
*/
/* Not API */
Rboolean StringFalse(const char *); // used by iotools
Rboolean StringTrue(const char *); // used by iotools
Rboolean isBlankString(const char *); // used by iotools and openxlsx2

/* These two are guaranteed to use '.' as the decimal point,
   and to accept "NA". Documented since 4.4.0 patched.
 */
double R_atof(const char *str);
double R_strtod(const char *c, char **end);

char *R_tmpnam(const char *prefix, const char *tempdir);
char *R_tmpnam2(const char *prefix, const char *tempdir, const char *fileext);
void R_free_tmpnam(char *name);

void R_CheckUserInterrupt(void);
void R_CheckStack(void);
void R_CheckStack2(R_SIZE_T);


/* ../../appl/interv.c: first and also in Applic.h 
   Both are API
*/
int findInterval(double *xt, int n, double x,
		 Rboolean rightmost_closed,  Rboolean all_inside, int ilo,
		 int *mflag);
int findInterval2(double *xt, int n, double x,
		  Rboolean rightmost_closed,  Rboolean all_inside, Rboolean left_open,
		  int ilo, int *mflag);
/* Removed in 4.5.0, but still API according too WRE */
#ifdef R_RS_H
// Was Rboolean*, but that is not possible in Fortran.
int F77_SUB(interv)(double *xt, int *n, double *x,
		    int *rightmost_closed, int *all_inside,
		    int *ilo, int *mflag);
#endif
 
/* not API, no longer in R
void find_interv_vec(double *xt, int *n,	double *x,   int *nx,
		     int *rightmost_closed, int *all_inside, int *indx);
*/

/* ../../appl/maxcol.c */
void R_max_col(double *matrix, int *nr, int *nc, int *maxes, int *ties_meth);

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_UTILS_H_ */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: Error.h
// @includes: R_ext.h, Print.h
// @depends: Print
// @provides: R_ERROR_H_, NORET, Rf_error, Rf_warning, UNIMPLEMENTED, WrongArgCount, error, warning, R_ShowMessage
// @used: error, warning
// @used_includes: Print.h, R_ext/Print.h
// @to_shim: error, warning
// @to_shim_deterministic: error, warning
// @to_shim_kind: error=define_object_identifier, warning=define_object_identifier
// @all_depends_count: 1
// @all_depends: Print
// @load_order: 17

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2025   The R Core Team
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Included by R.h: Part of the API. */

#ifndef R_ERROR_H_
#define R_ERROR_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
// openclport-disabled-include: #include <R_ext/Print.h> // for R_PRINTF_FORMAT

#ifdef  __cplusplus
extern "C" {
#endif

/*
 * As this is sometimes an attribute, it should precede 'static' in a
 * function declaration.
 * gcc 15 requires it to precede 'attribute_hidden'.
 * OTOH, '_Noreturn' is an obsolescent (in C23) function specifier.
 */
#if defined NORET
#elif (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 202301L)
// gcc 15 LLVM clang 19- and Apple clang 17
# define NORET [[noreturn]]
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201102L
# define NORET _Noreturn
#elif defined(__GNUC__) && __GNUC__ >= 3
// All platforms these days should be using C >= 11 but perhaps used for C++
# define NORET __attribute__((noreturn))
#else
// C++ and legacy
# define NORET
#endif

#ifdef __OPENCL_VERSION__
/* OpenCL C disallows variadic function declarations/prototypes. */
# undef Rf_error
# undef Rf_warning
# define Rf_error(...) ((void)0)
# define Rf_warning(...) ((void)0)
# define UNIMPLEMENTED(...) ((void)0)
# ifdef USE_BASE_R_SUPPORT
#  define WrongArgCount(...) ((void)0)
# endif
#else
#ifdef  __cplusplus
// Only supported in C++ >= 11, but that is all current R supports
// Defining NORET caused conflict in many C++-using packages
[[noreturn]] void Rf_error(const char *, ...) R_PRINTF_FORMAT(1, 2);

[[noreturn]] void UNIMPLEMENTED(const char *);
#ifdef USE_BASE_R_SUPPORT
[[noreturn]] void WrongArgCount(const char *);
#endif
#else
NORET void Rf_error(const char *, ...) R_PRINTF_FORMAT(1, 2);

NORET void UNIMPLEMENTED(const char *);
#ifdef USE_BASE_R_SUPPORT
NORET void WrongArgCount(const char *);
#endif
#endif

void Rf_warning(const char *, ...) R_PRINTF_FORMAT(1,2);
#endif

void R_ShowMessage(const char *s);

/* xerbla is a C function intended to be called from Fortran.
 * which forerly had a C declaration here.
 *
 * It wraps Rf_error, so use that directtly from C/C++
*/

#ifdef  __cplusplus
}
#endif

#ifndef R_NO_REMAP
#define error Rf_error
#define warning Rf_warning
#endif


#endif /* R_ERROR_H_ */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----


// @source_type: h
// @source_origin: R_ext.h
// @includes: Rconfig.h, Rinternals.h, Rdefines.h
// @provides: R_EXT_INTERNALS_AGGREGATE_H_
// @used_includes: Rconfig.h, Rinternals.h, Rdefines.h
// @all_depends_count: 0
// @load_order: 1

#ifndef R_EXT_INTERNALS_AGGREGATE_H_
#define R_EXT_INTERNALS_AGGREGATE_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <Rconfig.h>
// openclport-disabled-include: #include <Rinternals.h>
// openclport-disabled-include: #include <Rdefines.h>

#endif

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: Parse.h
// @includes: R_ext.h
// @depends: R_ext_internals
// @provides: R_EXT_PARSE_H_, PARSE_NULL, PARSE_OK, PARSE_INCOMPLETE, PARSE_ERROR, PARSE_EOF, R_ParseVector
// @all_depends_count: 1
// @all_depends: R_ext_internals
// @load_order: 11

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2006 R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* NOTE:
   This file exports a part of the current internal parse interface.
   It is subject to change at any minor (x.y.0) version of R.
   So not API.
 */

#ifndef R_EXT_PARSE_H_
#define R_EXT_PARSE_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
#ifdef __cplusplus
extern "C" {
#endif

/* PARSE_NULL will not be returned by R_ParseVector */
typedef enum {
    PARSE_NULL,
    PARSE_OK,
    PARSE_INCOMPLETE,
    PARSE_ERROR,
    PARSE_EOF
} ParseStatus;

SEXP R_ParseVector(SEXP, int, ParseStatus *, SEXP);

#ifdef __cplusplus
}
#endif

#endif

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----

// @source_type: h
// @source_origin: stats_package.h
// @includes: R_ext.h, Rconfig.h, Visibility.h
// @depends: R_ext_internals
// @provides: R_STATS_PACKAGE_H, NREG, OPT, F, F0, FDIF, G, HC, AI, AM, ALGSAV, COVMAT, COVPRT, COVREQ, DRADPR, DTYPE, IERR, INITH, INITS, IPIVOT, IVNEED, LASTIV, LASTV, LMAT, MXFCAL, MXITER, NEXTV, NFCALL, NFCOV, NFGCAL, NGCOV, NITER, NVDFLT, NVSAVE, OUTLEV, PARPRT, PARSAV, PERM, PRUNIT, QRTYP, RDREQ, RMAT, SOLPRT, STATPR, TOOBIG, VNEED, VSAVE, X0PRT
// @used: F, F0, G
// @used_includes: Rconfig.h, Visibility.h, R_ext/Visibility.h
// @to_shim: F, F0, G
// @to_shim_reason: F, F0, G
// @to_shim_kind: F=reason, F0=reason, G=reason
// @all_depends_count: 1
// @all_depends: R_ext_internals
// @load_order: 14

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2007--2025  The R Core Team.
 *
 *  This header file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This file is part of R. R is distributed under the terms of the
 *  GNU General Public License, either Version 2, June 1991 or Version 3,
 *  June 2007. See doc/COPYRIGHTS for details of the copyright status of R.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Not part of the API, not callable from C++ */

#ifndef R_STATS_PACKAGE_H
#define R_STATS_PACKAGE_H
// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <R_ext.h>
// openclport-disabled-include: #include <Rconfig.h>
// openclport-disabled-include: #include <R_ext/Visibility.h>

/*
#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif
*/

enum AlgType {NREG = 1, OPT = 2};
				/* 0-based indices into v */
enum  VPos {F = 9, F0 = 12, FDIF = 10, G = 27, HC = 70};
				/* 0-based indices into iv */
enum IVPos {AI = 90, AM = 94, ALGSAV = 50, COVMAT = 25,
	    COVPRT = 13, COVREQ = 14, DRADPR = 100,
	    DTYPE = 15, IERR = 74, INITH = 24, INITS = 24,
	    IPIVOT = 75, IVNEED =  2, LASTIV = 42, LASTV = 44,
	    LMAT =  41, MXFCAL = 16, MXITER = 17, NEXTV  = 46,
	    NFCALL =  5, NFCOV = 51, NFGCAL = 6, NGCOV = 52,
	    NITER = 30, NVDFLT = 49, NVSAVE = 8, OUTLEV = 18,
	    PARPRT = 19, PARSAV = 48, PERM = 57, PRUNIT = 20,
	    QRTYP = 79, RDREQ = 56, RMAT = 77, SOLPRT = 21,
	    STATPR = 22, TOOBIG = 1, VNEED = 3, VSAVE = 59,
	    X0PRT = 23};

attribute_hidden void
S_Rf_divset(int alg, int iv[], int liv, int lv, double v[]);

attribute_hidden void
S_nlsb_iterate(double b[], double d[], double dr[], int iv[],
	       int liv, int lv, int n, int nd, int p,
	       double r[], double rd[], double v[], double x[]);

attribute_hidden void
S_nlminb_iterate(double b[], double d[], double fx, double g[],
		 double h[], int iv[], int liv, int lv, int n,
		 double v[], double x[]);

attribute_hidden void
S_rcont2(int nrow, int ncol, const int nrowt[], const int ncolt[],
         int ntotal, const double fact[],
	 int jwork[], int matrix[]);

static R_INLINE int S_v_length(int alg, int n)
{
    return (alg - 1) ? (105 + (n * (2 * n + 20))) :
	(130 + (n * (n + 27))/2);
}

static R_INLINE int S_iv_length(int alg, int n)
{
    return (alg - 1) ? (82 + 4 * n) : (78 + 3 * n);
}

#endif /* R_STATS_PACKAGE_H */

// ---- BEGIN AUTO DETERMINISTIC SHIM ----
// (none)
// ---- END AUTO DETERMINISTIC SHIM ----

// ---- BEGIN MANUAL REASONED SHIM ----
// (none)
// ---- END MANUAL REASONED SHIM ----


// @source_type: h
// @source_origin: stdint_shim_opencl.h
// @provides: OPENCL_STDINT_SHIM_H, int64_t, uint64_t, int_least64_t, uint_least64_t
// @depends:
// @includes: System

#ifndef OPENCL_STDINT_SHIM_H
#define OPENCL_STDINT_SHIM_H

/* Minimal stdint shim for OpenCL nmath staging.
   Keep this intentionally small: add more typedefs only as needed. */

#if defined(__OPENCL_VERSION__) || defined(__OPENCL_C_VERSION__)
typedef long long int64_t;
typedef unsigned long long uint64_t;

typedef int64_t int_least64_t;
typedef uint64_t uint_least64_t;
#else
/* Host fallback if compiled outside OpenCL toolchain */
typedef long long int64_t;
typedef unsigned long long uint64_t;

typedef int64_t int_least64_t;
typedef uint64_t uint_least64_t;
#endif

#endif /* OPENCL_STDINT_SHIM_H */


// @source_type: h
// @source_origin: dpq.h
// @provides: give_log, R_D__0, R_D__1, R_D_Clog, R_D_Cval, R_D_exp, R_D_fexp, R_D_half, R_D_LExp, R_D_log, R_D_Lval, R_D_negInonint, R_D_nonint_check, R_D_qIv, R_D_rtxp, R_D_val, R_DT_0, R_DT_1, R_DT_Cexp, R_DT_CIv, R_DT_Clog, R_DT_Cval, R_DT_exp, R_DT_log, R_DT_Log, R_DT_qIv, R_DT_val, R_Log1_Exp, R_P_bounds_01, R_P_bounds_Inf_01, R_Q_P01_boundaries, R_Q_P01_check
// @all_depends_count: 0
// @load_order: 2

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000--2021 The  R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */
	/* Utilities for `dpq' handling (density/probability/quantile) */

/* give_log in "d";  log_p in "p" & "q" : */
#define give_log log_p
							/* "DEFAULT" */
							/* --------- */
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */
#define R_D_half (log_p ? -M_LN2 : 0.5)		// 1/2 (lower- or upper tail)


/* Use 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy */
#define R_D_Lval(p)	(lower_tail ? (p) : (0.5 - (p) + 0.5))	/*  p  */
#define R_D_Cval(p)	(lower_tail ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */

#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */
#define R_D_qIv(p)	(log_p	? exp(p) : (p))		/*  p  in qF(p,..) */
#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))	/* exp(x) */
#define R_D_log(p)	(log_p	?  (p)	 : log(p))	/* log(p) */
#define R_D_Clog(p)	(log_p	? log1p(-(p)) : (0.5 - (p) + 0.5)) /* [log](1-p) */

// log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x)) :
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))

/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : log1p(-x))

#define R_DT_val(x)	(lower_tail ? R_D_val(x)  : R_D_Clog(x))
#define R_DT_Cval(x)	(lower_tail ? R_D_Clog(x) : R_D_val(x))

/*#define R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		 *  p  in qF ! */
#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
			       : R_D_Lval(p))

/*#define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		 *  1 - p in qF */
#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) \
			       : R_D_Cval(p))

#define R_DT_exp(x)	R_D_exp(R_D_Lval(x))		/* exp(x) */
#define R_DT_Cexp(x)	R_D_exp(R_D_Cval(x))		/* exp(1 - x) */

#define R_DT_log(p)	(lower_tail? R_D_log(p) : R_D_LExp(p))/* log(p) in qF */
#define R_DT_Clog(p)	(lower_tail? R_D_LExp(p): R_D_log(p))/* log(1-p) in qF*/
#define R_DT_Log(p)	(lower_tail? (p) : R_Log1_Exp(p))
// ==   R_DT_log when we already "know" log_p == TRUE


#define R_Q_P01_check(p)			\
    if ((log_p	&& p > 0) ||			\
	(!log_p && (p < 0 || p > 1)) )		\
	ML_WARN_return_NAN

/* Do the boundaries exactly for q*() functions :
 * Often  _LEFT_ = ML_NEGINF , and very often _RIGHT_ = ML_POSINF;
 *
 * R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  :<==>
 *
 *     R_Q_P01_check(p);
 *     if (p == R_DT_0) return _LEFT_ ;
 *     if (p == R_DT_1) return _RIGHT_;
 *
 * the following implementation should be more efficient (less tests):
 */
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)		\
    if (log_p) {					\
	if(p > 0)					\
	    ML_WARN_return_NAN;				\
	if(p == 0) /* upper bound*/			\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
	if(p == ML_NEGINF)				\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
    }							\
    else { /* !log_p */					\
	if(p < 0 || p > 1)				\
	    ML_WARN_return_NAN;				\
	if(p == 0)					\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
	if(p == 1)					\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
    }

#define R_P_bounds_01(x, x_min, x_max)	\
    if(x <= x_min) return R_DT_0;		\
    if(x >= x_max) return R_DT_1
/* is typically not quite optimal for (-Inf,Inf) where
 * you'd rather have */
#define R_P_bounds_Inf_01(x)			\
    if(!R_FINITE(x)) {				\
	if (x > 0) return R_DT_1;		\
	/* x < 0 */return R_DT_0;		\
    }



/* additions for density functions (C.Loader) */
#define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
// version working with rf := sqrt(f) [avoiding overflow in computation of f in the caller]
#define R_D_rtxp(rf,x)    (give_log ? -log(rf)+(x) : exp(x)/(rf))

/* [neg]ative or [non int]eger : */
#define R_D_negInonint(x) (x < 0. || R_nonint(x))

// for discrete d<distr>(x, ...) :
#define R_D_nonint_check(x)				\
   if(R_nonint(x)) {					\
       MATHLIB_WARNING(_("non-integer x = %f"), x);	\
	return R_D__0;					\
   }


// @source_type: h
// @source_origin: refactored.h
// @provides: bessel_j_cycle_dependent, bessel_j_cycle_dependent_ex, bessel_j_cycle_free, bessel_j_cycle_free_ex, bessel_y_cycle_dependent, bessel_y_cycle_dependent_ex, bessel_y_cycle_free, bessel_y_cycle_free_ex, NMATH_REFACTORED_INTERNALS_H, stirlerr_cycle_dependent, stirlerr_cycle_free
// @all_depends_count: 0
// @load_order: 5

#ifndef NMATH_REFACTORED_INTERNALS_H
#define NMATH_REFACTORED_INTERNALS_H

/*
 * Internal declarations introduced by cycle-breaking refactors.
 * These are implementation details, not public API entry points.
 */

/* Stirlerr split internals */
attribute_hidden double stirlerr_cycle_free(double);
attribute_hidden double stirlerr_cycle_dependent(double);

/* Bessel split internals */
attribute_hidden double bessel_j_cycle_free(double, double);
attribute_hidden double bessel_j_cycle_free_ex(double, double, double *);
attribute_hidden double bessel_j_cycle_dependent(double, double);
attribute_hidden double bessel_j_cycle_dependent_ex(double, double, double *);

attribute_hidden double bessel_y_cycle_free(double, double);
attribute_hidden double bessel_y_cycle_free_ex(double, double, double *);
attribute_hidden double bessel_y_cycle_dependent(double, double);
attribute_hidden double bessel_y_cycle_dependent_ex(double, double, double *);

#endif /* NMATH_REFACTORED_INTERNALS_H */


// @source_type: h
// @source_origin: Rmath.h
// @includes: cmath, math.h, Boolean.h
// @provides: __STDC_WANT_IEC_60559_FUNCS_EXT__, bessel_i, bessel_i_ex, bessel_j, bessel_j_ex, bessel_k, bessel_k_ex, bessel_y, bessel_y_ex, beta, choose, cospi, dbeta, dbinom, dbinom_raw, dcauchy, dchisq, dexp, df, dgamma, dgeom, dhyper, digamma, dlnorm, dlogis, dnbeta, dnbinom, dnbinom_mu, dnchisq, dnf, dnorm, dnorm4, dnt, dpois, dpois_raw, dpsifn, dsignrank, dt, dtukey, dunif, dweibull, dwilcox, exp_rand , fmax2, fmin2, fprec, fround, fsign, ftrunc, gammafn, HAVE_EXPM1, HAVE_HYPOT, HAVE_LOG1P, HAVE_WORKING_LOG1P, imax2, imin2, lbeta, lchoose, lgamma1p, lgammafn, lgammafn_sign, log1mexp, log1p, log1pexp, log1pexp , log1pmx, log1pmx , logspace_add, logspace_sub, logspace_sum, M_1_PI, M_1_SQRT_2PI, M_2_PI, M_2_SQRTPI, M_2PI, M_E, M_LN_2PI, M_LN_SQRT_2PI, M_LN_SQRT_PI, M_LN_SQRT_PId2, M_LN10, M_LN2, M_LOG10_2, M_LOG10E, M_LOG2E, M_PI, M_PI_2, M_PI_4, M_SQRT_2dPI, M_SQRT_3, M_SQRT_32, M_SQRT_PI, M_SQRT1_2, M_SQRT2, norm_rand , pbeta, pbeta_raw, pbinom, pcauchy, pchisq, pentagamma, pexp, pf, pgamma, pgeom, phyper, plnorm, plogis, pnbeta, pnbinom, pnbinom_mu, pnchisq, pnf, pnorm, pnorm_both, pnorm5, pnt, pow1p, pow1p , ppois, psigamma, psignrank, pt, ptukey, punif, pweibull, pwilcox, qbeta, qbinom, qcauchy, qchisq, qchisq_appr, qexp, qf, qgamma, qgeom, qhyper, qlnorm, qlogis, qnbeta, qnbinom, qnbinom_mu, qnchisq, qnf, qnorm, qnorm5, qnt, qpois, qsignrank, qt, qtukey, qunif, qweibull, qwilcox, R_pow, R_pow_di, R_unif_index, R_VERSION_STRING, rbeta, rbinom, rcauchy, rchisq, rexp, rf, rgamma, rgeom, rhyper, rlnorm, Rlog1p, rlogis, RMATH_H, rmultinom, rnbeta, rnbinom, rnbinom_mu, rnchisq, rnf, rnorm, rnt, rpois, rsignrank, rt, Rtanpi , rtukey, runif, rweibull, rwilcox, sign, signrank_free , sinpi, tanpi, tetragamma, trigamma, unif_rand , wilcox_free 
// @all_depends_count: 0
// @load_order: 6

/* -*- C -*-
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2025  The R Core Team
 *  Copyright (C) 2004       The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *

 * Rmath.h  should contain ALL headers from R's C code in `src/nmath'
   -------  such that ``the Math library'' can be used by simply

   ``#include <Rmath.h> ''

   and nothing else.

   It is part of the API and supports 'standalone Rmath'.
   Some entries possibly are not yet documented in 'Writing R Extensions'.

*/
#ifndef RMATH_H
#define RMATH_H

/* needed for cospi etc */
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
# define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif

#if defined(__cplusplus) && !defined(DO_NOT_USE_CXX_HEADERS)
// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: # include <cmath>
#else
// openclport-disabled-include: # include <math.h>
#endif

#ifdef NO_C_HEADERS
# warning "use of NO_C_HEADERS is defunct and will be ignored"
#endif

/*  EXPAND CONFIGURE MACROS */
#define R_VERSION_STRING "4.6.0"

// Legacy defines -- C99 functions which R >= 3.5.0 requires
#ifndef HAVE_EXPM1
# define HAVE_EXPM1 1
#endif
#ifndef HAVE_HYPOT
# define HAVE_HYPOT 1
#endif
#ifndef HAVE_LOG1P
# define HAVE_LOG1P 1
#endif

#ifndef HAVE_WORKING_LOG1P
# define HAVE_WORKING_LOG1P 1
#endif

#if !defined(HAVE_WORKING_LOG1P)
/* remap to avoid problems with getting the right entry point */
extern "C" {
	double  Rlog1p(double);
}
#define log1p Rlog1p
#endif


/* ----- The following constants and entry points are part of the R API ---- */

/* 30 Decimal-place constants */
/* Computed with bc -l (scale=32; proper round) */

/* SVID & X/Open Constants */
/* Names from Solaris math.h */

#ifndef M_E
#define M_E		2.718281828459045235360287471353	/* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E		1.442695040888963407359924681002	/* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918917	/* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10		2.302585092994045684017991454684	/* ln(10) */
#endif

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.570796326794896619231321691640	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845820	/* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745	/* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490	/* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.128379167095512573896158903122	/* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.414213562373095048801688724210	/* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362105	/* 1/sqrt(2) */
#endif

/* R-Specific Constants */

#ifndef M_SQRT_3
#define M_SQRT_3	1.732050807568877293527446341506	/* sqrt(3) */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#endif


#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi))
								   == log(pi)/2 */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi))
								 == log(2*pi)/2 */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2))
								   == log(pi/2)/2 */
#endif


// openclport-disabled-include: # include <R_ext/Boolean.h>

#define bessel_i	Rf_bessel_i
#define bessel_j	Rf_bessel_j
#define bessel_k	Rf_bessel_k
#define bessel_y	Rf_bessel_y
#define bessel_i_ex	Rf_bessel_i_ex
#define bessel_j_ex	Rf_bessel_j_ex
#define bessel_k_ex	Rf_bessel_k_ex
#define bessel_y_ex	Rf_bessel_y_ex
#define beta		Rf_beta
#define choose		Rf_choose
#define dbeta		Rf_dbeta
#define dbinom		Rf_dbinom
#define dbinom_raw	Rf_dbinom_raw
#define dcauchy		Rf_dcauchy
#define dchisq		Rf_dchisq
#define dexp		Rf_dexp
#define df		Rf_df
#define dgamma		Rf_dgamma
#define dgeom		Rf_dgeom
#define dhyper		Rf_dhyper
#define digamma		Rf_digamma
#define dlnorm		Rf_dlnorm
#define dlogis		Rf_dlogis
#define dnbeta		Rf_dnbeta
#define dnbinom		Rf_dnbinom
#define dnbinom_mu	Rf_dnbinom_mu
#define dnchisq		Rf_dnchisq
#define dnf		Rf_dnf
#define dnorm4		Rf_dnorm4
#define dnt		Rf_dnt
#define dpois_raw	Rf_dpois_raw
#define dpois		Rf_dpois
#define dpsifn		Rf_dpsifn
#define dsignrank	Rf_dsignrank
#define dt		Rf_dt
#define dtukey		Rf_dtukey
#define dunif		Rf_dunif
#define dweibull	Rf_dweibull
#define dwilcox		Rf_dwilcox
#define fmax2		Rf_fmax2
#define fmin2		Rf_fmin2
#define fprec		Rf_fprec
#define fround		Rf_fround
#define ftrunc		Rf_ftrunc
#define fsign		Rf_fsign
#define gammafn		Rf_gammafn
#define imax2		Rf_imax2
#define imin2		Rf_imin2
#define lbeta		Rf_lbeta
#define lchoose		Rf_lchoose
#define lgammafn	Rf_lgammafn
#define lgammafn_sign	Rf_lgammafn_sign
#define lgamma1p	Rf_lgamma1p
#define pow1p		Rf_pow1p
#define log1mexp       	Rf_log1mexp
#define log1pexp       	Rf_log1pexp
#define log1pmx		Rf_log1pmx
#define logspace_add	Rf_logspace_add
#define logspace_sub	Rf_logspace_sub
#define logspace_sum	Rf_logspace_sum
#define pbeta		Rf_pbeta
#define pbeta_raw	Rf_pbeta_raw
#define pbinom		Rf_pbinom
#define pcauchy		Rf_pcauchy
#define pchisq		Rf_pchisq
#define pentagamma	Rf_pentagamma
#define pexp		Rf_pexp
#define pf		Rf_pf
#define pgamma		Rf_pgamma
#define pgeom		Rf_pgeom
#define phyper		Rf_phyper
#define plnorm		Rf_plnorm
#define plogis		Rf_plogis
#define pnbeta		Rf_pnbeta
#define pnbinom		Rf_pnbinom
#define pnbinom_mu     	Rf_pnbinom_mu
#define pnchisq		Rf_pnchisq
#define pnf		Rf_pnf
#define pnorm5		Rf_pnorm5
#define pnorm_both	Rf_pnorm_both
#define pnt		Rf_pnt
#define ppois		Rf_ppois
#define psignrank	Rf_psignrank
#define psigamma	Rf_psigamma
#define pt		Rf_pt
#define ptukey		Rf_ptukey
#define punif		Rf_punif
#define pweibull	Rf_pweibull
#define pwilcox		Rf_pwilcox
#define qbeta		Rf_qbeta
#define qbinom		Rf_qbinom
#define qcauchy		Rf_qcauchy
#define qchisq		Rf_qchisq
#define qchisq_appr	Rf_qchisq_appr
#define qexp		Rf_qexp
#define qf		Rf_qf
#define qgamma		Rf_qgamma
#define qgeom		Rf_qgeom
#define qhyper		Rf_qhyper
#define qlnorm		Rf_qlnorm
#define qlogis		Rf_qlogis
#define qnbeta		Rf_qnbeta
#define qnbinom		Rf_qnbinom
#define qnbinom_mu     	Rf_qnbinom_mu
#define qnchisq		Rf_qnchisq
#define qnf		Rf_qnf
#define qnorm5		Rf_qnorm5
#define qnt		Rf_qnt
#define qpois		Rf_qpois
#define qsignrank	Rf_qsignrank
#define qt		Rf_qt
#define qtukey		Rf_qtukey
#define qunif		Rf_qunif
#define qweibull	Rf_qweibull
#define qwilcox		Rf_qwilcox
#define rbeta		Rf_rbeta
#define rbinom		Rf_rbinom
#define rcauchy		Rf_rcauchy
#define rchisq		Rf_rchisq
#define rexp		Rf_rexp
#define rf		Rf_rf
#define rgamma		Rf_rgamma
#define rgeom		Rf_rgeom
#define rhyper		Rf_rhyper
#define rlnorm		Rf_rlnorm
#define rlogis		Rf_rlogis
#define rmultinom	Rf_rmultinom
#define rnbeta		Rf_rnbeta
#define rnbinom		Rf_rnbinom
#define rnbinom_mu     	Rf_rnbinom_mu
#define rnchisq		Rf_rnchisq
#define rnf		Rf_rnf
#define rnorm		Rf_rnorm
#define rnt		Rf_rnt
#define rpois		Rf_rpois
#define rsignrank	Rf_rsignrank
#define rt		Rf_rt
#define rtukey		Rf_rtukey
#define runif		Rf_runif
#define rweibull	Rf_rweibull
#define rwilcox		Rf_rwilcox
#define sign		Rf_sign
#define tetragamma	Rf_tetragamma
#define trigamma	Rf_trigamma

#define dnorm dnorm4
#define pnorm pnorm5
#define qnorm qnorm5

#ifdef  __cplusplus
extern "C" {
#endif
	/* R's versions with !R_FINITE checks */

	double R_pow(double x, double y);
	double R_pow_di(double, int);

	/* Random Number Generators */

	double	norm_rand(void); //not remapped
	double	unif_rand(void); //not remapped
	double  R_unif_index(double);
	double	exp_rand(void); //not remapped

	/* Normal Distribution */

	double	dnorm(double, double, double, int);
	double	pnorm(double, double, double, int, int);
	double	qnorm(double, double, double, int, int);
	double	rnorm(double, double);
	void	pnorm_both(double, double*, double*, int, int);/* both tails */

	/* Uniform Distribution */

	double	dunif(double, double, double, int);
	double	punif(double, double, double, int, int);
	double	qunif(double, double, double, int, int);
	double	runif(double, double);

	/* Gamma Distribution */

	double	dgamma(double, double, double, int);
	double	pgamma(double, double, double, int, int);
	double	qgamma(double, double, double, int, int);
	double	rgamma(double, double);

	double  log1pmx(double); /* Accurate log(1+x) - x, {care for small x} */
	double  log1pexp(double); // <-- ../nmath/plogis.c
	double  log1mexp(double);
	double  lgamma1p(double);/* accurate log(gamma(x+1)), small x (0 < x < 0.5) */

	double  pow1p(double, double); /* pow1p(x, y) := (1+x)^y  accurately also for |x| << 1 */

	/* Compute the log of a sum or difference from logs of terms, i.e.,
	 *
	 *     log (exp (logx) + exp (logy))
	 * or  log (exp (logx) - exp (logy))
	 *
	 * without causing overflows or throwing away too much accuracy:
	 */
	double  logspace_add(double logx, double logy);
	double  logspace_sub(double logx, double logy);
	double  logspace_sum(const double*, int);

	/* Beta Distribution */

	double	dbeta(double, double, double, int);
	double	pbeta(double, double, double, int, int);
	double	qbeta(double, double, double, int, int);
	double	rbeta(double, double);

	/* Lognormal Distribution */

	double	dlnorm(double, double, double, int);
	double	plnorm(double, double, double, int, int);
	double	qlnorm(double, double, double, int, int);
	double	rlnorm(double, double);

	/* Chi-squared Distribution */

	double	dchisq(double, double, int);
	double	pchisq(double, double, int, int);
	double	qchisq(double, double, int, int);
	double	rchisq(double);

	/* Non-central Chi-squared Distribution */

	double	dnchisq(double, double, double, int);
	double	pnchisq(double, double, double, int, int);
	double	qnchisq(double, double, double, int, int);
	double	rnchisq(double, double);

	/* F Distribution */

	double	df(double, double, double, int);
	double	pf(double, double, double, int, int);
	double	qf(double, double, double, int, int);
	double	rf(double, double);

	/* Student t Distribution */

	double	dt(double, double, int);
	double	pt(double, double, int, int);
	double	qt(double, double, int, int);
	double	rt(double);

	/* Binomial Distribution */

	double  dbinom_raw(double x, double n, double p, double q, int give_log);
	double	dbinom(double, double, double, int);
	double	pbinom(double, double, double, int, int);
	double	qbinom(double, double, double, int, int);
	double	rbinom(double, double);

	/* Multinomial Distribution */

	void	rmultinom(int, double*, int, int*);

	/* Cauchy Distribution */

	double	dcauchy(double, double, double, int);
	double	pcauchy(double, double, double, int, int);
	double	qcauchy(double, double, double, int, int);
	double	rcauchy(double, double);

	/* Exponential Distribution */

	double	dexp(double, double, int);
	double	pexp(double, double, int, int);
	double	qexp(double, double, int, int);
	double	rexp(double);

	/* Geometric Distribution */

	double	dgeom(double, double, int);
	double	pgeom(double, double, int, int);
	double	qgeom(double, double, int, int);
	double	rgeom(double);

	/* Hypergeometric Distribution */

	double	dhyper(double, double, double, double, int);
	double	phyper(double, double, double, double, int, int);
	double	qhyper(double, double, double, double, int, int);
	double	rhyper(double, double, double);

	/* Negative Binomial Distribution */

	double	dnbinom(double, double, double, int);
	double	pnbinom(double, double, double, int, int);
	double	qnbinom(double, double, double, int, int);
	double	rnbinom(double, double);

	double	dnbinom_mu(double, double, double, int);
	double	pnbinom_mu(double, double, double, int, int);
	double	qnbinom_mu(double, double, double, int, int);
	double	rnbinom_mu(double, double);

	/* Poisson Distribution */

	double	dpois_raw(double, double, int);
	double	dpois(double, double, int);
	double	ppois(double, double, int, int);
	double	qpois(double, double, int, int);
	double	rpois(double);

	/* Weibull Distribution */

	double	dweibull(double, double, double, int);
	double	pweibull(double, double, double, int, int);
	double	qweibull(double, double, double, int, int);
	double	rweibull(double, double);

	/* Logistic Distribution */

	double	dlogis(double, double, double, int);
	double	plogis(double, double, double, int, int);
	double	qlogis(double, double, double, int, int);
	double	rlogis(double, double);

	/* Non-central Beta Distribution */

	double	dnbeta(double, double, double, double, int);
	double	pnbeta(double, double, double, double, int, int);
	double	qnbeta(double, double, double, double, int, int);
#if defined(COMPILING_RCPP)
# if RCPP_VERSION <= Rcpp_Version(1,1,0)
	// Rcpp 1.1.0 relies on this declaration existing
	double rnbeta(double, double, double);
# endif
#elif defined(RcppCommon_h) || defined(Rcpp_Rmath_h) || defined(__cplusplus)
	// Rcpp revdeps may also need this declaration to exist for now
	double rnbeta(double, double, double);
# endif

	/* Non-central F Distribution */

	double  dnf(double, double, double, double, int);
	double	pnf(double, double, double, double, int, int);
	double	qnf(double, double, double, double, int, int);

	/* Non-central Student t Distribution */

	double	dnt(double, double, double, int);
	double	pnt(double, double, double, int, int);
	double	qnt(double, double, double, int, int);

	/* Studentized Range Distribution */

	double	ptukey(double, double, double, double, int, int);
	double	qtukey(double, double, double, double, int, int);

	/* Wilcoxon Rank Sum Distribution */

	double dwilcox(double, double, double, int);
	double pwilcox(double, double, double, int, int);
	double qwilcox(double, double, double, int, int);
	double rwilcox(double, double);
	void wilcox_free(void); // not remapped
	/* Wilcoxon Signed Rank Distribution */

	double dsignrank(double, double, int);
	double psignrank(double, double, int, int);
	double qsignrank(double, double, int, int);
	double rsignrank(double);
	void signrank_free(void); // not remapped

	/* Gamma and Related Functions */
	double	gammafn(double);
	double	lgammafn(double);
	double	lgammafn_sign(double, int*);
	void    dpsifn(double, int, int, int, double*, int*, int*);
	double	psigamma(double, double);
	double	digamma(double);
	double	trigamma(double);
	double	tetragamma(double);
	double	pentagamma(double);

	double	beta(double, double);
	double	lbeta(double, double);

	double	choose(double, double);
	double	lchoose(double, double);

	/* Bessel Functions */

	double	bessel_i(double, double, double);
	double	bessel_j(double, double);
	double	bessel_k(double, double, double);
	double	bessel_y(double, double);
	double	bessel_i_ex(double, double, double, double*);
	double	bessel_j_ex(double, double, double*);
	double	bessel_k_ex(double, double, double, double*);
	double	bessel_y_ex(double, double, double*);


	/* General Support Functions */

	int	imax2(int, int);
	int	imin2(int, int);
	double	fmax2(double, double);
	double	fmin2(double, double);
	double	sign(double);
	double	fprec(double, double);
	double	fround(double, double);
	double	fsign(double, double);
	double	ftrunc(double);

	/* More accurate cos(pi*x), sin(pi*x), tan(pi*x)

	   These declarations might clash with system headers if someone had
	   already included math.h with __STDC_WANT_IEC_60559_FUNCS_EXT__
	   defined (and we try, above).
	   We check for that via the value of __STDC_IEC_60559_FUNCS__
	*/
#if !defined(__OPENCL_VERSION__) && !defined(__OPENCL_C_VERSION__)
# if !(defined(__STDC_IEC_60559_FUNCS__) && __STDC_IEC_60559_FUNCS__ >= 201506L)
	double cospi(double);
	double sinpi(double);
	double tanpi(double);
# endif
	double Rtanpi(double); /* host-side own helper */
#endif


#if defined(R_NO_REMAP_RMATH)
#undef bessel_i
#undef bessel_j
#undef bessel_k
#undef bessel_y
#undef bessel_i_ex
#undef bessel_j_ex
#undef bessel_k_ex
#undef bessel_y_ex
#undef beta
#undef choose
#undef dbeta
#undef dbinom
#undef dbinom_raw
#undef dcauchy
#undef dchisq
#undef dexp
#undef df
#undef dgamma
#undef dgeom
#undef dhyper
#undef digamma
#undef dlnorm
#undef dlogis
#undef dnbeta
#undef dnbinom
#undef dnbinom_mu
#undef dnchisq
#undef dnf
#undef dnorm4
#undef dnt
#undef dpois_raw
#undef dpois
#undef dpsifn
#undef dsignrank
#undef dt
#undef dtukey
#undef dunif
#undef dweibull
#undef dwilcox
#undef fmax2
#undef fmin2
#undef fprec
#undef fround
#undef ftrunc
#undef fsign
#undef gammafn
#undef imax2
#undef imin2
#undef lbeta
#undef lchoose
#undef lgammafn
#undef lgammafn_sign
#undef lgamma1p
#undef pow1p
#undef log1mexp
#undef log1pexp
#undef log1pmx
#undef logspace_add
#undef logspace_sub
#undef logspace_sum
#undef pbeta
#undef pbeta_raw
#undef pbinom
#undef pcauchy
#undef pchisq
#undef pentagamma
#undef pexp
#undef pf
#undef pgamma
#undef pgeom
#undef phyper
#undef plnorm
#undef plogis
#undef pnbeta
#undef pnbinom
#undef pnbinom_mu
#undef pnchisq
#undef pnf
#undef pnorm5
#undef pnorm_both
#undef pnt
#undef ppois
#undef psignrank
#undef psigamma
#undef pt
#undef ptukey
#undef punif
#undef pweibull
#undef pwilcox
#undef qbeta
#undef qbinom
#undef qcauchy
#undef qchisq
#undef qchisq_appr
#undef qexp
#undef qf
#undef qgamma
#undef qgeom
#undef qhyper
#undef qlnorm
#undef qlogis
#undef qnbeta
#undef qnbinom
#undef qnbinom_mu
#undef qnchisq
#undef qnf
#undef qnorm5
#undef qnt
#undef qpois
#undef qsignrank
#undef qt
#undef qtukey
#undef qunif
#undef qweibull
#undef qwilcox
#undef rbeta
#undef rbinom
#undef rcauchy
#undef rchisq
#undef rexp
#undef rf
#undef rgamma
#undef rgeom
#undef rhyper
#undef rlnorm
#undef rlogis
#undef rmultinom
#undef rnbeta
#undef rnbinom
#undef rnbinom_mu
#undef rnchisq
#undef rnf
#undef rnorm
#undef rnt
#undef rpois
#undef rsignrank
#undef rt
#undef rtukey
#undef runif
#undef rweibull
#undef rwilcox
#undef sign
#undef tetragamma
#undef trigamma

#undef dnorm
#undef pnorm
#undef qnorm
#endif /* R_NO_REMAP_RMATH */

	/* ----------------- Private part of the header file ------------------- */


#ifdef  __cplusplus
}
#endif

#endif /* RMATH_H */


// @source_type: h
// @source_origin: nmath.h
// @includes: config.h, math.h, float.h, Rconfig.h, Rmath.h, RS.h, Print.h, Error.h, Arith.h, libintl.h
// @depends: Rmath
// @provides: _, __STDC_WANT_IEC_60559_FUNCS_EXT__, attribute_hidden, bd0, bratio, calloc, chebyshev_eval, chebyshev_init, ebd0, free, gammalims, LDOUBLE, lfastchoose, lgamma1p, lgammacor, lgammacor , log1pmx, MATHLIB_ERROR, MATHLIB_PRIVATE_H, MATHLIB_WARNING, MATHLIB_WARNING2, MATHLIB_WARNING3, MATHLIB_WARNING4, MATHLIB_WARNING5, MATHLIB_WARNING6, ME_DOMAIN, ME_NOCONV, ME_NONE, ME_PRECISION, ME_RANGE, ME_UNDERFLOW, ML_NAN, ML_NEGINF, ML_POSINF, ML_VALID, ML_WARN_return_NAN, ML_WARNING, pbeta_raw, pgamma_raw, pnbeta_raw, pnbeta2, pnchisq_raw, qchisq_appr, R_CheckUserInterrupt, R_forceint, R_nonint, Rf_d1mach, Rf_gamma_cody, Rf_i1mach, stirlerr, stirlerr  , WILCOX_MAX
// @all_depends_count: 1
// @all_depends: Rmath
// @load_order: 8

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2025  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Private header file for use during compilation of Mathlib */
#ifndef MATHLIB_PRIVATE_H
#define MATHLIB_PRIVATE_H

#ifdef HAVE_CONFIG_H
// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #  include <config.h>
#endif

/* Required by C99 but might be slow */
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#else
#  define LDOUBLE double
#endif

/* To ensure atanpi, cospi,  sinpi, tanpi are defined */
# ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
#  define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
# endif

// openclport-disabled-include: #include <math.h>
// openclport-disabled-include: #include <float.h> /* DBL_MIN etc */

// openclport-disabled-include: #include <Rconfig.h>
// openclport-disabled-include: #include <Rmath.h>

/* Used internally only */
double  Rf_d1mach(int);
double	Rf_gamma_cody(double);

// openclport-disabled-include: #include <R_ext/RS.h>

/* possibly needed for debugging */
// openclport-disabled-include: #include <R_ext/Print.h>

/* moved from dpq.h */
#if defined(__OPENCL_VERSION__) || defined(__OPENCL_C_VERSION__)
# define R_forceint(x)   rint(x)
#elif defined(HAVE_NEARBYINT)
# define R_forceint(x)   nearbyint(x)
#else
# define R_forceint(x)   round(x)
#endif
//R >= 3.1.0; previously: (fabs((x) - R_forceint(x)) > 1e-7)
//R >= 4.4.0; previously: (fabs((x) - R_forceint(x)) > 1e-7 * fmax2(1., fabs(x)))
# define R_nonint(x) 	  (fabs((x) - R_forceint(x)) > 1e-9 * fmax2(1., fabs(x)))
/*						       .... maybe change even to ~ 1e-11 or 12 */

// openclport-disabled-include: #include <R_ext/Error.h>
# define MATHLIB_ERROR(fmt,x)		error(fmt,x);
# define MATHLIB_WARNING(fmt,x)		warning(fmt,x)
# define MATHLIB_WARNING2(fmt,x,x2)	warning(fmt,x,x2)
# define MATHLIB_WARNING3(fmt,x,x2,x3)	warning(fmt,x,x2,x3)
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) warning(fmt,x,x2,x3,x4)
# define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) warning(fmt,x,x2,x3,x4,x5)
# define MATHLIB_WARNING6(fmt,x,x2,x3,x4,x5,x6) warning(fmt,x,x2,x3,x4,x5,x6)

// openclport-disabled-include: #include <R_ext/Arith.h>
#ifndef ML_POSINF
#define ML_POSINF	R_PosInf
#endif
#ifndef ML_NEGINF
#define ML_NEGINF	R_NegInf
#endif
#ifndef ML_NAN
#define ML_NAN		R_NaN
#endif


void R_CheckUserInterrupt(void);
/* Ei-ji Nakama reported that AIX 5.2 has calloc as a macro and objected
   to redefining it.  Tests added for 2.2.1 */
#ifdef calloc
# undef calloc
#endif
#define calloc R_chk_calloc
#ifdef free
# undef free
#endif
#define free R_chk_free

#ifdef ENABLE_NLS
// openclport-disabled-include: #include <libintl.h>
#define _(String) gettext (String)
#else
#define _(String) (String)
#endif

#define ML_VALID(x)	(!ISNAN(x))

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occurred (important for IEEE)*/


/* Device-only OpenCL ports: string literals live in __constant address space;
 * the upstream Mathlib ML_WARNING macro (char * + _("...")) is ill-formed here.
 * Toolchains sometimes omit __OPENCL_VERSION__, so branch on it is unreliable.
 * No-op warnings + ML_NAN-only return matches archived working packaged kernels.
 * Do not replicate this tweak in `<openclport>/nmath/` — only in this refactor tree. */
#define ML_WARNING(x, s) ((void)0)
#define ML_WARN_return_NAN { return ML_NAN; }

/* Wilcoxon Rank Sum Distribution */

#define WILCOX_MAX 50

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

/* Formerly private part of Mathlib.h */

/* always remap internal functions */
#define bd0       	Rf_bd0
#define ebd0       	Rf_ebd0
#define chebyshev_eval	Rf_chebyshev_eval
#define chebyshev_init	Rf_chebyshev_init
#define gammalims	Rf_gammalims
#define lfastchoose	Rf_lfastchoose
#define lgammacor	Rf_lgammacor
#define stirlerr       	Rf_stirlerr
#define pnchisq_raw   	Rf_pnchisq_raw
#define pgamma_raw   	Rf_pgamma_raw
#define pnbeta_raw   	Rf_pnbeta_raw
#define pnbeta2       	Rf_pnbeta2
#define bratio       	Rf_bratio

	/* Chebyshev Series */

attribute_hidden int chebyshev_init(double*, int, double);
attribute_hidden double chebyshev_eval(double, const double *, const int);

	/* Gamma and Related Functions */

attribute_hidden void gammalims(double*, double*);
attribute_hidden double lgammacor(double); /* log(gamma) correction */
attribute_hidden double stirlerr(double);  /* Stirling expansion "error" */

attribute_hidden double lfastchoose(double, double);

attribute_hidden double bd0(double, double);
attribute_hidden void ebd0(double, double, double*, double*);
attribute_hidden double log1pmx(double);
attribute_hidden double lgamma1p(double);

attribute_hidden double pnchisq_raw(double, double, double, double, double,
				     int, Rboolean, Rboolean);
attribute_hidden double pgamma_raw(double, double, int, int);
attribute_hidden double pbeta_raw(double, double, double, int, int);
attribute_hidden double qchisq_appr(double, double, double, int, int, double tol);
attribute_hidden LDOUBLE pnbeta_raw(double, double, double, double, double);
attribute_hidden double pnbeta2(double, double, double, double, double, int, int);

int	Rf_i1mach(int);

/* From toms708.c */
attribute_hidden void bratio(double a, double b, double x, double y,
	    		     double *w, double *w1, int *ierr, int log_p);


#endif /* MATHLIB_PRIVATE_H */


// @source_type: c
// @source_origin: stirlerr_cycle_free.c
// @includes: nmath.h
// @depends: nmath
// @provides: stirlerr_cycle_free
// @all_depends_count: 2
// @all_depends: Rmath, nmath
// @load_order: 22
// @local_macros: S0, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16

// openclport: macro hygiene pre-clean for concatenated translation units.
#ifdef S0
# undef S0
#endif
#ifdef S1
# undef S1
#endif
#ifdef S2
# undef S2
#endif
#ifdef S3
# undef S3
#endif
#ifdef S4
# undef S4
#endif
#ifdef S5
# undef S5
#endif
#ifdef S6
# undef S6
#endif
#ifdef S7
# undef S7
#endif
#ifdef S8
# undef S8
#endif
#ifdef S9
# undef S9
#endif
#ifdef S10
# undef S10
#endif
#ifdef S11
# undef S11
#endif
#ifdef S12
# undef S12
#endif
#ifdef S13
# undef S13
#endif
#ifdef S14
# undef S14
#endif
#ifdef S15
# undef S15
#endif
#ifdef S16
# undef S16
#endif

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */
#define S5 0.0019175269175269175269175262
#define S6 0.0064102564102564102564102561
#define S7 0.029550653594771241830065352
#define S8 0.17964437236883057316493850
#define S9 1.3924322169059011164274315
#define S10 13.402864044168391994478957
#define S11 156.84828462600201730636509
#define S12 2193.1033333333333333333333
#define S13 36108.771253724989357173269
#define S14 691472.26885131306710839498
#define S15 15238221.539407416192283370
#define S16 382900751.39141414141414141

static const double sferr_halves[31] = {
    0.0,
    0.1534264097200273452913848, 0.0810614667953272582196702,
    0.0548141210519176538961390, 0.0413406959554092940938221,
    0.03316287351993628748511048, 0.02767792568499833914878929,
    0.02374616365629749597132920, 0.02079067210376509311152277,
    0.01848845053267318523077934, 0.01664469118982119216319487,
    0.01513497322191737887351255, 0.01387612882307074799874573,
    0.01281046524292022692424986, 0.01189670994589177009505572,
    0.01110455975820691732662991, 0.010411265261972096497478567,
    0.009799416126158803298389475, 0.009255462182712732917728637,
    0.008768700134139385462952823, 0.008330563433362871256469318,
    0.007934114564314020547248100, 0.007573675487951840794972024,
    0.007244554301320383179543912, 0.006942840107209529865664152,
    0.006665247032707682442354394, 0.006408994188004207068439631,
    0.006171712263039457647532867, 0.005951370112758847735624416,
    0.005746216513010115682023589, 0.005554733551962801371038690
};

attribute_hidden double stirlerr_cycle_free(double n)
{
    double nn;

    if (n <= 23.5) {
        nn = n + n;
        if (n <= 15. && (nn == (int)nn)) return sferr_halves[(int)nn];
        if (n <= 5.25) {
            if (n >= 1.) {
                double l_n = log(n);
                return lgamma(n) + n * (1 - l_n) + ldexp(l_n - M_LN_2PI, -1);
            }
            /* For n < 1, wrapper routes to cycle_dependent branch. */
            return lgamma(n + 1.) - (n + 0.5) * log(n) + n - M_LN_SQRT_2PI;
        }
        nn = n * n;
        if (n > 12.8) return (S0-(S1-(S2-(S3-(S4-(S5 -S6/nn)/nn)/nn)/nn)/nn)/nn)/n;
        if (n > 12.3) return (S0-(S1-(S2-(S3-(S4-(S5-(S6 -S7/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
        if (n > 8.9) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7 -S8/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
        if (n > 7.3) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-S10/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
        if (n > 6.6) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-S12/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
        if (n > 6.1) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-S14/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
        return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-S16/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
    }

    nn = n * n;
    if (n > 15.7e6) return S0 / n;
    if (n > 6180) return (S0 - S1 / nn) / n;
    if (n > 205) return (S0 - (S1 - S2 / nn) / nn) / n;
    if (n > 86) return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
    if (n > 27) return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
    return (S0 - (S1 - (S2 - (S3 - (S4 - S5 / nn) / nn) / nn) / nn) / nn) / n;
}

// openclport: macro hygiene post-clean for concatenated translation units.
#undef S0
#undef S1
#undef S2
#undef S3
#undef S4
#undef S5
#undef S6
#undef S7
#undef S8
#undef S9
#undef S10
#undef S11
#undef S12
#undef S13
#undef S14
#undef S15
#undef S16


// @source_type: c
// @source_origin: chebyshev.c
// @includes: nmath.h
// @depends: nmath
// @provides: chebyshev_eval, chebyshev_init
// @all_depends_count: 2
// @all_depends: Rmath, nmath
// @load_order: 24

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    int chebyshev_init(double *dos, int nos, double eta)
 *    double chebyshev_eval(double x, double *a, int n)
 *
 *  DESCRIPTION
 *
 *    "chebyshev_init" determines the number of terms for the
 *    double precision orthogonal series "dos" needed to insure
 *    the error is no larger than "eta".  Ordinarily eta will be
 *    chosen to be one-tenth machine precision.
 *
 *    "chebyshev_eval" evaluates the n-term Chebyshev series
 *    "a" at "x".
 *
 *  NOTES
 *
 *    These routines are translations into C of Fortran routines
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    Based on the Fortran routine dcsevl by W. Fullerton.
 *    Adapted from R. Broucke, Algorithm 446, CACM., 16, 254 (1973).
 */

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"

/* NaNs propagated correctly */


attribute_hidden int chebyshev_init(double *dos, int nos, double eta)
{
    int i, ii;
    double err;

    if (nos < 1)
	return 0;

    err = 0.0;
    i = 0;			/* just to avoid compiler warnings */
    for (ii=1; ii<=nos; ii++) {
	i = nos - ii;
	err += fabs(dos[i]);
	if (err > eta) {
	    return i;
	}
    }
    return i;
}


attribute_hidden double chebyshev_eval(double x, const double *a, const int n)
{
    double b0, b1, b2, twox;
    int i;

    if (n < 1 || n > 1000) ML_WARN_return_NAN;

    if (x < -1.1 || x > 1.1) ML_WARN_return_NAN;

    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) {
	b2 = b1;
	b1 = b0;
	b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}


// @source_type: c
// @source_origin: cospi.c
// @includes: nmath.h
// @depends: nmath
// @provides: cospi, Rtanpi, sinpi, tanpi
// @all_depends_count: 2
// @all_depends: Rmath, nmath
// @load_order: 25

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2013-2022 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 */

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"

/* OpenCL C provides sinpi/cospi/tanpi builtins; avoid redeclaring/defining
   host-style variants here when this file is ported to .cl. */
#ifndef __OPENCL_VERSION__

/* HAVE_COSPI etc will not be defined in standalone-use: the
   intention is to make the versions here available in that case.

   The __cospi etc variants are from macOS (and perhaps other BSD-based systems).
*/

#ifdef HAVE_COSPI
#elif defined HAVE___COSPI
double cospi(double x) {
    return __cospi(x);
}
#else
// cos(pi * x)  -- exact when x = k/2  for all integer k
double cospi(double x) {
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x)) return x;
#endif
    if(!R_FINITE(x)) ML_WARN_return_NAN;

    x = fmod(fabs(x), 2.);// cos() symmetric; cos(pi(x + 2k)) == cos(pi x) for all integer k
    if(fmod(x, 1.) == 0.5) return 0.;
    if( x == 1.)	return -1.;
    if( x == 0.)	return  1.;
    // otherwise
    return cos(M_PI * x);
}
#endif

#ifdef HAVE_SINPI
#elif defined HAVE___SINPI
double sinpi(double x) {
    return __sinpi(x);
}
#else
// sin(pi * x)  -- exact when x = k/2  for all integer k
double sinpi(double x) {
#ifdef IEEE_754
    if (ISNAN(x)) return x;
#endif
    if(!R_FINITE(x)) ML_WARN_return_NAN;

    x = fmod(x, 2.); // sin(pi(x + 2k)) == sin(pi x)  for all integer k
    // map (-2,2) --> (-1,1] :
    if(x <= -1) x += 2.; else if (x > 1.) x -= 2.;
    if(x == 0. || x == 1.) return 0.;
    if(x ==  0.5)	return  1.;
    if(x == -0.5)	return -1.;
    // otherwise
    return sin(M_PI * x);
}
#endif

// tan(pi * x)  -- exact when x = k/4  for all integer k and half-values give NaN
// ----------- e.g. used in ../main/arithmetic.c : 
double Rtanpi(double x)
{
#ifdef IEEE_754
    if (ISNAN(x)) return x;
#endif
    if(!R_FINITE(x)) ML_WARN_return_NAN;

    x = fmod(x, 1.); // tan(pi(x + k)) == tan(pi x)  for all integer k
    // map (-1,1] --> (-1/2, 1/2] :
    if(x <= -0.5) x++; else if(x > 0.5) x--;
    return (x == 0.) ? 0. :
	((x ==  0.5 ) ? ML_NAN :
	((x ==  0.25) ?  1. :
	((x == -0.25) ? -1. :
			tan(M_PI * x)
	    )));
}

#if defined(HAVE_TANPI) || defined(HAVE___TANPI)
#else
double tanpi(double x) {
    return Rtanpi(x);
}
#endif

#if !defined(HAVE_TANPI) && defined(HAVE___TANPI)
double tanpi(double x) {
    return __tanpi(x);
}
/* #else tanpi() defined from C standard math lib */
#endif

#endif /* !__OPENCL_VERSION__ */


// @source_type: c
// @source_origin: fmax2.c
// @includes: nmath.h
// @depends: nmath
// @provides: fmax2
// @all_depends_count: 2
// @all_depends: Rmath, nmath
// @load_order: 34

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"

double fmax2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? y : x;
}


// @source_type: c
// @source_origin: gammalims.c
// @includes: nmath.h
// @depends: fmax2, nmath
// @provides: gammalims
// @all_depends_count: 3
// @all_depends: Rmath, nmath, fmax2
// @load_order: 41

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 1999-2025  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    void gammalims(double *xmin, double *xmax);
 *
 *  DESCRIPTION
 *
 *    This function calculates the minimum and maximum legal bounds
 *    for x in gammafn(x).  These are not the only bounds, but they
 *    are the only non-trivial ones to calculate.
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 */

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"

attribute_hidden void gammalims(double *xmin, double *xmax)
{
/* FIXME: Even better: If IEEE, #define these in nmath.h
	  and don't call gammalims() at all
*/
#ifdef IEEE_754
    *xmin = -170.5674972726612;
    *xmax =  171.61447887182298;/*(3 Intel/Sparc architectures)*/
#else
    double alnbig, alnsml, xln, xold;
    int i;

    alnsml = log(d1mach(1));
    *xmin = -alnsml;
    for (i=1; i<=10; ++i) {
	xold = *xmin;
	xln = log(*xmin);
	*xmin -= *xmin * ((*xmin + .5) * xln - *xmin - .2258 + alnsml) /
		(*xmin * xln + .5);
	if (fabs(*xmin - xold) < .005) {
	    *xmin = -(*xmin) + .01;
	    goto find_xmax;
	}
    }

    /* unable to find xmin */

    ML_WARNING(ME_NOCONV, "gammalims");
    *xmin = *xmax = ML_NAN;

find_xmax:

    alnbig = log(d1mach(2));
    *xmax = alnbig;
    for (i=1; i<=10; ++i) {
	xold = *xmax;
	xln = log(*xmax);
	*xmax -= *xmax * ((*xmax - .5) * xln - *xmax + .9189 - alnbig) /
		(*xmax * xln - .5);
	if (fabs(*xmax - xold) < .005) {
	    *xmax += -.01;
	    goto done;
	}
    }

    /* unable to find xmax */

    ML_WARNING(ME_NOCONV, "gammalims");
    *xmin = *xmax = ML_NAN;

done:
    *xmin = fmax2(*xmin, -(*xmax) + 1);
#endif
}



// @source_type: c
// @source_origin: lgammacor.c
// @includes: nmath.h
// @depends: chebyshev, nmath
// @provides: lgammacor
// @all_depends_count: 3
// @all_depends: Rmath, nmath, chebyshev
// @load_order: 45
// @local_macros: nalgm, xbig, xmax

// openclport: macro hygiene pre-clean for concatenated translation units.
#ifdef nalgm
# undef nalgm
#endif
#ifdef xbig
# undef xbig
#endif
#ifdef xmax
# undef xmax
#endif

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000-2021 The R Core Team
 *  Copyright (C) 1998 Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double lgammacor(double x);
 *
 *  DESCRIPTION
 *
 *    Compute the log gamma correction factor for x >= 10 so that
 *                                               ---------
 *
 *    log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x)
 *
 *    [ lgammacor(x) is called	Del(x)	in other contexts (e.g. dcdflib)], or  stirlerr(x)
 *				~~~~~~					       ~~~~~~~~~~~
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    written by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *  SEE ALSO
 *
 *    Loader(1999)'s stirlerr() {in ./stirlerr.c} is *very* similar in spirit,
 *    is faster and cleaner, but is only defined "fast" for half integers.
 */

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"

attribute_hidden double lgammacor(double x)
{
    const static double algmcs[15] = {  // below, nalgm = 5 ==> only the first 5 are used!
	+.1666389480451863247205729650822e+0,
	-.1384948176067563840732986059135e-4,
	+.9810825646924729426157171547487e-8,
	-.1809129475572494194263306266719e-10,
	+.6221098041892605227126015543416e-13,
	-.3399615005417721944303330599666e-15,
	+.2683181998482698748957538846666e-17,
	-.2868042435334643284144622399999e-19,
	+.3962837061046434803679306666666e-21,
	-.6831888753985766870111999999999e-23,
	+.1429227355942498147573333333333e-24,
	-.3547598158101070547199999999999e-26,
	+.1025680058010470912000000000000e-27,
	-.3401102254316748799999999999999e-29,
	+.1276642195630062933333333333333e-30
    };

/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 *   xbig = 2 ^ 26.5
 *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
#define nalgm 5
#define xbig  94906265.62425156

    if (x < 10) // possibly consider stirlerr()
	ML_WARN_return_NAN
#ifndef IEEE_754
#   define xmax  3.745194030963158e306
    else if (x >= xmax) {
	ML_WARNING(ME_UNDERFLOW, "lgammacor");
	/* allow to underflow below */
    }
#endif
    else if (x < xbig) {
	double tmp = 10 / x;
	return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
    }
    // x >= xbig
    return 1 / (x * 12);
}

// openclport: macro hygiene post-clean for concatenated translation units.
#undef nalgm
#undef xbig
#undef xmax


// @source_type: c
// @source_origin: log1p.c
// @includes: config.h, nmath.h
// @depends: chebyshev, nmath
// @provides: log1p, Rlog1p
// @all_depends_count: 3
// @all_depends: Rmath, nmath, chebyshev
// @load_order: 46
// @local_macros: nlnrel

// openclport: macro hygiene pre-clean for concatenated translation units.
#ifdef nlnrel
# undef nlnrel
#endif

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2018 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *	#include <Rmath.h>
 *	double log1p(double x);
 *
 *  DESCRIPTION
 *
 *	Compute the relative error logarithm.
 *
 *			log(1 + x)
 *
 *  NOTES
 *
 *	This code is a translation of the Fortran subroutine `dlnrel'
 *	written by W. Fullerton of Los Alamos Scientific Laboratory.
 */

/* Every currently known platform has log1p (which is C99), 
   but NetBSD/OpenBSD were at least at one time inaccurate */
#ifdef HAVE_CONFIG_H
// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: # include <config.h>
#endif
// openclport-disabled-include: #include "nmath.h"

#ifndef HAVE_WORKING_LOG1P
double Rlog1p(double x)
{
    /* series for log1p on the interval -.375 to .375
     *				     with weighted error   6.35e-32
     *				      log weighted error  31.20
     *			    significant figures required  30.93
     *				 decimal places required  32.01
     */
    const static double alnrcs[43] = {
	+.10378693562743769800686267719098e+1,
	-.13364301504908918098766041553133e+0,
	+.19408249135520563357926199374750e-1,
	-.30107551127535777690376537776592e-2,
	+.48694614797154850090456366509137e-3,
	-.81054881893175356066809943008622e-4,
	+.13778847799559524782938251496059e-4,
	-.23802210894358970251369992914935e-5,
	+.41640416213865183476391859901989e-6,
	-.73595828378075994984266837031998e-7,
	+.13117611876241674949152294345011e-7,
	-.23546709317742425136696092330175e-8,
	+.42522773276034997775638052962567e-9,
	-.77190894134840796826108107493300e-10,
	+.14075746481359069909215356472191e-10,
	-.25769072058024680627537078627584e-11,
	+.47342406666294421849154395005938e-12,
	-.87249012674742641745301263292675e-13,
	+.16124614902740551465739833119115e-13,
	-.29875652015665773006710792416815e-14,
	+.55480701209082887983041321697279e-15,
	-.10324619158271569595141333961932e-15,
	+.19250239203049851177878503244868e-16,
	-.35955073465265150011189707844266e-17,
	+.67264542537876857892194574226773e-18,
	-.12602624168735219252082425637546e-18,
	+.23644884408606210044916158955519e-19,
	-.44419377050807936898878389179733e-20,
	+.83546594464034259016241293994666e-21,
	-.15731559416479562574899253521066e-21,
	+.29653128740247422686154369706666e-22,
	-.55949583481815947292156013226666e-23,
	+.10566354268835681048187284138666e-23,
	-.19972483680670204548314999466666e-24,
	+.37782977818839361421049855999999e-25,
	-.71531586889081740345038165333333e-26,
	+.13552488463674213646502024533333e-26,
	-.25694673048487567430079829333333e-27,
	+.48747756066216949076459519999999e-28,
	-.92542112530849715321132373333333e-29,
	+.17578597841760239233269760000000e-29,
	-.33410026677731010351377066666666e-30,
	+.63533936180236187354180266666666e-31,
    };

#ifdef NOMORE_FOR_THREADS
    static int nlnrel = 0;
    static double xmin = 0.0;

    if (xmin == 0.0) xmin = -1 + sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)); */
    if (nlnrel == 0) /* initialize chebyshev coefficients */
	nlnrel = chebyshev_init(alnrcs, 43, DBL_EPSILON/20);/*was .1*d1mach(3)*/
#else
# define nlnrel 22
    const static double xmin = -0.999999985;
/* 22: for IEEE double precision where DBL_EPSILON =  2.22044604925031e-16 */
#endif

    if (x == 0.) return 0.;/* speed */
    if (x == -1) return(ML_NEGINF);
    if (x  < -1) ML_WARN_return_NAN;

    if (fabs(x) <= .375) {
        /* Improve on speed (only);
	   again give result accurate to IEEE double precision: */
	if(fabs(x) < .5 * DBL_EPSILON)
	    return x;

	if( (0 < x && x < 1e-8) || (-1e-9 < x && x < 0))
	    return x * (1 - .5 * x);
	/* else */
	return x * (1 - x * chebyshev_eval(x / .375, alnrcs, nlnrel));
    }
    /* else */
    if (x < xmin) {
	/* answer less than half precision because x too near -1 */
	ML_WARNING(ME_PRECISION, "log1p");
    }
    return log(1 + x);
}
#endif

// openclport: macro hygiene post-clean for concatenated translation units.
#undef nlnrel


// @source_type: c
// @source_origin: gamma.c
// @includes: nmath.h, refactored.h
// @depends: chebyshev, cospi, fmax2, gammalims, lgammacor, stirlerr_cycle_free, nmath, refactored
// @provides: gammafn
// @all_depends_count: 9
// @all_depends: refactored, Rmath, nmath, stirlerr_cycle_free, chebyshev, cospi, fmax2, gammalims, lgammacor
// @load_order: 61
// @local_macros: ngam, xmin, xmax, xsml, dxrel

// openclport: macro hygiene pre-clean for concatenated translation units.
#ifdef ngam
# undef ngam
#endif
#ifdef xmin
# undef xmin
#endif
#ifdef xmax
# undef xmax
#endif
#ifdef xsml
# undef xsml
#endif
#ifdef dxrel
# undef dxrel
#endif

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000-2024 The R Core Team
 *  Copyright (C) 2002-2024 The R Foundation
 *  Copyright (C) 1998 Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double gammafn(double x);
 *
 *  DESCRIPTION
 *
 *    This function computes the value of the gamma function.
 *
 *  NOTES
 *
 *    This function is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *    (e.g. http://www.netlib.org/slatec/fnlib/gamma.f)
 *
 *    The accuracy of this routine compares (very) favourably
 *    with those of the Sun Microsystems portable mathematical
 *    library.
 *
 *    MM specialized the case of  n!  for n < 50 - for even better precision
 */

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"
// openclport-disabled-include: #include "refactored.h"

double gammafn(double x)
{
    const static double gamcs[42] = {
	+.8571195590989331421920062399942e-2,
	+.4415381324841006757191315771652e-2,
	+.5685043681599363378632664588789e-1,
	-.4219835396418560501012500186624e-2,
	+.1326808181212460220584006796352e-2,
	-.1893024529798880432523947023886e-3,
	+.3606925327441245256578082217225e-4,
	-.6056761904460864218485548290365e-5,
	+.1055829546302283344731823509093e-5,
	-.1811967365542384048291855891166e-6,
	+.3117724964715322277790254593169e-7,
	-.5354219639019687140874081024347e-8,
	+.9193275519859588946887786825940e-9,
	-.1577941280288339761767423273953e-9,
	+.2707980622934954543266540433089e-10,
	-.4646818653825730144081661058933e-11,
	+.7973350192007419656460767175359e-12,
	-.1368078209830916025799499172309e-12,
	+.2347319486563800657233471771688e-13,
	-.4027432614949066932766570534699e-14,
	+.6910051747372100912138336975257e-15,
	-.1185584500221992907052387126192e-15,
	+.2034148542496373955201026051932e-16,
	-.3490054341717405849274012949108e-17,
	+.5987993856485305567135051066026e-18,
	-.1027378057872228074490069778431e-18,
	+.1762702816060529824942759660748e-19,
	-.3024320653735306260958772112042e-20,
	+.5188914660218397839717833550506e-21,
	-.8902770842456576692449251601066e-22,
	+.1527474068493342602274596891306e-22,
	-.2620731256187362900257328332799e-23,
	+.4496464047830538670331046570666e-24,
	-.7714712731336877911703901525333e-25,
	+.1323635453126044036486572714666e-25,
	-.2270999412942928816702313813333e-26,
	+.3896418998003991449320816639999e-27,
	-.6685198115125953327792127999999e-28,
	+.1146998663140024384347613866666e-28,
	-.1967938586345134677295103999999e-29,
	+.3376448816585338090334890666666e-30,
	-.5793070335782135784625493333333e-31
    };

#ifdef NOMORE_FOR_THREADS
    static int ngam = 0;
    static double xmin = 0, xmax = 0., xsml = 0., dxrel = 0.;

    /* Initialize machine dependent constants, the first time gamma() is called.
	FIXME for threads ! */
    if (ngam == 0) {
	ngam = chebyshev_init(gamcs, 42, DBL_EPSILON/20);/*was .1*d1mach(3)*/
	gammalims(&xmin, &xmax);/*-> ./gammalims.c */
	xsml = exp(fmax2(log(DBL_MIN), -log(DBL_MAX)) + 0.01);
	/*   = exp(.01)*DBL_MIN = 2.247e-308 for IEEE */
	dxrel = sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)) */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 * (xmin, xmax) are non-trivial, see ./gammalims.c
 * xsml = exp(.01)*DBL_MIN
 * dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
*/
# define ngam 22
# define xmin -170.5674972726612
# define xmax  171.61447887182298
# define xsml 2.2474362225598545e-308
# define dxrel 1.490116119384765696e-8
#endif

    if(ISNAN(x)) return x;

    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
    if (x == 0 || (x < 0 && x == round(x))) {
	ML_WARNING(ME_DOMAIN, "gammafn");
	return ML_NAN;
    }

    int i;
    double y = fabs(x), value;

    if (y <= 10) {

	/* Compute gamma(x) for -10 <= x <= 10
	 * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
	 * first of all. */

	int n = (int) x;
	if(x < 0) --n;
	y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
	--n;
	value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
	if (n == 0)
	    return value;/* x = 1.dddd = 1+y */

	if (n < 0) {
	    /* compute gamma(x) for -10 <= x < 1 */

	    /* exact 0 or "-n" checked already above */

	    /* The answer is less than half precision */
	    /* because x too near a negative integer. */
	    if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
		ML_WARNING(ME_PRECISION, "gammafn");
	    }

	    /* The argument is so close to 0 that the result would overflow. */
	    if (y < xsml) {
		ML_WARNING(ME_RANGE, "gammafn");
		if(x > 0) return ML_POSINF;
		else return ML_NEGINF;
	    }

	    n = -n;

	    for (i = 0; i < n; i++) {
		value /= (x + i);
	    }
	    return value;
	}
	else {
	    /* gamma(x) for 2 <= x <= 10 */

	    for (i = 1; i <= n; i++) {
		value *= (y + i);
	    }
	    return value;
	}
    }
    else {
	/* gamma(x) for	 y = |x| > 10. */

	if (x > xmax) {			/* Overflow */
	    // No warning: +Inf is the best answer
	    return ML_POSINF;
	}

	if (x < xmin) {			/* Underflow */
	    // No warning: 0 is the best answer
	    return 0.;
	}

	if(y <= 50 && y == (int)y) { /* compute (n - 1)! */
	    value = 1.;
	    for (i = 2; i < y; i++) value *= i;
	}
	else { /* normal case */
	    value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI +
			((2*y == (int)2*y) ? stirlerr_cycle_free(y) : lgammacor(y)));
	}

	if (x > 0)
	    return value;
	// else:  x < 0, not an integer :

	if (fabs((x - (int)(x - 0.5))/x) < dxrel) {
	    /* The answer is less than half precision because */
	    /* the argument is too near a negative integer. */

	    ML_WARNING(ME_PRECISION, "gammafn");
	}

	double sinpiy = sinpi(y);
	if (sinpiy == 0) {		/* Negative integer arg - overflow */
	    ML_WARNING(ME_RANGE, "gammafn");
	    return ML_POSINF;
	}

	return -M_PI / (y * sinpiy * value);
    }
}

// openclport: macro hygiene post-clean for concatenated translation units.
#undef ngam
#undef xmin
#undef xmax
#undef xsml
#undef dxrel


// @source_type: c
// @source_origin: lgamma.c
// @includes: nmath.h
// @depends: cospi, gamma, lgammacor, nmath
// @provides: lgammafn, lgammafn_sign
// @all_depends_count: 10
// @all_depends: refactored, Rmath, nmath, stirlerr_cycle_free, chebyshev, cospi, fmax2, gammalims, lgammacor, gamma
// @load_order: 62
// @local_macros: xmax, dxrel

// openclport: macro hygiene pre-clean for concatenated translation units.
#ifdef xmax
# undef xmax
#endif
#ifdef dxrel
# undef dxrel
#endif

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000-2020 The R Core Team
 *  Copyright (C) 1998 Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double lgammafn_sign(double x, int *sgn);
 *    double lgammafn(double x);
 *
 *  DESCRIPTION
 *
 *    The function lgammafn computes log|gamma(x)|.  The function
 *    lgammafn_sign in addition assigns the sign of the gamma function
 *    to the address in the second argument if this is not NULL.
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    The accuracy of this routine compares (very) favourably
 *    with those of the Sun Microsystems portable mathematical
 *    library.
 *
 *  ./toms708.c  has  gamln()
 */

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"

double lgammafn_sign(double x, int *sgn)
{
    double ans, y, sinpiy;

#ifdef NOMORE_FOR_THREADS
    static double xmax = 0.;
    static double dxrel = 0.;

    if (xmax == 0) {/* initialize machine dependent constants _ONCE_ */
	xmax = d1mach(2)/log(d1mach(2));/* = 2.533 e305	 for IEEE double */
	dxrel = sqrt (d1mach(4));/* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
   xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
   dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
 */
#define xmax  2.5327372760800758e+305
#define dxrel 1.490116119384765625e-8
#endif

    if (sgn != NULL) *sgn = 1;

#ifdef IEEE_754
    if(ISNAN(x)) return x;
#endif

    if (sgn != NULL && x < 0 && fmod(floor(-x), 2.) == 0)
	*sgn = -1;

    if (x <= 0 && x == trunc(x)) { /* Negative integer argument */
	// No warning: this is the best answer; was  ML_WARNING(ME_RANGE, "lgamma");
	return ML_POSINF;/* +Inf, since lgamma(x) = log|gamma(x)| */
    }

    y = fabs(x);

    if (y < 1e-306) return -log(y); // denormalized range, R change
    if (y <= 10) return log(fabs(gammafn(x)));
    /*
      ELSE  y = |x| > 10 ---------------------- */

    if (y > xmax) {
	// No warning: +Inf is the best answer
	return ML_POSINF;
    }

    if (x > 0) { /* i.e. y = x > 10 */
#ifdef IEEE_754
	if(x > 1e17)
	    return(x*(log(x) - 1.));
	else if(x > 4934720.)
	    return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
	else
#endif
	    return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
    }
    /* else: x < -10; y = -x */
    sinpiy = fabs(sinpi(y));

    if (sinpiy == 0) { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
	MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
	ML_WARN_return_NAN;
    }

    ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);

    if(fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) {

	/* The answer is less than half precision because
	 * the argument is too near a negative integer; e.g. for  lgamma(1e-7 - 11) */

	ML_WARNING(ME_PRECISION, "lgamma");
    }

    return ans;
}

double lgammafn(double x)
{
    return lgammafn_sign(x, NULL);
}

// openclport: macro hygiene post-clean for concatenated translation units.
#undef xmax
#undef dxrel


// @source_type: c
// @source_origin: pgamma_utils.c
// @includes: nmath.h
// @depends: lgamma, log1p, nmath
// @provides: lgamma1p, log1pmx
// @all_depends_count: 12
// @all_depends: refactored, Rmath, nmath, stirlerr_cycle_free, chebyshev, cospi, fmax2, gammalims, lgammacor, log1p, gamma, lgamma
// @load_order: 64
// @local_macros: SQR

// openclport: macro hygiene pre-clean for concatenated translation units.
#ifdef SQR
# undef SQR
#endif

/*
 * Utility helpers extracted from pgamma.c during cycle refactor:
 * - log1pmx()
 * - lgamma1p()
 *
 * This keeps these reusable numerical helpers separate from pgamma()
 * distribution-flow logic so callers such as bd0.c and stirlerr.c do not
 * depend on full pgamma.c.
 */

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"

#define SQR(x) ((x) * (x))
static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxiliary in log1pmx() and lgamma1p()
 */
static double
logcf(double x, double i, double d,
      double eps /* ~ relative tolerance */)
{
    double c1 = 2 * d;
    double c2 = i + d;
    double c4 = c2 + d;
    double a1 = c2;
    double b1 = i * (c2 - i * x);
    double b2 = d * d * x;
    double a2 = c4 * c2 - b2;

#if 0
    assert (i > 0);
    assert (d >= 0);
#endif

    b2 = c4 * b1 - i * b2;

    while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
        double c3 = c2 * c2 * x;
        c2 += d;
        c4 += d;
        a1 = c4 * a2 - c3 * a1;
        b1 = c4 * b2 - c3 * b1;

        c3 = c1 * c1 * x;
        c1 += d;
        c4 += d;
        a2 = c4 * a1 - c3 * a2;
        b2 = c4 * b1 - c3 * b2;

        if (fabs(b2) > scalefactor) {
            a1 /= scalefactor;
            b1 /= scalefactor;
            a2 /= scalefactor;
            b2 /= scalefactor;
        } else if (fabs(b2) < 1 / scalefactor) {
            a1 *= scalefactor;
            b1 *= scalefactor;
            a2 *= scalefactor;
            b2 *= scalefactor;
        }
    }

    return a2 / b2;
}

/* Accurate calculation of log(1+x)-x, particularly for small x.  */
double log1pmx(double x)
{
    static const double minLog1Value = -0.79149064;

    if (x > 1 || x < minLog1Value)
        return log1p(x) - x;
    else { /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
            * log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
            * ---------------------------------------------
            * S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
           */
        double r = x / (2 + x), y = r * r;
        if (fabs(x) < 1e-2) {
            static const double two = 2;
            return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
                         two / 3) *
                            y -
                        x);
        } else {
            static const double tol_logcf = 1e-14;
            return r * (2 * y * logcf(y, 3, 2, tol_logcf) - x);
        }
    }
}

/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
double lgamma1p(double a)
{
    if (fabs(a) >= 0.5)
        return lgammafn(a + 1);

    const double eulers_const = 0.5772156649015328606065120900824024;

    /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
    const int N = 40;
    static const double coeffs[40] = {
        0.3224670334241132182362075833230126e-0, /* = (zeta(2)-1)/2 */
        0.6735230105319809513324605383715000e-1, /* = (zeta(3)-1)/3 */
        0.2058080842778454787900092413529198e-1,
        0.7385551028673985266273097291406834e-2,
        0.2890510330741523285752988298486755e-2,
        0.1192753911703260977113935692828109e-2,
        0.5096695247430424223356548135815582e-3,
        0.2231547584535793797614188036013401e-3,
        0.9945751278180853371459589003190170e-4,
        0.4492623673813314170020750240635786e-4,
        0.2050721277567069155316650397830591e-4,
        0.9439488275268395903987425104415055e-5,
        0.4374866789907487804181793223952411e-5,
        0.2039215753801366236781900709670839e-5,
        0.9551412130407419832857179772951265e-6,
        0.4492469198764566043294290331193655e-6,
        0.2120718480555466586923135901077628e-6,
        0.1004322482396809960872083050053344e-6,
        0.4769810169363980565760193417246730e-7,
        0.2271109460894316491031998116062124e-7,
        0.1083865921489695409107491757968159e-7,
        0.5183475041970046655121248647057669e-8,
        0.2483674543802478317185008663991718e-8,
        0.1192140140586091207442548202774640e-8,
        0.5731367241678862013330194857961011e-9,
        0.2759522885124233145178149692816341e-9,
        0.1330476437424448948149715720858008e-9,
        0.6422964563838100022082448087644648e-10,
        0.3104424774732227276239215783404066e-10,
        0.1502138408075414217093301048780668e-10,
        0.7275974480239079662504549924814047e-11,
        0.3527742476575915083615072228655483e-11,
        0.1711991790559617908601084114443031e-11,
        0.8315385841420284819798357793954418e-12,
        0.4042200525289440065536008957032895e-12,
        0.1966475631096616490411045679010286e-12,
        0.9573630387838555763782200936508615e-13,
        0.4664076026428374224576492565974577e-13,
        0.2273736960065972320633279596737272e-13,
        0.1109139947083452201658320007192334e-13 /* = (zeta(40+1)-1)/(40+1) */
    };

    const double c = 0.2273736845824652515226821577978691e-12; /* zeta(N+2)-1 */
    const double tol_logcf = 1e-14;

    /* Abramowitz & Stegun 6.1.33 : for |x| < 2,
     * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
     * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
     *
     * Here, another convergence acceleration trick is used to compute
     * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
     */
    double lgam = c * logcf(-a / 2, N + 2, 1, tol_logcf);
    for (int i = N - 1; i >= 0; i--)
        lgam = coeffs[i] - a * lgam;

    return (a * lgam - eulers_const) * a - log1pmx(a);
} /* lgamma1p */

// openclport: macro hygiene post-clean for concatenated translation units.
#undef SQR


// @source_type: c
// @source_origin: stirlerr_cycle_dependent.c
// @includes: nmath.h
// @depends: pgamma_utils, nmath
// @provides: stirlerr_cycle_dependent
// @all_depends_count: 13
// @all_depends: refactored, Rmath, nmath, stirlerr_cycle_free, chebyshev, cospi, fmax2, gammalims, lgammacor, log1p, gamma, lgamma, pgamma_utils
// @load_order: 73

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"

/* Branch isolated for cycle diagnostics:
 * this is the only stirlerr split file that depends on lgamma1p().
 */
attribute_hidden double stirlerr_cycle_dependent(double n)
{
    if (n <= 0.) {
        ML_WARN_return_NAN;
    }
    return lgamma1p(n) - (n + 0.5) * log(n) + n - M_LN_SQRT_2PI;
}


// @source_type: c
// @source_origin: bd0.c
// @includes: nmath.h
// @depends: pgamma_utils, nmath
// @provides: bd0, ebd0
// @all_depends_count: 13
// @all_depends: refactored, Rmath, nmath, stirlerr_cycle_free, chebyshev, cospi, fmax2, gammalims, lgammacor, log1p, gamma, lgamma, pgamma_utils
// @load_order: 74
// @local_macros: lg_x_n, ADD1

// openclport: macro hygiene pre-clean for concatenated translation units.
#ifdef lg_x_n
# undef lg_x_n
#endif
#ifdef ADD1
# undef ADD1
#endif

/*
 *  AUTHORS
 *	Catherine Loader, catherine@research.bell-labs.com, October 23, 2000. [ bd0() ]
 *	Morten Welinder, see Bugzilla PR#15628, 2014                          [ebd0() ]
 *
 *  Merge in to R (and much more):
 *
 *	Copyright (C) 2000-2025 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 *  DESCRIPTION
 *	Evaluates the "deviance part"
 *	bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =
 *		  =  x * log(x/M) + M - x
 *	where M = E[X] = n*p (or = lambda), for	  x, M > 0
 *
 *	in a manner that should be stable (with small relative error)
 *	for all x and M=np. In particular for x/np close to 1, direct
 *	evaluation fails, and evaluation is based on the Taylor series
 *	of log((1+v)/(1-v)) with v = (x-M)/(x+M) = (x-np)/(x+np).
 *
 * Martyn Plummer had the nice idea to use log1p() and Martin Maechler
 * emphasized the extra need to control cancellation.
 *
 * MP:   t := (x-M)/M  ( <==> 1+t = x/M  ==>
 *
 * bd0 = M*[ x/M * log(x/M) + 1 - (x/M) ] = M*[ (1+t)*log1p(t) + 1 - (1+t) ]
 *     = M*[ (1+t)*log1p(t) - t ] =: M * p1log1pm(t) =: M * p1l1(t)
 * MM: The above is very nice, as the "simple" p1l1() function would be useful
 *    to have available in a fast numerical stable way more generally.
 */
// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"

attribute_hidden double bd0(double x, double np)
{
    if(!R_FINITE(x) || !R_FINITE(np) || np == 0.0) ML_WARN_return_NAN;

    if (fabs(x-np) < 0.1*(x+np)) {
    	double d = x - np,
	    v = d/(x+np);
	if((d != 0.0)  && (v == 0.0)) {  // v has underflown to 0 (as  x+np = inf)
	    double
		x_ = ldexp(x, -2),
		n_ = ldexp(np,-2);
	    v = (x_ - n_)/(x_ + n_);
	}
	double s = ldexp(d, -1) * v; // was d * v
	if(fabs(ldexp(s, 1)) < DBL_MIN) return ldexp(s, 1);
	double ej = x * v; // as 2*x*v could overflow:  v > 1/2  <==> ej = 2xv > x
	v *= v; // "v = v^2"
	for (int j = 1; j < 1000; j++) { /* Taylor series; 1000: no infinite loop
					    as |v| < .1,  v^2000 is "zero" */
	    ej *= v;// = x v^(2j+1)
	    double s_ = s;
	    s += ej/((j<<1)+1);
	    if (s == s_) { /* last term was effectively 0 */
#ifdef DEBUG_bd0
		REprintf("bd0(%g, %g): T.series w/ %d terms -> bd0=%g\n", x, np, j, ldexp(s, 1));
#endif
		return ldexp(s, 1); // 2*s ; as we dropped '2 *' above
	    }
	}
	/* ---- the following should _never_ happen ------------ */
	MATHLIB_WARNING4("bd0(%g, %g): T.series failed to converge in 1000 it.; s=%g, ej/(2j+1)=%g\n",
			 x, np, s, ej/((1000<<1)+1));
    }
    /* else:  | x - np |  is not too small */
    /* NB: x/np |--> inf (overflow)  doesn't happen when called from rpois_raw() */
#define lg_x_n (R_FINITE(x/np) ? log(x/np) : (log(x) - log(np)))
    return (x > np) ? x*(lg_x_n -1.) + np
	            : x* lg_x_n + np -x;
#undef lg_x_n
}


// ebd0(): R Bugzilla PR#15628 -- proposed accuracy improvement by Morten Welinder

/*
 * A table of logs for scaling purposes.  Each value has four parts with
 * 23 bits in each.  That means each part can be multiplied by a double
 * with at most 30 bits set and not have any rounding error.  Note, that
 * the first entry is log(2).
 *
 * Entry i is associated with the value r = 0.5 + i / 256.0.  The
 * argument to log is p/q where q=1024 and p=floor(q / r + 0.5).
 * Thus r*p/q is close to 1.
 */
static const float bd0_scale[128 + 1][4] = {
	{ +0x1.62e430p-1, -0x1.05c610p-29, -0x1.950d88p-54, +0x1.d9cc02p-79 }, /* 128: log(2048/1024.) */
	{ +0x1.5ee02cp-1, -0x1.6dbe98p-25, -0x1.51e540p-50, +0x1.2bfa48p-74 }, /* 129: log(2032/1024.) */
	{ +0x1.5ad404p-1, +0x1.86b3e4p-26, +0x1.9f6534p-50, +0x1.54be04p-74 }, /* 130: log(2016/1024.) */
	{ +0x1.570124p-1, -0x1.9ed750p-25, -0x1.f37dd0p-51, +0x1.10b770p-77 }, /* 131: log(2001/1024.) */
	{ +0x1.5326e4p-1, -0x1.9b9874p-25, -0x1.378194p-49, +0x1.56feb2p-74 }, /* 132: log(1986/1024.) */
	{ +0x1.4f4528p-1, +0x1.aca70cp-28, +0x1.103e74p-53, +0x1.9c410ap-81 }, /* 133: log(1971/1024.) */
	{ +0x1.4b5bd8p-1, -0x1.6a91d8p-25, -0x1.8e43d0p-50, -0x1.afba9ep-77 }, /* 134: log(1956/1024.) */
	{ +0x1.47ae54p-1, -0x1.abb51cp-25, +0x1.19b798p-51, +0x1.45e09cp-76 }, /* 135: log(1942/1024.) */
	{ +0x1.43fa00p-1, -0x1.d06318p-25, -0x1.8858d8p-49, -0x1.1927c4p-75 }, /* 136: log(1928/1024.) */
	{ +0x1.3ffa40p-1, +0x1.1a427cp-25, +0x1.151640p-53, -0x1.4f5606p-77 }, /* 137: log(1913/1024.) */
	{ +0x1.3c7c80p-1, -0x1.19bf48p-34, +0x1.05fc94p-58, -0x1.c096fcp-82 }, /* 138: log(1900/1024.) */
	{ +0x1.38b320p-1, +0x1.6b5778p-25, +0x1.be38d0p-50, -0x1.075e96p-74 }, /* 139: log(1886/1024.) */
	{ +0x1.34e288p-1, +0x1.d9ce1cp-25, +0x1.316eb8p-49, +0x1.2d885cp-73 }, /* 140: log(1872/1024.) */
	{ +0x1.315124p-1, +0x1.c2fc60p-29, -0x1.4396fcp-53, +0x1.acf376p-78 }, /* 141: log(1859/1024.) */
	{ +0x1.2db954p-1, +0x1.720de4p-25, -0x1.d39b04p-49, -0x1.f11176p-76 }, /* 142: log(1846/1024.) */
	{ +0x1.2a1b08p-1, -0x1.562494p-25, +0x1.a7863cp-49, +0x1.85dd64p-73 }, /* 143: log(1833/1024.) */
	{ +0x1.267620p-1, +0x1.3430e0p-29, -0x1.96a958p-56, +0x1.f8e636p-82 }, /* 144: log(1820/1024.) */
	{ +0x1.23130cp-1, +0x1.7bebf4p-25, +0x1.416f1cp-52, -0x1.78dd36p-77 }, /* 145: log(1808/1024.) */
	{ +0x1.1faa34p-1, +0x1.70e128p-26, +0x1.81817cp-50, -0x1.c2179cp-76 }, /* 146: log(1796/1024.) */
	{ +0x1.1bf204p-1, +0x1.3a9620p-28, +0x1.2f94c0p-52, +0x1.9096c0p-76 }, /* 147: log(1783/1024.) */
	{ +0x1.187ce4p-1, -0x1.077870p-27, +0x1.655a80p-51, +0x1.eaafd6p-78 }, /* 148: log(1771/1024.) */
	{ +0x1.1501c0p-1, -0x1.406cacp-25, -0x1.e72290p-49, +0x1.5dd800p-73 }, /* 149: log(1759/1024.) */
	{ +0x1.11cb80p-1, +0x1.787cd0p-25, -0x1.efdc78p-51, -0x1.5380cep-77 }, /* 150: log(1748/1024.) */
	{ +0x1.0e4498p-1, +0x1.747324p-27, -0x1.024548p-51, +0x1.77a5a6p-75 }, /* 151: log(1736/1024.) */
	{ +0x1.0b036cp-1, +0x1.690c74p-25, +0x1.5d0cc4p-50, -0x1.c0e23cp-76 }, /* 152: log(1725/1024.) */
	{ +0x1.077070p-1, -0x1.a769bcp-27, +0x1.452234p-52, +0x1.6ba668p-76 }, /* 153: log(1713/1024.) */
	{ +0x1.04240cp-1, -0x1.a686acp-27, -0x1.ef46b0p-52, -0x1.5ce10cp-76 }, /* 154: log(1702/1024.) */
	{ +0x1.00d22cp-1, +0x1.fc0e10p-25, +0x1.6ee034p-50, -0x1.19a2ccp-74 }, /* 155: log(1691/1024.) */
	{ +0x1.faf588p-2, +0x1.ef1e64p-27, -0x1.26504cp-54, -0x1.b15792p-82 }, /* 156: log(1680/1024.) */
	{ +0x1.f4d87cp-2, +0x1.d7b980p-26, -0x1.a114d8p-50, +0x1.9758c6p-75 }, /* 157: log(1670/1024.) */
	{ +0x1.ee1414p-2, +0x1.2ec060p-26, +0x1.dc00fcp-52, +0x1.f8833cp-76 }, /* 158: log(1659/1024.) */
	{ +0x1.e7e32cp-2, -0x1.ac796cp-27, -0x1.a68818p-54, +0x1.235d02p-78 }, /* 159: log(1649/1024.) */
	{ +0x1.e108a0p-2, -0x1.768ba4p-28, -0x1.f050a8p-52, +0x1.00d632p-82 }, /* 160: log(1638/1024.) */
	{ +0x1.dac354p-2, -0x1.d3a6acp-30, +0x1.18734cp-57, -0x1.f97902p-83 }, /* 161: log(1628/1024.) */
	{ +0x1.d47424p-2, +0x1.7dbbacp-31, -0x1.d5ada4p-56, +0x1.56fcaap-81 }, /* 162: log(1618/1024.) */
	{ +0x1.ce1af0p-2, +0x1.70be7cp-27, +0x1.6f6fa4p-51, +0x1.7955a2p-75 }, /* 163: log(1608/1024.) */
	{ +0x1.c7b798p-2, +0x1.ec36ecp-26, -0x1.07e294p-50, -0x1.ca183cp-75 }, /* 164: log(1598/1024.) */
	{ +0x1.c1ef04p-2, +0x1.c1dfd4p-26, +0x1.888eecp-50, -0x1.fd6b86p-75 }, /* 165: log(1589/1024.) */
	{ +0x1.bb7810p-2, +0x1.478bfcp-26, +0x1.245b8cp-50, +0x1.ea9d52p-74 }, /* 166: log(1579/1024.) */
	{ +0x1.b59da0p-2, -0x1.882b08p-27, +0x1.31573cp-53, -0x1.8c249ap-77 }, /* 167: log(1570/1024.) */
	{ +0x1.af1294p-2, -0x1.b710f4p-27, +0x1.622670p-51, +0x1.128578p-76 }, /* 168: log(1560/1024.) */
	{ +0x1.a925d4p-2, -0x1.0ae750p-27, +0x1.574ed4p-51, +0x1.084996p-75 }, /* 169: log(1551/1024.) */
	{ +0x1.a33040p-2, +0x1.027d30p-29, +0x1.b9a550p-53, -0x1.b2e38ap-78 }, /* 170: log(1542/1024.) */
	{ +0x1.9d31c0p-2, -0x1.5ec12cp-26, -0x1.5245e0p-52, +0x1.2522d0p-79 }, /* 171: log(1533/1024.) */
	{ +0x1.972a34p-2, +0x1.135158p-30, +0x1.a5c09cp-56, +0x1.24b70ep-80 }, /* 172: log(1524/1024.) */
	{ +0x1.911984p-2, +0x1.0995d4p-26, +0x1.3bfb5cp-50, +0x1.2c9dd6p-75 }, /* 173: log(1515/1024.) */
	{ +0x1.8bad98p-2, -0x1.1d6144p-29, +0x1.5b9208p-53, +0x1.1ec158p-77 }, /* 174: log(1507/1024.) */
	{ +0x1.858b58p-2, -0x1.1b4678p-27, +0x1.56cab4p-53, -0x1.2fdc0cp-78 }, /* 175: log(1498/1024.) */
	{ +0x1.7f5fa0p-2, +0x1.3aaf48p-27, +0x1.461964p-51, +0x1.4ae476p-75 }, /* 176: log(1489/1024.) */
	{ +0x1.79db68p-2, -0x1.7e5054p-26, +0x1.673750p-51, -0x1.a11f7ap-76 }, /* 177: log(1481/1024.) */
	{ +0x1.744f88p-2, -0x1.cc0e18p-26, -0x1.1e9d18p-50, -0x1.6c06bcp-78 }, /* 178: log(1473/1024.) */
	{ +0x1.6e08ecp-2, -0x1.5d45e0p-26, -0x1.c73ec8p-50, +0x1.318d72p-74 }, /* 179: log(1464/1024.) */
	{ +0x1.686c80p-2, +0x1.e9b14cp-26, -0x1.13bbd4p-50, -0x1.efeb1cp-78 }, /* 180: log(1456/1024.) */
	{ +0x1.62c830p-2, -0x1.a8c70cp-27, -0x1.5a1214p-51, -0x1.bab3fcp-79 }, /* 181: log(1448/1024.) */
	{ +0x1.5d1bdcp-2, -0x1.4fec6cp-31, +0x1.423638p-56, +0x1.ee3feep-83 }, /* 182: log(1440/1024.) */
	{ +0x1.576770p-2, +0x1.7455a8p-26, -0x1.3ab654p-50, -0x1.26be4cp-75 }, /* 183: log(1432/1024.) */
	{ +0x1.5262e0p-2, -0x1.146778p-26, -0x1.b9f708p-52, -0x1.294018p-77 }, /* 184: log(1425/1024.) */
	{ +0x1.4c9f08p-2, +0x1.e152c4p-26, -0x1.dde710p-53, +0x1.fd2208p-77 }, /* 185: log(1417/1024.) */
	{ +0x1.46d2d8p-2, +0x1.c28058p-26, -0x1.936284p-50, +0x1.9fdd68p-74 }, /* 186: log(1409/1024.) */
	{ +0x1.41b940p-2, +0x1.cce0c0p-26, -0x1.1a4050p-50, +0x1.bc0376p-76 }, /* 187: log(1402/1024.) */
	{ +0x1.3bdd24p-2, +0x1.d6296cp-27, +0x1.425b48p-51, -0x1.cddb2cp-77 }, /* 188: log(1394/1024.) */
	{ +0x1.36b578p-2, -0x1.287ddcp-27, -0x1.2d0f4cp-51, +0x1.38447ep-75 }, /* 189: log(1387/1024.) */
	{ +0x1.31871cp-2, +0x1.2a8830p-27, +0x1.3eae54p-52, -0x1.898136p-77 }, /* 190: log(1380/1024.) */
	{ +0x1.2b9304p-2, -0x1.51d8b8p-28, +0x1.27694cp-52, -0x1.fd852ap-76 }, /* 191: log(1372/1024.) */
	{ +0x1.265620p-2, -0x1.d98f3cp-27, +0x1.a44338p-51, -0x1.56e85ep-78 }, /* 192: log(1365/1024.) */
	{ +0x1.211254p-2, +0x1.986160p-26, +0x1.73c5d0p-51, +0x1.4a861ep-75 }, /* 193: log(1358/1024.) */
	{ +0x1.1bc794p-2, +0x1.fa3918p-27, +0x1.879c5cp-51, +0x1.16107cp-78 }, /* 194: log(1351/1024.) */
	{ +0x1.1675ccp-2, -0x1.4545a0p-26, +0x1.c07398p-51, +0x1.f55c42p-76 }, /* 195: log(1344/1024.) */
	{ +0x1.111ce4p-2, +0x1.f72670p-37, -0x1.b84b5cp-61, +0x1.a4a4dcp-85 }, /* 196: log(1337/1024.) */
	{ +0x1.0c81d4p-2, +0x1.0c150cp-27, +0x1.218600p-51, -0x1.d17312p-76 }, /* 197: log(1331/1024.) */
	{ +0x1.071b84p-2, +0x1.fcd590p-26, +0x1.a3a2e0p-51, +0x1.fe5ef8p-76 }, /* 198: log(1324/1024.) */
	{ +0x1.01ade4p-2, -0x1.bb1844p-28, +0x1.db3cccp-52, +0x1.1f56fcp-77 }, /* 199: log(1317/1024.) */
	{ +0x1.fa01c4p-3, -0x1.12a0d0p-29, -0x1.f71fb0p-54, +0x1.e287a4p-78 }, /* 200: log(1311/1024.) */
	{ +0x1.ef0adcp-3, +0x1.7b8b28p-28, -0x1.35bce4p-52, -0x1.abc8f8p-79 }, /* 201: log(1304/1024.) */
	{ +0x1.e598ecp-3, +0x1.5a87e4p-27, -0x1.134bd0p-51, +0x1.c2cebep-76 }, /* 202: log(1298/1024.) */
	{ +0x1.da85d8p-3, -0x1.df31b0p-27, +0x1.94c16cp-57, +0x1.8fd7eap-82 }, /* 203: log(1291/1024.) */
	{ +0x1.d0fb80p-3, -0x1.bb5434p-28, -0x1.ea5640p-52, -0x1.8ceca4p-77 }, /* 204: log(1285/1024.) */
	{ +0x1.c765b8p-3, +0x1.e4d68cp-27, +0x1.5b59b4p-51, +0x1.76f6c4p-76 }, /* 205: log(1279/1024.) */
	{ +0x1.bdc46cp-3, -0x1.1cbb50p-27, +0x1.2da010p-51, +0x1.eb282cp-75 }, /* 206: log(1273/1024.) */
	{ +0x1.b27980p-3, -0x1.1b9ce0p-27, +0x1.7756f8p-52, +0x1.2ff572p-76 }, /* 207: log(1266/1024.) */
	{ +0x1.a8bed0p-3, -0x1.bbe874p-30, +0x1.85cf20p-56, +0x1.b9cf18p-80 }, /* 208: log(1260/1024.) */
	{ +0x1.9ef83cp-3, +0x1.2769a4p-27, -0x1.85bda0p-52, +0x1.8c8018p-79 }, /* 209: log(1254/1024.) */
	{ +0x1.9525a8p-3, +0x1.cf456cp-27, -0x1.7137d8p-52, -0x1.f158e8p-76 }, /* 210: log(1248/1024.) */
	{ +0x1.8b46f8p-3, +0x1.11b12cp-30, +0x1.9f2104p-54, -0x1.22836ep-78 }, /* 211: log(1242/1024.) */
	{ +0x1.83040cp-3, +0x1.2379e4p-28, +0x1.b71c70p-52, -0x1.990cdep-76 }, /* 212: log(1237/1024.) */
	{ +0x1.790ed4p-3, +0x1.dc4c68p-28, -0x1.910ac8p-52, +0x1.dd1bd6p-76 }, /* 213: log(1231/1024.) */
	{ +0x1.6f0d28p-3, +0x1.5cad68p-28, +0x1.737c94p-52, -0x1.9184bap-77 }, /* 214: log(1225/1024.) */
	{ +0x1.64fee8p-3, +0x1.04bf88p-28, +0x1.6fca28p-52, +0x1.8884a8p-76 }, /* 215: log(1219/1024.) */
	{ +0x1.5c9400p-3, +0x1.d65cb0p-29, -0x1.b2919cp-53, +0x1.b99bcep-77 }, /* 216: log(1214/1024.) */
	{ +0x1.526e60p-3, -0x1.c5e4bcp-27, -0x1.0ba380p-52, +0x1.d6e3ccp-79 }, /* 217: log(1208/1024.) */
	{ +0x1.483bccp-3, +0x1.9cdc7cp-28, -0x1.5ad8dcp-54, -0x1.392d3cp-83 }, /* 218: log(1202/1024.) */
	{ +0x1.3fb25cp-3, -0x1.a6ad74p-27, +0x1.5be6b4p-52, -0x1.4e0114p-77 }, /* 219: log(1197/1024.) */
	{ +0x1.371fc4p-3, -0x1.fe1708p-27, -0x1.78864cp-52, -0x1.27543ap-76 }, /* 220: log(1192/1024.) */
	{ +0x1.2cca10p-3, -0x1.4141b4p-28, -0x1.ef191cp-52, +0x1.00ee08p-76 }, /* 221: log(1186/1024.) */
	{ +0x1.242310p-3, +0x1.3ba510p-27, -0x1.d003c8p-51, +0x1.162640p-76 }, /* 222: log(1181/1024.) */
	{ +0x1.1b72acp-3, +0x1.52f67cp-27, -0x1.fd6fa0p-51, +0x1.1a3966p-77 }, /* 223: log(1176/1024.) */
	{ +0x1.10f8e4p-3, +0x1.129cd8p-30, +0x1.31ef30p-55, +0x1.a73e38p-79 }, /* 224: log(1170/1024.) */
	{ +0x1.08338cp-3, -0x1.005d7cp-27, -0x1.661a9cp-51, +0x1.1f138ap-79 }, /* 225: log(1165/1024.) */
	{ +0x1.fec914p-4, -0x1.c482a8p-29, -0x1.55746cp-54, +0x1.99f932p-80 }, /* 226: log(1160/1024.) */
	{ +0x1.ed1794p-4, +0x1.d06f00p-29, +0x1.75e45cp-53, -0x1.d0483ep-78 }, /* 227: log(1155/1024.) */
	{ +0x1.db5270p-4, +0x1.87d928p-32, -0x1.0f52a4p-57, +0x1.81f4a6p-84 }, /* 228: log(1150/1024.) */
	{ +0x1.c97978p-4, +0x1.af1d24p-29, -0x1.0977d0p-60, -0x1.8839d0p-84 }, /* 229: log(1145/1024.) */
	{ +0x1.b78c84p-4, -0x1.44f124p-28, -0x1.ef7bc4p-52, +0x1.9e0650p-78 }, /* 230: log(1140/1024.) */
	{ +0x1.a58b60p-4, +0x1.856464p-29, +0x1.c651d0p-55, +0x1.b06b0cp-79 }, /* 231: log(1135/1024.) */
	{ +0x1.9375e4p-4, +0x1.5595ecp-28, +0x1.dc3738p-52, +0x1.86c89ap-81 }, /* 232: log(1130/1024.) */
	{ +0x1.814be4p-4, -0x1.c073fcp-28, -0x1.371f88p-53, -0x1.5f4080p-77 }, /* 233: log(1125/1024.) */
	{ +0x1.6f0d28p-4, +0x1.5cad68p-29, +0x1.737c94p-53, -0x1.9184bap-78 }, /* 234: log(1120/1024.) */
	{ +0x1.60658cp-4, -0x1.6c8af4p-28, +0x1.d8ef74p-55, +0x1.c4f792p-80 }, /* 235: log(1116/1024.) */
	{ +0x1.4e0110p-4, +0x1.146b5cp-29, +0x1.73f7ccp-54, -0x1.d28db8p-79 }, /* 236: log(1111/1024.) */
	{ +0x1.3b8758p-4, +0x1.8b1b70p-28, -0x1.20aca4p-52, -0x1.651894p-76 }, /* 237: log(1106/1024.) */
	{ +0x1.28f834p-4, +0x1.43b6a4p-30, -0x1.452af8p-55, +0x1.976892p-80 }, /* 238: log(1101/1024.) */
	{ +0x1.1a0fbcp-4, -0x1.e4075cp-28, +0x1.1fe618p-52, +0x1.9d6dc2p-77 }, /* 239: log(1097/1024.) */
	{ +0x1.075984p-4, -0x1.4ce370p-29, -0x1.d9fc98p-53, +0x1.4ccf12p-77 }, /* 240: log(1092/1024.) */
	{ +0x1.f0a30cp-5, +0x1.162a68p-37, -0x1.e83368p-61, -0x1.d222a6p-86 }, /* 241: log(1088/1024.) */
	{ +0x1.cae730p-5, -0x1.1a8f7cp-31, -0x1.5f9014p-55, +0x1.2720c0p-79 }, /* 242: log(1083/1024.) */
	{ +0x1.ac9724p-5, -0x1.e8ee08p-29, +0x1.a7de04p-54, -0x1.9bba74p-78 }, /* 243: log(1079/1024.) */
	{ +0x1.868a84p-5, -0x1.ef8128p-30, +0x1.dc5eccp-54, -0x1.58d250p-79 }, /* 244: log(1074/1024.) */
	{ +0x1.67f950p-5, -0x1.ed684cp-30, -0x1.f060c0p-55, -0x1.b1294cp-80 }, /* 245: log(1070/1024.) */
	{ +0x1.494accp-5, +0x1.a6c890p-32, -0x1.c3ad48p-56, -0x1.6dc66cp-84 }, /* 246: log(1066/1024.) */
	{ +0x1.22c71cp-5, -0x1.8abe2cp-32, -0x1.7e7078p-56, -0x1.ddc3dcp-86 }, /* 247: log(1061/1024.) */
	{ +0x1.03d5d8p-5, +0x1.79cfbcp-31, -0x1.da7c4cp-58, +0x1.4e7582p-83 }, /* 248: log(1057/1024.) */
	{ +0x1.c98d18p-6, +0x1.a01904p-31, -0x1.854164p-55, +0x1.883c36p-79 }, /* 249: log(1053/1024.) */
	{ +0x1.8b31fcp-6, -0x1.356500p-30, +0x1.c3ab48p-55, +0x1.b69bdap-80 }, /* 250: log(1049/1024.) */
	{ +0x1.3cea44p-6, +0x1.a352bcp-33, -0x1.8865acp-57, -0x1.48159cp-81 }, /* 251: log(1044/1024.) */
	{ +0x1.fc0a8cp-7, -0x1.e07f84p-32, +0x1.e7cf6cp-58, +0x1.3a69c0p-82 }, /* 252: log(1040/1024.) */
	{ +0x1.7dc474p-7, +0x1.f810a8p-31, -0x1.245b5cp-56, -0x1.a1f4f8p-80 }, /* 253: log(1036/1024.) */
	{ +0x1.fe02a8p-8, -0x1.4ef988p-32, +0x1.1f86ecp-57, +0x1.20723cp-81 }, /* 254: log(1032/1024.) */
	{ +0x1.ff00acp-9, -0x1.d4ef44p-33, +0x1.2821acp-63, +0x1.5a6d32p-87 }, /* 255: log(1028/1024.) */
	{ 0, 0, 0, 0 } /* log(1024/1024) = log(1) = 0 */
};


/*
 * Compute x * log (x / M) + (M - x)
 * aka -x * log1pmx ((M - x) / x)
 *
 * Deliver the result back in two parts, *yh and *yl.
 */
attribute_hidden void ebd0(double x, double M, double *yh, double *yl)
{
	const int Sb = 10;
	const double S = 1u << Sb; // = 2^10 = 1024
	const int N = 128; // == ? == G_N_ELEMENTS(bd0_scale) - 1; <<<< FIXME:

	*yl = *yh = 0;

	if (x == M) return;
	if (x == 0) { *yh = M;         return; }
	if (M == 0) { *yh = ML_POSINF; return; }

	if (M/x == ML_POSINF) { *yh = M; return; }//  as when (x == 0)

	int e;
	// NB: M/x overflow handled above; underflow should be handled by fg = Inf
	double r = frexp (M / x, &e); // => r in  [0.5, 1) and 'e' (int) such that  M/x = r * 2^e

	// prevent later overflow
	if (M_LN2 * ((double) -e)  > 1. + DBL_MAX / x) { *yh = ML_POSINF; return; }

	int i = (int) floor ((r - 0.5) * (2 * N) + 0.5);
	// now,  0 <= i <= N
	double f = floor (S / (0.5 + i / (2.0 * N)) + 0.5);
	double fg = ldexp (f, -(e + Sb)); // ldexp(f, E) := f * 2^E
#ifdef DEBUG_bd0
	REprintf("ebd0(x=%g, M=%g): M/x = (r=%.15g) * 2^(e=%d); i=%d,\n  f=%g, fg=f*2^-(e+%d)=%g\n",
		 x, M, r,e, i, f, Sb, fg);
	if (fg == ML_POSINF) {
	    REprintf(" --> fg = +Inf --> return( +Inf )\n");
	    *yh = fg; return;
	}
	REprintf("     bd0_sc[0][0..3]= ("); for(int j=0; j < 4; j++) REprintf("%g ", bd0_scale[0][j]); REprintf(")\n");
	REprintf("i -> bd0_sc[i][0..3]= ("); for(int j=0; j < 4; j++) REprintf("%g ", bd0_scale[i][j]); REprintf(")\n");
	REprintf( "  small(?)  (M*fg-x)/x = (M*fg)/x - 1 = %.16g\n", (M*fg-x)/x);
#else
	if (fg == ML_POSINF) {
	    *yh = fg; return;
	}
#endif
	/* We now have (M * fg / x) close to 1.  */

	/*
	 * We need to compute this:
	 * (x/M)^x * exp(M-x) =
	 * (M/x)^-x * exp(M-x) =
	 * (M*fg/x)^-x * (fg)^x * exp(M-x) =
	 * (M*fg/x)^-x * (fg)^x * exp(M*fg-x) * exp(M-M*fg)
	 *
	 * In log terms:
	 * log((x/M)^x * exp(M-x)) =
	 * log((M*fg/x)^-x * (fg)^x * exp(M*fg-x) * exp(M-M*fg)) =
	 * log((M*fg/x)^-x * exp(M*fg-x)) + x*log(fg) + (M-M*fg) =
	 * -x*log1pmx((M*fg-x)/x) + x*log(fg) + M - M*fg =
	 *
	 * Note, that fg has at most 10 bits.  If M and x are suitably
	 * "nice" -- such as being integers or half-integers -- then
	 * we can compute M*fg as well as x * bd0_scale[.][.] without
	 * rounding errors.
	 */

#define ADD1(d_) do {				\
   volatile double d = (d_);			\
	    double d1 = floor (d + 0.5);	\
	    double d2 = d - d1;/* in [-.5,.5) */ \
	    *yh += d1;				\
	    *yl += d2;				\
	} while(0)

#ifdef DEBUG_bd0
	{
	    double log1__ = log1pmx((M * fg - x) / x),
		xl = -x * log1__;
	    REprintf(" 1a. before adding  -x * log1pmx(.) = -x * %g = %g\n", log1__, xl);
	    ADD1(xl);
	    REprintf(" 1. after A.(-x*l..):       yl,yh = (%13g, %13g); yl+yh= %g\n",
		     *yl, *yh, (*yl)+(*yh));
	}
        if(fg == 1) {
            REprintf("___ fg = 1 ___ skipping further steps\n");
            return;
        }
	// else  [ fg != 1 ]
	REprintf(" 2:  A(x*b[i,j]) and A(-x*e*b[0,j]), j=1:4:\n");
	for (int j = 0; j < 4; j++) {
 	    ADD1( x * bd0_scale[i][j]);     // handles  x*log(fg*2^e)
	    REprintf(" j=%d: (%13g, %13g);", j, *yl, *yh);
	    ADD1(-x * bd0_scale[0][j] * e); // handles  x*log(1/ 2^e)
	    REprintf(" (%13g, %13g); yl+yh= %g\n", *yl, *yh, (*yl)+(*yh));
            if(!R_FINITE(*yh)) {
                REprintf(" non-finite yh --> return((yh=Inf, yl=0))\n");
		*yh = ML_POSINF; *yl = 0; return;
            }
	}
#else
	ADD1(-x * log1pmx ((M * fg - x) / x));
        if(fg == 1) return;
	// else (fg != 1) :
	for (int j = 0; j < 4; j++) {
	    ADD1( x * bd0_scale[i][j]);     // handles  x*log(fg*2^e)
	    ADD1(-x * bd0_scale[0][j] * e); // handles  x*log(1/ 2^e)
	    //                        ^^^ at end prevents overflow in  ebd0(1e307, 1e300)
            if(!R_FINITE(*yh)) { *yh = ML_POSINF; *yl = 0; return; }
	}
#endif

	ADD1(M);
#ifdef DEBUG_bd0
	REprintf(" 3. after ADD1(M):            yl,yh = (%13g, %13g); yl+yh= %g\n", *yl, *yh, (*yl)+(*yh));
#endif
	ADD1(-M * fg);
#ifdef DEBUG_bd0
	REprintf(" 4. after ADD1(- M*fg):       yl,yh = (%13g, %13g); yl+yh= %g\n\n", *yl, *yh, (*yl)+(*yh));
#endif
}

#undef ADD1

// openclport: macro hygiene post-clean for concatenated translation units.
#undef lg_x_n
#undef ADD1


// @source_type: c
// @source_origin: stirlerr.c
// @includes: nmath.h, refactored.h
// @depends: stirlerr_cycle_dependent, stirlerr_cycle_free, nmath, refactored
// @provides: stirlerr
// @all_depends_count: 14
// @all_depends: refactored, Rmath, nmath, stirlerr_cycle_free, chebyshev, cospi, fmax2, gammalims, lgammacor, log1p, gamma, lgamma, pgamma_utils, stirlerr_cycle_dependent
// @load_order: 81

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"
// openclport-disabled-include: #include "refactored.h"

/* Wrapper split for cycle refactor:
 * - stirlerr_cycle_free handles n >= 1 without lgamma1p
 * - stirlerr_cycle_dependent handles n < 1 using lgamma1p
 */
attribute_hidden double stirlerr(double n)
{
    if (n < 1.) {
        return stirlerr_cycle_dependent(n);
    }
    return stirlerr_cycle_free(n);
}


// @source_type: c
// @source_origin: dbinom.c
// @includes: nmath.h, dpq.h
// @depends: bd0, log1p, stirlerr, nmath, dpq
// @provides: dbinom, dbinom_raw, pow1p
// @all_depends_count: 17
// @all_depends: dpq, refactored, Rmath, nmath, stirlerr_cycle_free, chebyshev, cospi, fmax2, gammalims, lgammacor, log1p, gamma, lgamma, pgamma_utils, stirlerr_cycle_dependent, bd0, stirlerr
// @load_order: 85

/*
 * AUTHOR
 *   Catherine Loader, catherine@research.bell-labs.com.
 *   October 23, 2000.
 *
 *  Merge in to R and further tweaks :
 *  notably using log1p() and pow1p(), thanks to Morten Welinder, PR#18642
 *
 *	Copyright (C) 2000-2025 The R Core Team
 *	Copyright (C) 2008 The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 * DESCRIPTION
 *
 *   To compute the binomial probability, call dbinom(x,n,p).
 *   This checks for argument validity, and calls dbinom_raw().
 *
 *   dbinom_raw() does the actual computation; note this is called by
 *   other functions in addition to dbinom().
 *     (1) dbinom_raw() has both p and q arguments, when one may be represented
 *         more accurately than the other (in particular, in df()).
 *     (2) dbinom_raw() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *         -- but is not the case at all when called e.g., from df() or dbeta() !
 *     (3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's.
 *         Do this in the calling function.
 */

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"
// openclport-disabled-include: #include "dpq.h"

/* Compute  (1+x)^y  accurately also for |x| << 1  */
double pow1p(double x, double y)
{
    if(isnan(y))
	return (x == 0) ? 1. : y; // (0+1)^NaN := 1  by standards
    if(0 <= y && y == trunc(y) && y <= 4.) {
	switch((int)y) {
	case 0: return 1;
	case 1: return x + 1.;
	case 2: return x*(x + 2.) + 1.;
	case 3: return x*(x*(x + 3.) + 3.) + 1.;
	case 4: return x*(x*(x*(x + 4.) + 6.) + 4.) + 1.;
	}
    }
    /* naive algorithm in two cases: (1) when 1+x is exact (compiler should not over-optimize !),
     * and (2) when |x| > 1/2 and we have no better algorithm.
     */
    volatile double xp1 = x + 1., x_ = xp1 - 1.; // compiler should *not* optimize these
    if (x_ == x || fabs(x) > 0.5 || isnan(x)) {
	return pow(xp1, y);
    } else { /* not perfect, e.g., for small |x|, non-huge y, use
	    binom expansion 1 + y*x + y(y-1)/2 x^2 + .. */
	return exp(y * log1p(x));
    }
}

double dbinom_raw(double x, double n, double p, double q, int give_log)
{
    if (p == 0) return((x == 0) ? R_D__1 : R_D__0);
    if (q == 0) return((x == n) ? R_D__1 : R_D__0);

    // NB: The smaller of p and q is the most accurate
    if (x == 0) {
	if(n == 0) return R_D__1;
	if (p > q)
	    return give_log ? n * log(q)    : pow(q, n);
	else // 0 < p <= 1/2
	    return give_log ? n * log1p(-p) : pow1p(-p, n);
    }
    if (x == n) { // r = p^x = p^n  -- accurately
	if (p > q)
	    return give_log ? n * log1p(-q) : pow1p(-q, n);
	else
	    return give_log ? n * log (p)   : pow (p, n);
    }
    if (x < 0 || x > n) return( R_D__0 );

    if(!R_FINITE(n)) {
	if(R_FINITE(x)) return( R_D__0 ); /* finite x << n = Inf */
	else n = DBL_MAX; // helps ? extreme dnbinom() cases
    }

// TODO?  Improve accuracy in these cases:
#ifdef _NO_LOG_DBINOM_
    if(!give_log) { // more accurate *not* going via log when result is much much smaller than 1
	if (x <= M || n-x <= M) { /* use "recursive" direct formula with
				     k := min(x, n-x) multiplications */
	}
    }
#endif

    /* n*p or n*q can underflow to zero if n and p or q are small.  This
       used to occur in dbeta, and gives NaN as from R 2.3.0.  */
    double lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) - bd0(x,n*p) - bd0(n-x,n*q);

    /* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
    /* Upto R 2.7.1:
     * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
     * -- following is much better for  x << n : */
    double lf = M_LN_2PI + log(x) + log1p(- x/n);

    return R_D_exp(lc - 0.5*lf);
}

double dbinom(double x, double n, double p, int give_log)
{
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x) || ISNAN(n) || ISNAN(p)) return x + n + p;
#endif

    if (p < 0 || p > 1 || R_D_negInonint(n))
	ML_WARN_return_NAN;
    R_D_nonint_check(x);
    if (x < 0 || !R_FINITE(x)) return R_D__0;

    n = R_forceint(n);
    x = R_forceint(x);

    return dbinom_raw(x, n, p, 1-p, give_log);
}



// @library_deps: nmath
// @calls_nmath: dbinom_raw
// @depends_nmath: dbinom
// @calls_opencl_builtin: (none)
// @all_depends_nmath_count: 18
// @all_depends_nmath: dpq, refactored, Rmath, nmath, stirlerr_cycle_free, chebyshev, cospi, fmax2, gammalims, lgammacor, log1p, gamma, lgamma, pgamma_utils, stirlerr_cycle_dependent, bd0, stirlerr, dbinom

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable

#define MAX_L2 64

static inline double nll_binomial_glmb_ocl(double y_prop, double wt, double mean_p_raw) {
    int trials  = (int)round(wt);
    int success = (int)round(y_prop * wt);
    double p = fmin(1.0, fmax(0.0, mean_p_raw));
    double q = 1.0 - p;
    double logpmf = dbinom_raw((double)success, (double)trials, p, q, 1);
    return -logpmf;
}


__kernel void f2_f3_binomial_cloglog(
    __global const double* X,
    __global const double* B,
    __global const double* mu,
    __global const double* P,
    __global const double* alpha,
    __global const double* y,
    __global const double* wt,
    __global double*       qf,
    __global double*       grad,
    const int l1,
    const int l2,
    const int m1
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ℓ = 0; ℓ < l2; ++ℓ) {
            acc += P[k*l2 + ℓ] * (B[j*l2 + ℓ] - mu[ℓ]);
        }
        tmp[k] = acc;
    }

    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        qsum += (B[j*l2 + k] - mu[k]) * tmp[k];
    }
    double res_acc = 0.5 * qsum;

    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        g_loc[k] = tmp[k];
    }

    for (int i = 0; i < l1; ++i) {
        double dot = alpha[i];
        for (int k = 0; k < l2; ++k) {
            dot += X[k*l1 + i] * B[j*l2 + k];
        }

        double p1 = -expm1(-exp(dot));
        double p2 = exp(-exp(dot));
        double atemp = exp(dot - exp(dot));

        res_acc += nll_binomial_glmb_ocl(y[i], wt[i], p1);

        double resid = ((y[i] * atemp / p1) - ((1.0 - y[i]) * atemp / p2)) * wt[i];
        for (int k = 0; k < l2; ++k) {
            g_loc[k] -= X[k*l1 + i] * resid;
        }
    }

    qf[j] = res_acc;

    for (int k = 0; k < l2; ++k) {
        grad[k * m1 + j] = g_loc[k];
    }
}

