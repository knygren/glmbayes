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
// @source_origin: dnorm.c
// @includes: nmath.h, dpq.h
// @depends: nmath, dpq
// @provides: dnorm, dnorm4
// @all_depends_count: 3
// @all_depends: dpq, Rmath, nmath
// @load_order: 31

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2014 The R Core Team
 *  Copyright (C) 2003	    The R Foundation
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
 *	double dnorm4(double x, double mu, double sigma, int give_log)
 *	      {dnorm (..) is synonymous and preferred inside R}
 *
 *  DESCRIPTION
 *
 *	Compute the density of the normal distribution.
 */

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include "nmath.h"
// openclport-disabled-include: #include "dpq.h"

double dnorm4(double x, double mu, double sigma, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
	return x + mu + sigma;
#endif
    if (sigma < 0) ML_WARN_return_NAN;
    if(!R_FINITE(sigma)) return R_D__0;
    if(!R_FINITE(x) && mu == x) return ML_NAN;/* x-mu is NaN */
    if (sigma == 0) 
	return (x == mu) ? ML_POSINF : R_D__0;
    x = (x - mu) / sigma;

    if(!R_FINITE(x)) return R_D__0;

    x = fabs (x);
    if (x >= 2 * sqrt(DBL_MAX)) return R_D__0;
    if (give_log)
	return -(M_LN_SQRT_2PI + 0.5 * x * x + log(sigma));
    //  M_1_SQRT_2PI = 1 / sqrt(2 * pi)
#ifdef MATHLIB_FAST_dnorm
    // and for R <= 3.0.x and R-devel upto 2014-01-01:
    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;
#else
    // more accurate, less fast :
    if (x < 5)    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;

    /* ELSE:

     * x*x  may lose upto about two digits accuracy for "large" x
     * Morten Welinder's proposal for PR#15620
     * https://bugs.r-project.org/show_bug.cgi?id=15620

     * -- 1 --  No hoop jumping when we underflow to zero anyway:

     *  -x^2/2 <         log(2)*.Machine$double.min.exp  <==>
     *     x   > sqrt(-2*log(2)*.Machine$double.min.exp) =IEEE= 37.64031
     * but "thanks" to denormalized numbers, underflow happens a bit later,
     *  effective.D.MIN.EXP <- with(.Machine, double.min.exp + double.ulp.digits)
     * for IEEE, DBL_MIN_EXP is -1022 but "effective" is -1074
     * ==> boundary = sqrt(-2*log(2)*(.Machine$double.min.exp + .Machine$double.ulp.digits))
     *              =IEEE=  38.58601
     * [on one x86_64 platform, effective boundary a bit lower: 38.56804]
     */
    if (x > sqrt(-2*M_LN2*(DBL_MIN_EXP + 1-DBL_MANT_DIG))) return 0.;

    /* Now, to get full accuracy, split x into two parts,
     *  x = x1+x2, such that |x2| <= 2^-16.
     * Assuming that we are using IEEE doubles, that means that
     * x1*x1 is error free for x<1024 (but we have x < 38.6 anyway).

     * If we do not have IEEE this is still an improvement over the naive formula.
     */
    double x1 = //  R_forceint(x * 65536) / 65536 =
	ldexp( R_forceint(ldexp(x, 16)), -16);
    double x2 = x - x1;
    return M_1_SQRT_2PI / sigma *
	(exp(-0.5 * x1 * x1) * exp( (-0.5 * x2 - x1) * x2 ) );
#endif
}



// @library_deps: nmath
// @calls_nmath: dnorm4
// @depends_nmath: dnorm
// @calls_opencl_builtin: (none)
// @all_depends_nmath_count: 4
// @all_depends_nmath: dpq, Rmath, nmath, dnorm

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed

__kernel void f2_f3_gaussian(
    __global const double* X,    // l1 x l2, column-major: X[k*l1 + i]
    __global const double* B,    // m1 x l2, row-major:    B[j*l2 + k]
    __global const double* mu,   // l2
    __global const double* P,    // l2 x l2, row-major
    __global const double* alpha,// l1
    __global const double* y,    // l1
    __global const double* wt,   // l1 (precision)
    __global double*       qf,   // m1
    __global double*       grad, // m1 x l2, column-major: grad[k*m1 + j]
    const int l1,
    const int l2,
    const int m1
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    // tmp = P * (B_j - mu)
    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ell = 0; ell < l2; ++ell) {
            acc += P[k*l2 + ell] * (B[j*l2 + ell] - mu[ell]);
        }
        tmp[k] = acc;
    }

    // objective: 0.5*(B_j - mu)'*tmp
    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        double dk = B[j*l2 + k] - mu[k];
        qsum += dk * tmp[k];
    }
    double res_acc = 0.5 * qsum;

    // gradient starts with prior term
    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) g_loc[k] = tmp[k];

    // data term
    for (int i = 0; i < l1; ++i) {
        // dot = alpha[i] + X[i,*] %* B_j
        double dot = alpha[i];
        for (int k = 0; k < l2; ++k) {
            dot += X[k*l1 + i] * B[j*l2 + k];
        }

        // objective: -log dnorm(y | mean=dot, sd=1/sqrt(wt))
        double wi   = wt[i];
        double sd_i = 1.0 / sqrt(wi);
        double ll   = dnorm4(y[i], dot, sd_i, 1);
        res_acc    += -ll;

        // gradient contribution: Xᵀ * (wt * (dot - y))
        double resid = wi * (dot - y[i]);
        for (int k = 0; k < l2; ++k) {
            g_loc[k] += X[k*l1 + i] * resid;
        }
    }

    // write back
    qf[j] = res_acc;
    for (int k = 0; k < l2; ++k) {
        grad[k*m1 + j] = g_loc[k];
    }
}

