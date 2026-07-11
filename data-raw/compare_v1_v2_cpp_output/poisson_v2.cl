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



// @library_deps: (none)
// @calls_nmath: (none)
// @calls_opencl_builtin: lgamma

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed


__kernel void f2_f3_poisson(
    __global const double* X,      // design matrix,   l1 × l2, col-major
    __global const double* B,      // grid points,     m1 × l2, row-major per grid
    __global const double* mu,     // prior mean,      length = l2
    __global const double* P,      // prior precision, l2 × l2, row-major
    __global const double* alpha,  // offsets,         length = l1
    __global const double* y,      // counts,          length = l1
    __global const double* wt,     // weights,         length = l1
    __global double*       qf,     // out: neg-log-post, length = m1
//    __global double*       xb,     // out: μ = exp(η),    size = m1 × l1
    __global double*       grad,   // out: ∂(neg-log-post)/∂B, size = l2 × m1 (col-major)
    const int l1,
    const int l2,
    const int m1
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    // 1) Prior term: tmp = P * (B_j - mu)
    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ell = 0; ell < l2; ++ell) {
            acc += P[k*l2 + ell] * (B[j*l2 + ell] - mu[ell]);
        }
        tmp[k] = acc;
    }

    // 2) quadratic form qf[j] = 0.5 * (B_j - mu)' * tmp
    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        qsum += (B[j*l2 + k] - mu[k]) * tmp[k];
    }
    qf[j] = 0.5 * qsum;
 

    // 3) Initialize gradient with prior part
    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        g_loc[k] = tmp[k];
    }

    // 4) Data term: Poisson log-link
    int base = j * l1;
    for (int i = 0; i < l1; ++i) {
        // linear predictor η = α[i] + X_i·B_j
        double dot = alpha[i];
        for (int k = 0; k < l2; ++k) {
            dot += X[k*l1 + i] * B[j*l2 + k];
        }

        // fitted mean
        double mui = exp(dot);
//        xb[base + i] = mui;

        // negate log-likelihood contribution

        // Non-integer values requires replacement of dpois with version using lgamma function

//        double logp = dpois(y[i], mui, 1);
          double logp = -mui + y[i] * log(mui) - lgamma(y[i] + 1.0);
        
        qf[j] -= wt[i] * logp;

        // gradient contribution: -(y - μ) * X * wt
        double resid = (y[i] - mui) * wt[i];
        for (int k = 0; k < l2; ++k) {
            g_loc[k] -= X[k*l1 + i] * resid;
        }
    }



    // 5) Write back
  
    for (int k = 0; k < l2; ++k) {
        grad[k*m1 + j] = g_loc[k];
    }
}
