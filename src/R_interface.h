/**
 * @file R_interface.h
 * @brief Centralized, static‑cached access to R functions used throughout the
 *        C++ backend of glmbayes, providing a clean and explicit C++ → R
 *        boundary.
 *
 * @namespace glmbayes_R
 * @brief Cached R‑function accessors for non‑family‑specific utilities such as
 *        coercion, grid construction, optimization, system queries, and
 *        distribution helpers.
 *
 * @section Context
 *   This header defines a stable interface for calling R functions from C++.
 *   Each accessor returns a *static‑cached* `Rcpp::Function`, ensuring:
 *     - no repeated dynamic lookups,
 *     - no hidden R‑level dependencies scattered across source files,
 *     - a single, explicit boundary between C++ and R.
 *
 *   Only global, non‑family‑specific R functions belong here.
 *   Family‑dependent functions (e.g., f2, f3, glmbfamfunc) must remain runtime
 *   arguments and must NOT be wrapped here, as they vary by model and are part
 *   of the statistical interface rather than the infrastructure layer.
 *
 * @section ImplementedIn
 *   These declarations are implemented inline within this header.
 *
 * @section UsedBy
 *   These accessors are consumed by:
 *     - envelope construction routines (EnvelopeBuild, EnvelopeSort, EnvelopeOpt),
 *     - simulation and sampling routines (Normal, Normal–Gamma, GLM samplers),
 *     - optimization and model‑fitting components,
 *     - R‑facing wrappers that require consistent access to R utilities.
 *
 * @section Responsibilities
 *   Provides:
 *     - time and formatting helpers (`format`, `Sys.time`),
 *     - coercion utilities (`as.matrix`, `as.vector`, `as.numeric`),
 *     - grid and envelope helpers (`expand.grid`, `EnvelopeOpt`, `EnvelopeSort`),
 *     - interactive utilities (`interactive`, `readline`),
 *     - distribution helpers (`qgamma`, `rgamma_ct`, `runif`),
 *     - optimization and model‑fitting helpers (`optim`, `try`, `lm.fit`,
 *       `lm.wfit`, `gaussian`, `rNormal_reg.wfit`),
 *     - system utilities (`system.file`).
 *
 *   This module:
 *     - ensures consistent, centralized access to R functions,
 *     - avoids repeated lookup overhead,
 *     - keeps the C++ codebase modular, explicit, and maintainable.
 */


#ifndef R_INTERFACE_H
#define R_INTERFACE_H

#include <Rcpp.h>
#include "package_ns.h"

// -----------------------------------------------------------------------------
//  R_interface.h
//
//  Centralized, static‑cached accessors for R functions used across the C++
//  codebase.  This header provides a clean, explicit C++ → R boundary and
//  eliminates repeated dynamic lookups scattered throughout the source files.
//
//  Only *global, non‑family‑specific* R functions belong here.
//  Family‑dependent functions (f2, f3, glmbfamfunc, etc.) must remain
//  runtime arguments and should NOT be wrapped here.
// -----------------------------------------------------------------------------

namespace glmbayes_R {

// Registered-namespace environment for package-local R callbacks (EnvelopeOpt,
// EnvelopeSort, glmbfamfunc, rNormal_reg.wfit, rgamma_ct).  Base-R utilities
// below still resolve from the search path.
inline Rcpp::Environment pkg_env() {
  static Rcpp::Environment env = Rcpp::Environment::namespace_env(GLMBAYES_R_NS);
  return env;
}

// -----------------------------------------------------------------------------
//  Time / formatting utilities
// -----------------------------------------------------------------------------

inline Rcpp::Function r_format() {
  static Rcpp::Function fn("format");
  return fn;
}

inline Rcpp::Function r_sys_time() {
  static Rcpp::Function fn("Sys.time");
  return fn;
}

// inline std::string timestamp() {
//   return Rcpp::as<std::string>( r_format()( r_sys_time()() ) );
// }


// -----------------------------------------------------------------------------
//  Basic coercion helpers
// -----------------------------------------------------------------------------

inline Rcpp::Function r_as_matrix() {
  static Rcpp::Function fn("as.matrix");
  return fn;
}

inline Rcpp::Function r_as_vector() {
  static Rcpp::Function fn("as.vector");
  return fn;
}

inline Rcpp::Function r_as_numeric() {
  static Rcpp::Function fn("as.numeric");
  return fn;
}


// -----------------------------------------------------------------------------
//  Grid / envelope helpers
// -----------------------------------------------------------------------------

inline Rcpp::Function r_expand_grid() {
  static Rcpp::Function fn("expand.grid");
  return fn;
}

inline Rcpp::Function r_envelope_opt() {
  static Rcpp::Function fn = pkg_env()["EnvelopeOpt"];
  return fn;
}

inline Rcpp::Function r_envelope_sort() {
  static Rcpp::Function fn = pkg_env()["EnvelopeSort"];
  return fn;
}

inline Rcpp::Function r_glmbfamfunc() {
  static Rcpp::Function fn = pkg_env()["glmbfamfunc"];
  return fn;
}


// -----------------------------------------------------------------------------
//  Interactive / readline utilities
// -----------------------------------------------------------------------------

inline Rcpp::Function r_interactive() {
  static Rcpp::Function fn("interactive");
  return fn;
}

inline Rcpp::Function r_readline() {
  static Rcpp::Function fn("readline");
  return fn;
}


// -----------------------------------------------------------------------------
//  Distribution helpers (Gamma/Gaussian samplers)
// -----------------------------------------------------------------------------

inline Rcpp::Function r_qgamma() {
  static Rcpp::Function fn("qgamma");
  return fn;
}

inline Rcpp::Function r_rgamma_ct() {
  static Rcpp::Function fn = pkg_env()["rgamma_ct"];
  return fn;
}

inline Rcpp::Function r_runif() {
  static Rcpp::Function fn("runif");
  return fn;
}


// -----------------------------------------------------------------------------
//  Optimization / model‑fitting helpers
// -----------------------------------------------------------------------------

inline Rcpp::Function r_optim() {
  static Rcpp::Function fn("optim");
  return fn;
}

inline Rcpp::Function r_try() {
  static Rcpp::Function fn("try");
  return fn;
}

inline Rcpp::Function r_lm_fit() {
  static Rcpp::Function fn("lm.fit");
  return fn;
}

inline Rcpp::Function r_lm_wfit() {
  static Rcpp::Function fn("lm.wfit");
  return fn;
}

inline Rcpp::Function r_gaussian() {
  static Rcpp::Function fn("gaussian");
  return fn;
}

inline Rcpp::Function r_rNormal_reg_wfit() {
  static Rcpp::Function fn = pkg_env()["rNormal_reg.wfit"];
  return fn;
}


// -----------------------------------------------------------------------------
//  System utilities
// -----------------------------------------------------------------------------

inline Rcpp::Function r_system_file() {
  static Rcpp::Function fn("system.file");
  return fn;
}

} // namespace glmbayes_R


#endif
