/**
 * @file openclPort.h
 * @brief Public OpenCL interface for glmbayes, including kernel loading,
 *        device discovery, capability probing, and Rcpp-to-std::vector
 *        conversion helpers.
 *
 * @namespace openclPort
 * @brief Lightweight OpenCL utility layer providing kernel management and
 *        device‑level information for optional GPU acceleration.
 *
 * @section ImplementedIn
 *   These declarations are implemented in:
 *     - OpenCL_helper.cpp
 *     - opencl_detect.cpp
 *     - kernel_loader.cpp
 *     - (optional) additional OpenCL backend files guarded by USE_OPENCL
 *
 * @section UsedBy
 *   These functions are consumed by:
 *     - Envelope construction routines (EnvelopeBuild, EnvelopeEval,
 *       EnvelopeDispersionBuild) when OpenCL acceleration is enabled
 *     - R wrappers that expose GPU availability and kernel loading to users
 *
 * @section Responsibilities
 *   Provides:
 *     - Rcpp → std::vector conversion utilities for kernel argument buffers
 *     - GPU/device enumeration and capability checks (gpu_names, has_opencl)
 *     - Kernel source and library loading from inst/cl/ directories
 *       (generic loaders; full likelihood-subgradient program assembly is in
 *       glmbayes::opencl::load_likelihood_subgradient_program — see inst/cl/README.md)
 *     - Conditional OpenCL configuration and build‑option generation
 *
 *   This module:
 *     - is optional and only active when compiled with USE_OPENCL,
 *     - isolates all OpenCL‑specific logic from the statistical code,
 *     - ensures safe fallback to CPU execution when no GPU is available.
 */


#ifndef OPENCLPORT_H
#define OPENCLPORT_H

#include <RcppArmadillo.h>
#include <string>
#include <vector>

#ifdef USE_OPENCL

// Ensure OpenCL types are available
#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>
#include <string>
#endif 

#ifdef __linux__
#include <stdio.h>
#include <stdlib.h>
#endif

using namespace Rcpp;

// Dependencies:

// 1) OpenCL_helper.cpp
// 2) 



//
// -----------------------------------------------------------------------------
// openclPort: Public API for OpenCL kernel loading, device utilities,
//             and Rcpp → std::vector conversion helpers.
// -----------------------------------------------------------------------------
// Everything a user needs to write OpenCL-enabled wrappers lives here.
// -----------------------------------------------------------------------------
namespace openclPort {

// -------------------------------------------------------------------------
// Rcpp → std::vector conversion utilities
// -------------------------------------------------------------------------
std::vector<double> flattenMatrix(const Rcpp::NumericMatrix& mat);
std::vector<double> copyVector(const Rcpp::NumericVector& vec);

// -------------------------------------------------------------------------
// Device / OpenCL utilities
// -------------------------------------------------------------------------
Rcpp::CharacterVector gpu_names();

// Internal-only GPU detection (used by envelope scaling)
int detect_num_gpus_internal();


// -------------------------------------------------------------------------
// R-facing wrappers for OpenCL build probes
// -------------------------------------------------------------------------

bool has_opencl();
int get_opencl_core_count();


// -------------------------------------------------------------------------
// Conditional declarations: only available when USE_OPENCL is defined
// -------------------------------------------------------------------------
#ifdef USE_OPENCL

// Ensure OpenCL types are available
#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>
#include <string>


// Load a single .cl kernel file from inst/cl/<relative_path>
std::string load_kernel_source(
    const std::string& relative_path,
    const std::string& package = "glmbayes"
);

// Load and concatenate all .cl files in a subdirectory (inst/cl/<subdir>/)
std::string load_kernel_library(
    const std::string& subdir,
    const std::string& package = "glmbayes",
    bool verbose = false
);

// Minimal library subset from a kernel annotation; kernel and library may
// live in different installed packages (paths relative to inst/cl/).
std::string load_library_for_kernel(
    const std::string& kernel_relative_path,
    const std::string& library_subdir,
    const std::string& package,
    const std::string& depends_tag);

std::string load_library_for_kernel_cross_package(
    const std::string& kernel_relative_path,
    const std::string& kernel_package,
    const std::string& library_subdir,
    const std::string& library_package,
    const std::string& depends_tag);

enum class ProgramPreloadKind { file, library };

struct ProgramPreloadEntry {
  int rank;
  ProgramPreloadKind kind;
  std::string rel_path;
};

std::vector<ProgramPreloadEntry> read_program_preload_manifest(
    const std::string& manifest_relative_path,
    const std::string& source_package);

std::string load_program_preload_entries(
    const std::vector<ProgramPreloadEntry>& entries,
    const std::string& source_package,
    bool verbose = false);

std::string load_program_preload(
    const std::string& manifest_relative_path = "program_preload_manifest.tsv",
    const std::string& source_package = "nmathopencl",
    bool verbose = false);

struct OpenCLConfig {
  bool have_expm1;
  bool have_log1p;
  std::string buildOptions;
};

// Probe OpenCL device capabilities and construct build options
OpenCLConfig configureOpenCL(cl_context context,
                             cl_device_id device);

#endif // USE_OPENCL



} // namespace openclPort

#endif // OPENCLPORT_H


