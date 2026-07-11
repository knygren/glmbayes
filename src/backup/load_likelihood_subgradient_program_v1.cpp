// ARCHIVED — not compiled. See src/backup/README.md.
//
// Former glmbayes::opencl::load_likelihood_subgradient_program (v1):
// all OpenCL pieces loaded from a single package (typically "glmbayes").

#ifdef USE_OPENCL

#include "openclPort.h"
#include "opencl.h"

namespace glmbayes {
namespace opencl {
namespace backup {

std::string load_likelihood_subgradient_program_v1(
    const std::string& family,
    const std::string& link,
    const std::string& package)
{
  const std::string kernel_file = openclPort::resolve_kernel_path(family, link);

  std::string opencl_source = openclPort::load_kernel_source("OPENCL.cl", package);
  std::string libr_shims_source =
      openclPort::load_kernel_library("libR_shims", package, false);
  std::string r_ext_types_source =
      openclPort::load_kernel_library("R_ext_types", package, false);
  std::string r_shims_source =
      openclPort::load_kernel_library("R_shims", package, false);
  std::string r_ext_runtime_source =
      openclPort::load_kernel_library("R_ext_runtime", package, false);
  std::string r_ext_internals_source =
      openclPort::load_kernel_library("R_ext_internals", package, false);
  std::string system_source =
      openclPort::load_kernel_library("System", package, false);
  std::string nmath_source = openclPort::load_library_for_kernel(
      kernel_file, "nmath", package, "all_depends_nmath");
  std::string ksrc = openclPort::load_kernel_source(kernel_file, package);

  return opencl_source +
    "\n" + libr_shims_source +
    "\n" + r_ext_types_source +
    "\n" + r_shims_source +
    "\n" + r_ext_runtime_source +
    "\n" + r_ext_internals_source +
    "\n" + system_source +
    "\n" + nmath_source +
    "\n" + ksrc;
}

} // namespace backup
} // namespace opencl
} // namespace glmbayes

#endif // USE_OPENCL
