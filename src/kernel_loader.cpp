#include <RcppArmadillo.h>
#include "openclPort.h"
#include "opencl.h"

#include <string>
#include <stdexcept>

#include <opencltools/opencltools_capi.h>

namespace openclPort {

namespace {

std::string opencltools_take_cstr(const char* p) {
  if (p == nullptr) {
    return std::string();
  }
  std::string out(p);
  opencltools_free_cstr(p);
  return out;
}

} // namespace

#ifdef USE_OPENCL

std::string load_kernel_source(const std::string& relative_path,
                               const std::string& package) {
  return opencltools_take_cstr(
      opencltools_load_kernel_source(relative_path.c_str(), package.c_str()));
}

std::string load_kernel_library(const std::string& subdir,
                                const std::string& package,
                                bool verbose) {
  return opencltools_take_cstr(
      opencltools_load_kernel_library(subdir.c_str(), package.c_str(),
                                      verbose ? 1 : 0));
}

std::string load_library_for_kernel(
    const std::string& kernel_relative_path,
    const std::string& library_subdir,
    const std::string& package,
    const std::string& depends_tag)
{
  return opencltools_take_cstr(opencltools_load_library_for_kernel(
      kernel_relative_path.c_str(), library_subdir.c_str(), package.c_str(),
      depends_tag.c_str()));
}

std::string load_library_for_kernel_cross_package(
    const std::string& kernel_relative_path,
    const std::string& kernel_package,
    const std::string& library_subdir,
    const std::string& library_package,
    const std::string& depends_tag)
{
  return opencltools_take_cstr(opencltools_load_library_for_kernel_cross_package(
      kernel_relative_path.c_str(), kernel_package.c_str(),
      library_subdir.c_str(), library_package.c_str(), depends_tag.c_str()));
}

std::string load_program_preload(
    const std::string& manifest_relative_path,
    const std::string& source_package,
    bool verbose)
{
  return opencltools_take_cstr(opencltools_load_program_preload(
      manifest_relative_path.c_str(), source_package.c_str(),
      verbose ? 1 : 0));
}

std::string resolve_kernel_path(
    const std::string& family,
    const std::string& link)
{
  if (family == "binomial" || family == "quasibinomial") {
    if (link == "logit") {
      return "src/f2_f3_binomial_logit.cl";
    }
    if (link == "probit") {
      return "src/f2_f3_binomial_probit.cl";
    }
    if (link == "cloglog") {
      return "src/f2_f3_binomial_cloglog.cl";
    }
    throw std::runtime_error(
        "Unsupported link function for binomial family: " + link);
  }
  if (family == "poisson" || family == "quasipoisson") {
    return "src/f2_f3_poisson.cl";
  }
  if (family == "Gamma") {
    return "src/f2_f3_gamma.cl";
  }
  if (family == "gaussian") {
    return "src/f2_f3_gaussian.cl";
  }
  throw std::runtime_error("Unsupported family: " + family);
}

#endif // USE_OPENCL

} // namespace openclPort

#ifdef USE_OPENCL
namespace glmbayes {
namespace opencl {

std::string load_likelihood_subgradient_program(
    const std::string& family,
    const std::string& link,
    const std::string& app_package,
    const std::string& nmath_package)
{
  const std::string kernel_file = openclPort::resolve_kernel_path(family, link);

  std::string preload_source = openclPort::load_program_preload(
      "program_preload_manifest.tsv", nmath_package);
  std::string nmath_source = openclPort::load_library_for_kernel_cross_package(
      kernel_file,
      app_package,
      "nmath",
      nmath_package,
      "all_depends_nmath");
  std::string ksrc =
      openclPort::load_kernel_source(kernel_file, app_package);

  return preload_source + "\n" + nmath_source + "\n" + ksrc;
}

} // namespace opencl
} // namespace glmbayes
#endif // USE_OPENCL

namespace openclPort {

int get_opencl_core_count() {
#ifdef USE_OPENCL
  return std::max(1, detect_num_gpus_internal());
#else
  return 1;
#endif
}

} // namespace openclPort
