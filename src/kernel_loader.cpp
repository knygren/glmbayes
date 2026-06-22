#include <RcppArmadillo.h>
#include "openclPort.h"
#include "opencl.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
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

namespace {

std::vector<std::string> parse_cl_tag(
    const std::vector<std::string>& lines,
    const std::string& tag)
{
  std::vector<std::string> result;
  std::string pattern = "@" + tag;
  for (const auto& line : lines) {
    auto pos = line.find(pattern);
    if (pos == std::string::npos) continue;
    auto colon = line.find(':', pos + pattern.size());
    if (colon == std::string::npos) continue;
    std::istringstream ss(line.substr(colon + 1));
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      tok.erase(0, tok.find_first_not_of(" \t\r\n"));
      auto last = tok.find_last_not_of(" \t\r\n");
      if (last != std::string::npos) tok.erase(last + 1);
      if (!tok.empty()) result.push_back(tok);
    }
  }
  return result;
}

struct KernelDepIndex {
  std::vector<std::string> stems_ordered;
  std::unordered_map<std::string, std::vector<std::string>> all_depends;
};

KernelDepIndex read_tsv_index(const std::string& tsv_path)
{
  KernelDepIndex idx;
  std::ifstream f(tsv_path);
  if (!f.is_open()) {
    throw std::runtime_error(
        "kernel_dependency_index.tsv not found: " + tsv_path);
  }
  std::string line;
  bool header = true;
  while (std::getline(f, line)) {
    if (header) { header = false; continue; }
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;

    auto tab = line.find('\t');
    std::string stem = (tab == std::string::npos) ? line : line.substr(0, tab);
    if (stem.empty()) continue;
    idx.stems_ordered.push_back(stem);

    std::vector<std::string> deps;
    if (tab != std::string::npos && tab + 1 < line.size()) {
      std::istringstream ss(line.substr(tab + 1));
      std::string tok;
      while (std::getline(ss, tok, ',')) {
        tok.erase(0, tok.find_first_not_of(" \t\r\n"));
        auto last = tok.find_last_not_of(" \t\r\n");
        if (last != std::string::npos) tok.erase(last + 1);
        if (!tok.empty()) deps.push_back(tok);
      }
    }
    idx.all_depends[stem] = std::move(deps);
  }
  return idx;
}

} // namespace

std::string load_library_for_kernel_cross_package(
    const std::string& kernel_relative_path,
    const std::string& kernel_package,
    const std::string& library_subdir,
    const std::string& library_package,
    const std::string& depends_tag)
{
  std::string kernel_path = Rcpp::as<std::string>(
      Rcpp::Function("system.file")(
          "cl", kernel_relative_path,
          Rcpp::Named("package") = kernel_package));
  if (kernel_path.empty()) {
    throw std::runtime_error(
        "Kernel file not found via system.file: " + kernel_relative_path +
        " (package=" + kernel_package + ")");
  }

  std::string lib_dir = Rcpp::as<std::string>(
      Rcpp::Function("system.file")(
          "cl", library_subdir,
          Rcpp::Named("package") = library_package));
  if (lib_dir.empty()) {
    throw std::runtime_error(
        "Library directory not found via system.file: " + library_subdir +
        " (package=" + library_package + ")");
  }

  std::ifstream kf(kernel_path);
  if (!kf.is_open()) {
    throw std::runtime_error("Cannot open kernel file: " + kernel_path);
  }
  std::vector<std::string> klines;
  {
    std::string kl;
    while (std::getline(kf, kl)) klines.push_back(kl);
  }
  kf.close();

  std::vector<std::string> needed_stems = parse_cl_tag(klines, depends_tag);
  if (needed_stems.empty()) {
    return "";
  }

  std::string tsv_path = lib_dir + "/kernel_dependency_index.tsv";
  KernelDepIndex idx = read_tsv_index(tsv_path);

  std::unordered_set<std::string> needed_set(needed_stems.begin(), needed_stems.end());
  std::vector<std::string> to_load;
  to_load.reserve(needed_set.size());
  for (const auto& stem : idx.stems_ordered) {
    if (needed_set.count(stem)) to_load.push_back(stem);
  }

  std::string combined;
  for (const auto& stem : to_load) {
    std::string cl_path = lib_dir + "/" + stem + ".cl";
    std::ifstream cf(cl_path);
    if (!cf.is_open()) {
      throw std::runtime_error(
          "Library file not found for stem '" + stem + "': " + cl_path);
    }
    std::ostringstream oss;
    oss << cf.rdbuf();
    combined += oss.str() + "\n\n";
  }

  return combined;
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

#endif // USE_OPENCL

} // namespace openclPort

#ifdef USE_OPENCL
namespace glmbayes {
namespace opencl {

std::string load_likelihood_subgradient_program(
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

std::string load_likelihood_subgradient_program_v2(
    const std::string& family,
    const std::string& link,
    const std::string& app_package,
    const std::string& nmath_package)
{
  const std::string kernel_file = openclPort::resolve_kernel_path(family, link);

  std::string opencl_source =
      openclPort::load_kernel_source("OPENCL.cl", nmath_package);
  std::string libr_shims_source =
      openclPort::load_kernel_library("libR_shims", nmath_package, false);
  std::string r_ext_types_source =
      openclPort::load_kernel_library("R_ext_types", nmath_package, false);
  std::string r_shims_source =
      openclPort::load_kernel_library("R_shims", nmath_package, false);
  std::string r_ext_runtime_source =
      openclPort::load_kernel_library("R_ext_runtime", nmath_package, false);
  std::string r_ext_internals_source =
      openclPort::load_kernel_library("R_ext_internals", nmath_package, false);
  std::string system_source =
      openclPort::load_kernel_library("System", nmath_package, false);
  std::string nmath_source = openclPort::load_library_for_kernel_cross_package(
      kernel_file,
      app_package,
      "nmath",
      nmath_package,
      "all_depends_nmath");
  std::string ksrc =
      openclPort::load_kernel_source(kernel_file, app_package);

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

std::string load_kernel_source_wrapper(std::string relative_path,
                                       std::string package) {
#ifdef USE_OPENCL
  return load_kernel_source(relative_path, package);
#else
  Rcpp::stop("OpenCL support is not available in this build of glmbayes.");
#endif
}

std::string load_kernel_library_wrapper(std::string subdir,
                                        std::string package,
                                        bool verbose) {
#ifdef USE_OPENCL
  return load_kernel_library(subdir, package, verbose);
#else
  Rcpp::stop("OpenCL support is not available in this build of glmbayes.");
#endif
}

} // namespace openclPort
