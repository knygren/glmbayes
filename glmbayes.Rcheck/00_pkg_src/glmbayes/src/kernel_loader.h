#ifdef USE_OPENCL


#pragma once
#include <string>

#ifdef USE_OPENCL
std::string load_kernel_source(const std::string& relative_path, const std::string& package = "glmbayes");
#endif

#ifdef USE_OPENCL
std::string load_kernel_library(const std::string& subdir, const std::string& package = "glmbayes");
#endif


#endif
