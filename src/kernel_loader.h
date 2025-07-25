#pragma once
#include <string>

std::string load_kernel_source(const std::string& relative_path, const std::string& package = "glmbayes");
std::string load_kernel_library(const std::string& subdir, const std::string& package = "glmbayes");


