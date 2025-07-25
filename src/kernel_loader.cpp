#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>  // C++17

namespace fs = std::filesystem;

// Load a single file like "nmath/bd0.cl"





std::string load_kernel_source(const std::string& relative_path,
                               const std::string& package = "glmbayes") {
  // Retrieve full path via system.file()
  std::string path = Rcpp::as<std::string>(
    Rcpp::Function("system.file")("cl", relative_path,
                   Rcpp::Named("package") = package)
  );
  
  // Check for empty path returned by system.file (means file not found)
  if (path.empty()) {
    throw std::runtime_error("Kernel source not found via system.file: " + relative_path);
  }
  
  // Attempt to open the file
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open kernel source: " + path);
  }
  
  // Read file contents
  std::ostringstream oss;
  oss << file.rdbuf();
  return oss.str();
}




///////////////////////////////////////




std::string load_kernel_library(const std::string& subdir, const std::string& package = "glmbayes") {
  std::string dir_path = Rcpp::as<std::string>(
    Rcpp::Function("system.file")("cl", subdir, Rcpp::Named("package") = package)
  );
  
  std::vector<fs::path> cl_files;
  for (const auto& entry : fs::directory_iterator(dir_path)) {
    if (entry.path().extension() == ".cl") {
      cl_files.push_back(entry.path());
    }
  }
  
  // Sort alphabetically by filename
  std::sort(cl_files.begin(), cl_files.end(), [](const fs::path& a, const fs::path& b) {
    return a.filename().string() < b.filename().string();
  });
  
  std::string combined_source;
  for (const auto& path : cl_files) {
    std::string rel_path = subdir + "/" + path.filename().string();
    combined_source += load_kernel_source(rel_path, package) + "\n";
  }
  
  return combined_source;
}



