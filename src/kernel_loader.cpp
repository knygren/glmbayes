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

std::string load_kernel_library(const std::string& subdir, const std::string& package = "glmbayes") {
  std::string dir_path = Rcpp::as<std::string>(
    Rcpp::Function("system.file")("cl", subdir, Rcpp::Named("package") = package)
  );
  
  std::map<std::string, std::set<std::string>> provides_map;
  std::map<std::string, std::set<std::string>> depends_map;
  std::map<std::string, fs::path> file_map;
  
  // Parse metadata
  for (const auto& entry : fs::directory_iterator(dir_path)) {
    if (entry.path().extension() == ".cl") {
      std::ifstream infile(entry.path());
      std::string line;
      std::set<std::string> provides, depends;
      
      while (std::getline(infile, line)) {
        if (line.find("@provides") != std::string::npos) {
          std::stringstream ss(line.substr(line.find("@provides") + 9));
          std::string item;
          while (ss >> item) provides.insert(item);
        } else if (line.find("@depends") != std::string::npos) {
          std::stringstream ss(line.substr(line.find("@depends") + 9));
          std::string item;
          while (ss >> item) depends.insert(item);
        }
      }
      
      std::string file_id = entry.path().filename().string();
      file_map[file_id] = entry.path();
      provides_map[file_id] = provides;
      depends_map[file_id] = depends;
    }
  }
  
  // Track fulfilled provides during sort
  std::set<std::string> fulfilled;
  std::set<std::string> unsorted;
  for (const auto& [file, _] : file_map) {
    unsorted.insert(file);
  }
  
  std::vector<std::string> sorted;
  bool progress = true;
  
  while (progress) {
    progress = false;
    
    for (auto it = unsorted.begin(); it != unsorted.end(); ) {
      const std::string& file = *it;
      const std::set<std::string>& deps = depends_map[file];
      
      bool all_satisfied = std::all_of(deps.begin(), deps.end(), [&](const std::string& dep) {
        return fulfilled.count(dep);
      });
      
      if (all_satisfied) {
        sorted.push_back(file);
        fulfilled.insert(provides_map[file].begin(), provides_map[file].end());
        it = unsorted.erase(it);
        progress = true;
      } else {
        ++it;
      }
    }
  }
  
  // Anything remaining is unresolved—append at end
  for (const std::string& file : unsorted) {
    sorted.push_back(file);
  }
  
  // 📤 Debug: Final sort order
  std::cerr << "\n🔗 Final kernel load order (resolved first, unresolved last):\n";
  for (const auto& file : sorted) {
    std::cerr << " - " << file << "\n";
  }
  
  // Concatenate sources
  std::string combined_source;
  for (const auto& file : sorted) {
    std::string rel_path = subdir + "/" + file;
    combined_source += load_kernel_source(rel_path, package) + "\n";
  }
  
  return combined_source;
}