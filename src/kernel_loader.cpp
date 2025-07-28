#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <filesystem>  // C++17
#include <vector>
#include <map>
#include <set>
#include <string>
#include <stdexcept>
#include <algorithm>


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



/////////////////////////////

// [[Rcpp::export]]
std::string load_kernel_library(const std::string& subdir, const std::string& package = "glmbayes") {
  std::string dir_path = Rcpp::as<std::string>(
    Rcpp::Function("system.file")("cl", subdir, Rcpp::Named("package") = package)
  );
  
  std::map<std::string, std::set<std::string>> provides_map;
  std::map<std::string, std::set<std::string>> depends_map;
  std::map<std::string, std::filesystem::path> file_map;
  
  std::cerr << "\n📂 Files found in '" << subdir << "':\n";
  for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
    if (entry.path().extension() == ".cl") {
      std::string file_id = entry.path().stem().string();
      std::cerr << " - " << file_id << "\n";
      
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
          while (ss >> item) {
            // Remove only ‘,’ characters
            item.erase(std::remove(item.begin(), item.end(), ','),item.end()  );
            //  item.erase(std::remove_if(item.begin(), item.end(), ::ispunct), item.end());
            depends.insert(item);
            
            
          }
        }
      }
      
      file_map[file_id] = entry.path();
      provides_map[file_id] = provides;
      depends_map[file_id] = depends;
    }
  }
  
  std::vector<std::string> sorted;
  std::set<std::string> sorted_set;
  std::set<std::string> unsorted_set;
  
  std::cerr << "\n📤 Files with no dependencies:\n";
  for (const auto& [file, _] : file_map) {
    if (depends_map[file].empty()) {
      sorted.push_back(file);
      sorted_set.insert(file);
      std::cerr << " + " << file << "\n";
    } else {
      unsorted_set.insert(file);
    }
  }
  
  std::cerr << "\n🧪 Unsorted files:\n";
  for (const auto& file : unsorted_set) {
    std::cerr << " - " << file << "\n";
  }
  
  int pass_count = 0;
  while (!unsorted_set.empty()) {
    ++pass_count;
    std::cerr << "\n🔁 While Loop Pass #" << pass_count << " — Remaining unsorted: " << unsorted_set.size() << "\n";
    
    std::vector<std::string> newly_sorted;
    bool progress_made = false;
    int file_counter = 0;
    
    for (const std::string& file : unsorted_set) {
      ++file_counter;
      std::cerr << "   🔍 File #" << file_counter << ": " << file << "\n";
      
      const auto& deps = depends_map[file];
      int depends_counter = static_cast<int>(deps.size());
      std::cerr << "      📦 Dependency Count: " << depends_counter << "\n";
      
      int found_counter = 0;
      int dep_index = 0;
      for (const std::string& dep : deps) {
        ++dep_index;
        std::cerr << "         🔎 Checking classified #" << dep_index << ": " << dep << "\n";
        
        auto it = sorted_set.find(dep);
        if (it != sorted_set.end()) {
          std::cerr << "            ➤ Found in sorted? ✅ Yes\n";
          ++found_counter;
        } else {
          std::cerr << "            ➤ Found in sorted? ❌ No\n";
        }
      }
      
      std::cerr << "      🔍 Found count: " << found_counter << "\n";
      if (found_counter == depends_counter) {
        sorted.push_back(file);
        sorted_set.insert(file);
        newly_sorted.push_back(file);
        progress_made = true;
        std::cerr << " ✅ Promoted to Sorted: " << file << "\n";
      }
    }
    
    for (const std::string& file : newly_sorted) {
      unsorted_set.erase(file);
    }
    
    if (!progress_made) {
      std::cerr << "\n❌ No files promoted on pass #" << pass_count << "; possible circular or missing dependencies:\n";
      for (const std::string& file : unsorted_set) {
        std::cerr << " - " << file << "\n";
      }
      throw std::runtime_error("Dependency sort failed: unresolved dependencies remain.");
    }
  }
  
  std::cerr << "\n🔗 Final Sorted Load Order:\n";
  for (const auto& file : sorted) {
    std::cerr << " - " << file << "\n";
  }
  
  std::string combined_source;
  for (const auto& file : sorted) {
    std::string rel_path = subdir + "/" + file + ".cl";
    combined_source += load_kernel_source(rel_path, package) + "\n";
  }
  
  return combined_source;
}