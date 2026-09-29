#include <RcppArmadillo.h>
#include "openclPort.h"
using namespace Rcpp;

#ifdef __linux__
#include <stdio.h>
#include <stdlib.h>
#endif



namespace openclPort {
Rcpp::CharacterVector gpu_names() {
#ifdef __linux__
  FILE *fp;
  char buffer[256];
  Rcpp::CharacterVector gpus;
  
  // Run nvidia-smi with query options
  fp = popen("nvidia-smi --query-gpu=name --format=csv,noheader", "r");
  if (fp == NULL) {
    gpus.push_back("Failed to run nvidia-smi");
    return gpus;
  }
  
  while (fgets(buffer, sizeof(buffer), fp) != NULL) {
    std::string line(buffer);
    if (!line.empty() && line.back() == '\n') line.pop_back();
    gpus.push_back(line);
  }
  
  pclose(fp);
  return gpus;
#else
  return Rcpp::CharacterVector::create("CUDA probe not supported on this platform.");
#endif
}

}