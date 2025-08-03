// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "rng_utils.h"
#include <random>

// Thread-local RNG and distribution
thread_local std::mt19937 safe_rng_engine(std::random_device{}());
thread_local std::uniform_real_distribution<> safe_rng_dist(0.0, 1.0);

// Core sampling function
double safe_runif() {
  return safe_rng_dist(safe_rng_engine);
}
