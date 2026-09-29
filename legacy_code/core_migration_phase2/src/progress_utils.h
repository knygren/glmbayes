/**
 * @file progress_utils.h
 * @brief Timing, timestamping, and console‑progress utilities used throughout glmbayes.
 *
 * @namespace glmbayes::progress
 * @brief Lightweight helpers for wall‑clock timing, human‑readable duration
 *        formatting, timestamp generation, and progress‑bar output during
 *        long‑running computations.
 *
 * @section ImplementedIn
 *   These declarations are implemented in:
 *     - progress_utils.cpp        (progress_bar implementation)
 *     - inline definitions within this header (timestamps, Timer, formatting)
 *
 * @section UsedBy
 *   These functions are consumed by:
 *     - envelope construction routines (EnvelopeBuild, EnvelopeDispersionBuild)
 *     - simulation and sampling routines (Normal, Normal–Gamma, GLM samplers)
 *     - R‑facing wrappers that report timing or progress to the console
 *     - diagnostic and benchmarking utilities across the glmbayes backend
 *
 * @section Responsibilities
 *   Provides:
 *     - timestamp helpers (`now_hms`, `timestamp_cpp`) for consistent logging,
 *     - a lightweight `Timer` struct for measuring elapsed wall‑clock time,
 *     - duration formatting utilities (`format_hms` overloads),
 *     - comma‑formatted integer printing for large iteration counts,
 *     - a console progress bar (`progress_bar`) for iterative algorithms.
 *
 *   These utilities:
 *     - avoid external dependencies and remain cross‑platform,
 *     - use Rcpp::Rcout for seamless integration with R console output,
 *     - are designed to be lightweight enough for use inside tight loops.
 */

#ifndef PROGRESS_UTILS_H
#define PROGRESS_UTILS_H

#include <string>
#include <tuple>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <RcppArmadillo.h>

namespace glmbayes {
namespace progress {

// Existing now_hms()
inline std::string now_hms() {
  std::time_t t = std::time(nullptr);
  char buf[16];
  std::strftime(buf, sizeof(buf), "%H:%M:%S", std::localtime(&t));
  return std::string(buf);
}

// New: full timestamp helper (replaces format(Sys.time()))
inline std::string timestamp_cpp() {
  using clock = std::chrono::system_clock;
  auto now = clock::now();
  std::time_t t = clock::to_time_t(now);
  
  std::tm tm;
#ifdef _WIN32
  localtime_s(&tm, &t);
#else
  localtime_r(&t, &tm);
#endif
  
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
  return oss.str();
}

struct Timer {
  std::chrono::steady_clock::time_point start;
  void begin() { start = std::chrono::steady_clock::now(); }
  std::tuple<int,int,int> hms() const {
    auto dur = std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::steady_clock::now() - start
    ).count();
    int h = static_cast<int>(dur / 3600);
    int m = static_cast<int>((dur - h*3600) / 60);
    int s = static_cast<int>(dur - h*3600 - m*60);
    return {h,m,s};
  }
};

// Format from total seconds
inline std::string format_hms(long total_seconds) {
  long h = total_seconds / 3600;
  long m = (total_seconds % 3600) / 60;
  long s = total_seconds % 60;
  
  std::ostringstream oss;
  oss << h << "h " << m << "m " << s << "s";
  return oss.str();
}

// Format from (h, m, s) tuple
inline std::string format_hms(int h, int m, int s) {
  std::ostringstream oss;
  oss << h << "h " << m << "m " << s << "s";
  return oss.str();
}

inline void print_completed(const char* prefix, const Timer& tm) {
  auto [h,m,s] = tm.hms();
  Rcpp::Rcout << prefix << " completed in: " << h << "h " << m << "m " << s << "s.\n";
}

struct comma_numpunct : std::numpunct<char> {
protected:
  char do_thousands_sep() const override { return ','; }
  std::string do_grouping() const override { return "\3"; }
};

inline std::string format_int_with_commas(long long value) {
  std::ostringstream oss;
  oss.imbue(std::locale(std::locale::classic(), new comma_numpunct));
  oss << value;
  return oss.str();
}

void progress_bar(double x, double N, const std::string& prefix = "");

// Finalize an in-progress bar; when \a newline is false, leave the line for \r overwrite.
void progress_bar_finish(bool newline = true);

} // namespace progress
} // namespace glmbayes

#endif