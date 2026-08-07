#include "progress_utils.h"

#include <cmath>
#include <iomanip>

using namespace Rcpp;
using namespace glmbayes::progress;

namespace glmbayes {

namespace progress {

void progress_bar(double x, double N, const std::string& prefix)
{
  if (N <= 0.0 || !std::isfinite(N) || !std::isfinite(x)) {
    return;
  }
  int totaldotz = 40;
  double fraction = x / N;
  if (!std::isfinite(fraction)) {
    return;
  }
  int dotz = static_cast<int>(std::round(fraction * totaldotz));
  if (dotz < 0) dotz = 0;
  if (dotz > totaldotz) dotz = totaldotz;

  Rcpp::Rcout.precision(3);
  Rcout << "\r" << std::string(100, ' ') << "\r" << std::flush;
  if (!prefix.empty()) {
    Rcout << prefix << std::flush;
  }
  Rcout << std::fixed << std::setprecision(0) << fraction * 100.0 << std::flush;
  Rcout << "% [" << std::flush;
  for (int ii = 0; ii < dotz; ++ii) {
    Rcout << "=" << std::flush;
  }
  for (int ii = dotz; ii < totaldotz; ++ii) {
    Rcout << " " << std::flush;
  }
  Rcout << "]" << std::flush;
}

void progress_bar_finish(bool newline)
{
  if (newline) {
    Rcpp::Rcout << std::endl;
  }
}

}
}
