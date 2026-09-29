// rNormalGLMBlocks.cpp
// C++ counterpart to R rNormalGLM_reg_block(): serial loop over blocks, each
// iteration calls the existing rNormalGLM() on sliced (y_b, x_b, ...).
// Not exported to R yet (Phase 2).

#include "simfuncs.h"

#include <string>

namespace glmbayes {
namespace sim {

namespace {

using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::Function;
using Rcpp::CharacterVector;

inline int check_index_1based(int idx, int n, const char* what) {
  if (idx < 1 || idx > n) {
    Rcpp::stop("%s index %d out of range [1, %d]", what, idx, n);
  }
  return idx - 1;
}

NumericVector slice_numeric(const NumericVector& v, const IntegerVector& rows) {
  const int m = rows.size();
  NumericVector out(m);
  for (int i = 0; i < m; ++i) {
    out[i] = v[check_index_1based(rows[i], v.size(), "row")];
  }
  return out;
}

NumericMatrix slice_matrix_rows(const NumericMatrix& x, const IntegerVector& rows) {
  const int m = rows.size();
  const int l1 = x.ncol();
  NumericMatrix out(m, l1);
  for (int i = 0; i < m; ++i) {
    const int r = check_index_1based(rows[i], x.nrow(), "row");
    for (int j = 0; j < l1; ++j) {
      out(i, j) = x(r, j);
    }
  }
  return out;
}

NumericVector mu_for_block(const NumericMatrix& mu, int b, bool prior_by_block) {
  if (mu.ncol() < 1) {
    Rcpp::stop("mu must have at least one column");
  }
  if (!prior_by_block || mu.ncol() == 1) {
    return mu(Rcpp::_, 0);
  }
  if (b < 0 || b >= mu.ncol()) {
    Rcpp::stop("block index %d out of range for mu with %d columns", b + 1, mu.ncol());
  }
  return mu(Rcpp::_, b);
}

NumericMatrix P_for_block(const List& P_blocks, int b, bool prior_by_block) {
  if (P_blocks.size() < 1) {
    Rcpp::stop("P_blocks must be a non-empty list");
  }
  SEXP p_sexp;
  if (!prior_by_block || P_blocks.size() == 1) {
    p_sexp = P_blocks[0];
  } else {
    if (b < 0 || b >= P_blocks.size()) {
      Rcpp::stop("block index %d out of range for P_blocks of length %d", b + 1, P_blocks.size());
    }
    p_sexp = P_blocks[b];
  }
  return Rcpp::as<NumericMatrix>(p_sexp);
}

double dispersion_for_block(const NumericVector& dispersion, int b, bool prior_by_block) {
  if (dispersion.size() == 1) {
    return dispersion[0];
  }
  if (!prior_by_block) {
    return dispersion[0];
  }
  if (b < 0 || b >= dispersion.size()) {
    Rcpp::stop("block index %d out of range for dispersion of length %d", b + 1, dispersion.size());
  }
  return dispersion[b];
}

void validate_row_blocks(const List& row_blocks, int l2) {
  const int k = row_blocks.size();
  if (k < 1) {
    Rcpp::stop("row_blocks must contain at least one block");
  }
  int total = 0;
  for (int b = 0; b < k; ++b) {
    IntegerVector rows = row_blocks[b];
    if (rows.size() < 1) {
      Rcpp::stop("row_blocks[[%d]] is empty", b + 1);
    }
    total += rows.size();
  }
  if (total != l2) {
    Rcpp::stop(
      "row_blocks cover %d rows but length(y) = %d",
      total, l2
    );
  }
}

NumericMatrix coefficients_row_from_block(const List& out_b, int l1) {
  NumericMatrix coef = out_b["coefficients"];
  if (coef.nrow() < 1) {
    Rcpp::stop("block result has no coefficient draws");
  }
  NumericVector row = coef(0, Rcpp::_);
  if (row.size() != l1) {
    Rcpp::stop("expected %d coefficients in block draw, got %d", l1, row.size());
  }
  NumericMatrix out(1, l1);
  out(0, Rcpp::_) = row;
  return out;
}

NumericMatrix coef_mode_row_from_block(const List& out_b, int l1) {
  SEXP cm_sexp = out_b["coef.mode"];
  if (Rcpp::is<NumericMatrix>(cm_sexp)) {
    NumericMatrix cm = Rcpp::as<NumericMatrix>(cm_sexp);
    if (cm.nrow() == 1 && cm.ncol() == l1) {
      return cm;
    }
    if (cm.ncol() == 1 && cm.nrow() == l1) {
      return Rcpp::transpose(cm);
    }
    Rcpp::stop("coef.mode matrix has unexpected shape");
  }
  NumericVector cm = Rcpp::as<NumericVector>(cm_sexp);
  if (cm.size() != l1) {
    Rcpp::stop("expected %d coef.mode values, got %d", l1, cm.size());
  }
  NumericMatrix out(1, l1);
  out(0, Rcpp::_) = cm;
  return out;
}

} // namespace

// Definition must match simfuncs.h (see rNormalGLM for same pass-by-value style).
Rcpp::List rNormalGLMBlocks(
    int n,
    NumericVector y,
    NumericMatrix x,
    NumericVector offset,
    NumericVector wt,
    NumericVector dispersion,
    NumericMatrix mu,
    List P_blocks,
    bool prior_by_block,
    List row_blocks,
    Function f2,
    Function f3,
    std::string family,
    std::string link,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose
) {
  if (n < 1) {
    Rcpp::stop("n must be at least 1");
  }

  const int l2 = y.size();
  const int l1 = x.ncol();
  if (x.nrow() != l2) {
    Rcpp::stop("nrow(x) must equal length(y)");
  }
  if (mu.nrow() != l1) {
    Rcpp::stop("nrow(mu) must equal ncol(x) (= l1)");
  }
  if (offset.size() != l2 && offset.size() != 1) {
    Rcpp::stop("length(offset) must be 1 or length(y)");
  }
  if (wt.size() != l2 && wt.size() != 1) {
    Rcpp::stop("length(wt) must be 1 or length(y)");
  }
  if (dispersion.size() < 1) {
    Rcpp::stop("dispersion must have length at least 1");
  }

  validate_row_blocks(row_blocks, l2);

  const int k = row_blocks.size();
  if (prior_by_block && mu.ncol() > 1 && mu.ncol() != k) {
    Rcpp::stop("ncol(mu) must be 1 or number of blocks (%d)", k);
  }
  if (prior_by_block && dispersion.size() > 1 && dispersion.size() != k) {
    Rcpp::stop("length(dispersion) must be 1 or number of blocks (%d)", k);
  }

  NumericVector offset_full = offset;
  if (offset_full.size() == 1) {
    offset_full = Rcpp::rep(offset_full[0], l2);
  }
  NumericVector wt_full = wt;
  if (wt_full.size() == 1) {
    wt_full = Rcpp::rep(wt_full[0], l2);
  }

  NumericMatrix coef_draw(k, l1);
  NumericMatrix coef_mode_mat(k, l1);
  NumericVector dispersion_out(k);
  List block_results(k);
  List envelope_blocks(k);

  CharacterVector family_cv(1);
  family_cv[0] = family;
  CharacterVector link_cv(1);
  link_cv[0] = link;

  const int n_envopt_b = (n_envopt < 0) ? n : n_envopt;
  const bool use_parallel_b = use_parallel && (n > 1);

  for (int b = 0; b < k; ++b) {
    IntegerVector rows = row_blocks[b];
    NumericVector y_b = slice_numeric(y, rows);
    NumericMatrix x_b = slice_matrix_rows(x, rows);
    NumericVector offset_b = slice_numeric(offset_full, rows);
    NumericVector wt_b = slice_numeric(wt_full, rows);

    NumericVector mu_b = mu_for_block(mu, b, prior_by_block);
    NumericMatrix P_b = P_for_block(P_blocks, b, prior_by_block);
    double disp_b = dispersion_for_block(dispersion, b, prior_by_block);

    NumericVector start_b = Rcpp::clone(mu_b);

    List out_b = rNormalGLM(
      n,
      y_b,
      x_b,
      mu_b,
      P_b,
      offset_b,
      wt_b,
      disp_b,
      f2,
      f3,
      start_b,
      family,
      link,
      Gridtype,
      n_envopt_b,
      use_parallel_b,
      use_opencl,
      verbose && (b == 0)
    );

    block_results[b] = out_b;
    envelope_blocks[b] = out_b["Envelope"];

    NumericMatrix draw_row = coefficients_row_from_block(out_b, l1);
    NumericMatrix mode_row = coef_mode_row_from_block(out_b, l1);
    coef_draw(b, Rcpp::_) = draw_row(0, Rcpp::_);
    coef_mode_mat(b, Rcpp::_) = mode_row(0, Rcpp::_);
    dispersion_out[b] = Rcpp::as<double>(out_b["dispersion"]);
  }

  return Rcpp::List::create(
    Rcpp::Named("coefficients") = coef_draw,
    Rcpp::Named("coef.mode") = coef_mode_mat,
    Rcpp::Named("dispersion") = dispersion_out,
    Rcpp::Named("n") = n,
    Rcpp::Named("k") = k,
    Rcpp::Named("l1") = l1,
    Rcpp::Named("l2") = l2,
    Rcpp::Named("block_results") = block_results,
    Rcpp::Named("Envelope") = envelope_blocks,
    Rcpp::Named("y") = y,
    Rcpp::Named("x") = x,
    Rcpp::Named("offset") = offset,
    Rcpp::Named("prior.weights") = wt,
    Rcpp::Named("family") = family_cv,
    Rcpp::Named("link") = link_cv
  );
}

} // namespace sim
} // namespace glmbayes
