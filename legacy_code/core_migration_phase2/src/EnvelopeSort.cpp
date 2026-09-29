#include <Rcpp.h>
#include <algorithm>

namespace glmbayes {

namespace env {

using namespace Rcpp;

Rcpp::List EnvelopeSort_cpp(
    int l1,
    int l2,
    const NumericMatrix& GIndex,      // GridIndex: l2 x l1
    const NumericMatrix& G3,          // thetabars: l2 x l1
    const NumericMatrix& cbars,       // cbars: l2 x l1
    const NumericMatrix& logU,        // l2 x l1 (or l2 x ncol)
    const NumericMatrix& logrt,       // l2 x l1
    const NumericMatrix& loglt,       // l2 x l1
    const NumericMatrix& logP,        // l2 x 2 (or l2 x ncol)
    const NumericMatrix& LLconst,     // l2 x 1 (or l2 x ncol)
    const NumericVector& PLSD,        // length l2
    const NumericVector& a1,          // length l1
    double E_draws,                   // scalar
    const Nullable<NumericVector>& lg_prob_factor,
    const Nullable<NumericVector>& UB2min
) {
  // --- Dimension checks ---
  if (GIndex.nrow() != l2 || G3.nrow() != l2 || cbars.nrow() != l2 ||
      logrt.nrow()  != l2 || loglt.nrow() != l2 ||
      logU.nrow()   != l2 || LLconst.nrow() != l2 ||
      PLSD.size()   != l2 || logP.nrow()    != l2) {
    stop("EnvelopeSort_cpp: dimension mismatch with l2.");
  }
  
  if (GIndex.ncol() != l1 || G3.ncol() != l1 || cbars.ncol() != l1 ||
      logrt.ncol()  != l1 || loglt.ncol() != l1) {
    stop("EnvelopeSort_cpp: dimension mismatch with l1.");
  }
  
  if (a1.size() != l1) {
    stop("EnvelopeSort_cpp: a1 length must equal l1.");
  }
  
  // Optional vectors
  NumericVector lg_prob_factor_vec;
  bool has_lg_prob_factor = false;
  if (lg_prob_factor.isNotNull()) {
    lg_prob_factor_vec = as<NumericVector>(lg_prob_factor);
    if (lg_prob_factor_vec.size() != l2)
      stop("EnvelopeSort_cpp: lg_prob_factor length != l2.");
    has_lg_prob_factor = true;
  }
  
  NumericVector UB2min_vec;
  bool has_UB2min = false;
  if (UB2min.isNotNull()) {
    UB2min_vec = as<NumericVector>(UB2min);
    if (UB2min_vec.size() != l2)
      stop("EnvelopeSort_cpp: UB2min length != l2.");
    has_UB2min = true;
  }
  
  // --- Build order indices (decreasing PLSD) ---
  IntegerVector ord(l2);
  for (int i = 0; i < l2; ++i) ord[i] = i;
  
  std::sort(ord.begin(), ord.end(),
            [&](int i, int j) {
              return PLSD[i] > PLSD[j];   // decreasing
            });
  
  // --- Allocate output objects ---
  NumericMatrix GIndex_out(l2, l1);
  NumericMatrix G3_out(l2, l1);
  NumericMatrix cbars_out(l2, l1);
  NumericMatrix logrt_out(l2, l1);
  NumericMatrix loglt_out(l2, l1);
  
  NumericMatrix logU_out(logU.nrow(), logU.ncol());
  NumericMatrix LLconst_out(LLconst.nrow(), LLconst.ncol());
  NumericVector PLSD_out(l2);
  NumericMatrix logP_out(logP.nrow(), logP.ncol());
  
  NumericVector lg_prob_factor_out;
  if (has_lg_prob_factor)
    lg_prob_factor_out = NumericVector(l2);
  
  NumericVector UB2min_out;
  if (has_UB2min)
    UB2min_out = NumericVector(l2);
  
  // --- Reorder by ord ---
  for (int r = 0; r < l2; ++r) {
    int src = ord[r];
    
    for (int c = 0; c < l1; ++c) {
      GIndex_out(r, c) = GIndex(src, c);
      G3_out(r, c)     = G3(src, c);
      cbars_out(r, c)  = cbars(src, c);
      logrt_out(r, c)  = logrt(src, c);
      loglt_out(r, c)  = loglt(src, c);
    }
    for (int c = 0; c < logU.ncol(); ++c)
      logU_out(r, c) = logU(src, c);
    for (int c = 0; c < logP.ncol(); ++c)
      logP_out(r, c) = logP(src, c);
    for (int c = 0; c < LLconst.ncol(); ++c)
      LLconst_out(r, c) = LLconst(src, c);
    PLSD_out[r]     = PLSD[src];
    
    if (has_lg_prob_factor)
      lg_prob_factor_out[r] = lg_prob_factor_vec[src];
    
    if (has_UB2min)
      UB2min_out[r] = UB2min_vec[src];
  }
  
  // --- Build output list (logP/logU: downstream expects vectors for single-column case) ---
  // Downstream expects logP as vector (first column)
  NumericVector logP_col0(l2);
  for (int i = 0; i < l2; ++i) logP_col0[i] = logP_out(i, 0);

  List out = List::create(
    _["GridIndex"] = GIndex_out,
    _["thetabars"] = G3_out,
    _["cbars"]     = cbars_out,
    _["logU"]      = logU_out,
    _["logrt"]     = logrt_out,
    _["loglt"]     = loglt_out,
    _["LLconst"]   = LLconst_out,
    _["logP"]      = logP_col0,
    _["PLSD"]      = PLSD_out,
    _["a1"]        = a1,
    _["E_draws"]   = E_draws
  );
  
  if (has_lg_prob_factor)
    out["lg_prob_factor"] = lg_prob_factor_out;
  
  if (has_UB2min)
    out["UB2min"] = UB2min_out;
  
  return out;
}

} // namespace env

} // namespace glmbayes