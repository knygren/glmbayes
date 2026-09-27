// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "Envelopefuncs.h"
#include "progress_utils.h"
#include "package_ns.h"
#include "R_interface.h"

using namespace glmbayes::env;


namespace glmbayes {

namespace env {

Rcpp::List EnvelopeOrchestrator(
    NumericVector bstar2,
    NumericMatrix A,
    NumericVector y,
    NumericMatrix x2,
    NumericMatrix mu2,
    NumericMatrix P2,
    NumericVector alpha,
    NumericVector wt,
    
    int n,
    int Gridtype,
    Nullable<int> n_envopt,
    
    double shape,
    double rate,
    double RSS_Post2,
    double RSS_ML,
    
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    
    bool use_parallel,
    bool use_opencl,
    bool verbose
) {
  check_disp_bounds_or_stop(disp_lower, disp_upper, "EnvelopeOrchestrator (entry)");

  // Unknown dispersion always uses full envelope size (3^p grid); smaller grids not supported.
  Gridtype = 3;

  int disp_grid_type=2;
  
  if(use_parallel) disp_grid_type=2;
  
  
    
  int n_obs=y.size();
  
  double n_w = 0.0;
  for (int i = 0; i < wt.size(); ++i) n_w += wt[i];
  
  // Step 1: Posterior Gamma parameters (precision prior)
  double shape2 = shape + n_w / 2.0;
  double rate3  = rate  + RSS_Post2 / 2.0;
  
  double d1_star = rate3 / (shape2 - 1.0);
  
  
  // Rcpp::Rcout << "d1_star - Envelopeorchestrator= " << d1_star << "\n";
  
  Rcpp::NumericVector wt2(n_obs);
  for (int i = 0; i < n_obs; ++i)    wt2[i] = wt[i] / d1_star;
  
  
  
  // sortgrid: false when n==1, true when n>1 (matches rNormalGLM logic).
  // Saves redundant first sort; EnvelopeDispersionBuild + second EnvelopeSort
  // still run. Unknown-dispersion has poorer acceptance than fixed-dispersion
  // (~279 candidates/draw vs theoretical bound), but for n==1 (e.g., Gibbs)
  // the sort cost can outweigh the benefit.
  if (n == 1 && verbose) {
    Rcpp::Rcout << "[EnvelopeOrchestrator] n=1: skipping envelope sort (use n>1 for iid)\n";
  }

  // --- Step 1: EnvelopeBuild (direct C++ call) ---
  Rcpp::List Env2 = EnvelopeBuild(
    bstar2,                      // NumericVector
    A,                           // NumericMatrix
    y,                           // NumericVector
    x2,                          // NumericMatrix
    mu2,                         // NumericMatrix (p x 1)
    P2,                          // NumericMatrix
    alpha,                       // NumericVector
    wt2,                         // NumericVector
    "gaussian",                  // family
    "identity",                  // link
    Gridtype,                    // int
    n,                           // int
    n_envopt.isNull() ? -1 : Rcpp::as<int>(n_envopt),
    false,                       // sortgrid: skip redundant first sort; EnvelopeDispersionBuild
                                // overwrites PLSD, so we sort once in Step 3 with correct weights
    use_opencl,                  // bool
    verbose                      // bool
);


 
// --- Step 2: EnvelopeDispersionBuild (direct C++ call; always run, optionally timed) ---
Rcpp::List disp_env_out;
glmbayes::progress::Timer t_disp;
if (verbose) {
  t_disp.begin();
  Rcpp::Rcout << "[EnvelopeDipsersionBuild] Entering: "
              << glmbayes::progress::timestamp_cpp()
              << "\n";
}
// Call EnvelopeDispersionBuild once, regardless of verbosity
disp_env_out = EnvelopeDispersionBuild(
  Env2,                        // List Env
  shape,                       // double Shape (prior shape; EnvelopeDispersionBuild forms shape2 internally)
  rate,                        // double Rate
  P2,                          // NumericMatrix P
  y,                           // NumericVector y
  x2,                          // NumericMatrix x
  alpha,                       // NumericVector alpha
  y.size(),                    // int n_obs
  RSS_Post2,                   // double RSS_post
  RSS_ML,                      // double RSS_ML
  mu2,                         // NumericMatrix mu
  wt,                          // NumericVector wt
  max_disp_perc,               // double max_disp_perc
  disp_lower,                  // Nullable<double> disp_lower
  disp_upper,                  // Nullable<double> disp_upper
  verbose,                     // bool verbose
  use_parallel                 // bool use_parallel
);
if (verbose) {
  Rcpp::Rcout << "[EnvelopeDipsersionBuild] Exiting: "
              << glmbayes::progress::timestamp_cpp()
              << "\n";
  glmbayes::progress::print_completed("[EnvelopeDipsersionBuild]", t_disp);
}

// --- Step 3: Call R's EnvelopeSort ---

// Extract Env_out and UB_list as proper Lists
Rcpp::List Env3_raw    = disp_env_out["Env_out"];
Rcpp::List UB_list_new = disp_env_out["UB_list"];
//Rcpp::List UB_list_new = disp_env_out["UB_list"];
Rcpp::List gamma_list_new = disp_env_out["gamma_list"];






// Extract cbars and its dimensions
Rcpp::NumericMatrix cbars = Env3_raw["cbars"];
int l1 = cbars.ncol();
int l2 = cbars.nrow();

// logP is ALWAYS a NumericVector because EnvelopeBuild returns logP(_,0)
// Wrap it as an l2 x 1 matrix for EnvelopeSort (values unchanged)
Rcpp::NumericVector logP_vec = Env3_raw["logP"];
Rcpp::NumericMatrix logP_mat(logP_vec.size(), 1, logP_vec.begin());

// Look up EnvelopeSort in the glmbayesCore namespace
Rcpp::Function EnvelopeSort = glmbayes_R::r_envelope_sort();



// C++ sort commented out (slower than R; R result used downstream)
// {
//   Rcpp::NumericMatrix logU_disp = Rcpp::as<Rcpp::NumericMatrix>(Env3_raw["logU"]);
//   Rcpp::NumericMatrix LLconst_disp = Rcpp::as<Rcpp::NumericMatrix>(Env3_raw["LLconst"]);
//   if (verbose) {
//     Rcpp::Rcout << "[EnvelopeSort_cpp] Entering: "
//                 << glmbayes::progress::timestamp_cpp()
//                 << "\n";
//   }
//   (void) glmbayes::env::EnvelopeSort_cpp(
//     l1, l2,
//     Rcpp::as<Rcpp::NumericMatrix>(Env3_raw["GridIndex"]),
//     Rcpp::as<Rcpp::NumericMatrix>(Env3_raw["thetabars"]),
//     cbars,
//     logU_disp,
//     Rcpp::as<Rcpp::NumericMatrix>(Env3_raw["logrt"]),
//     Rcpp::as<Rcpp::NumericMatrix>(Env3_raw["loglt"]),
//     logP_mat,
//     LLconst_disp,
//     Rcpp::as<Rcpp::NumericVector>(Env3_raw["PLSD"]),
//     Rcpp::as<Rcpp::NumericVector>(Env3_raw["a1"]),
//     Rcpp::as<double>(Env3_raw["E_draws"]),
//     UB_list_new.containsElementNamed("lg_prob_factor")
//       ? Rcpp::Nullable<Rcpp::NumericVector>(UB_list_new["lg_prob_factor"])
//       : Rcpp::Nullable<Rcpp::NumericVector>(),
//     UB_list_new.containsElementNamed("UB2min")
//       ? Rcpp::Nullable<Rcpp::NumericVector>(UB_list_new["UB2min"])
//       : Rcpp::Nullable<Rcpp::NumericVector>()
//   );
//   if (verbose) {
//     Rcpp::Rcout << "[EnvelopeSort_cpp] Exiting: "
//                 << glmbayes::progress::timestamp_cpp()
//                 << "\n";
//   }
// }

Rcpp::List Env3;

// Call EnvelopeSort with the same arguments as the R orchestrator

if (verbose) {
  Rcpp::Rcout << "[EnvelopeSort] Entering: "
              << glmbayes::progress::timestamp_cpp()
              << "\n";
}

if (disp_grid_type == 1) {
  Env3 = EnvelopeSort(
      Rcpp::_["l1"]      = l1,
      Rcpp::_["l2"]      = l2,
      Rcpp::_["GIndex"]  = Env3_raw["GridIndex"],
      Rcpp::_["G3"]      = Env3_raw["thetabars"],
      Rcpp::_["cbars"]   = cbars,
      Rcpp::_["logU"]    = Env3_raw["logU"],
      Rcpp::_["logrt"]   = Env3_raw["logrt"],
      Rcpp::_["loglt"]   = Env3_raw["loglt"],
      Rcpp::_["logP"]    = logP_mat,   // <-- correct shape, untouched values
      Rcpp::_["LLconst"] = Env3_raw["LLconst"],
      Rcpp::_["PLSD"]    = Env3_raw["PLSD"],
      Rcpp::_["a1"]      = Env3_raw["a1"],
      Rcpp::_["E_draws"] = Env3_raw["E_draws"],
      Rcpp::_["lg_prob_factor"] = UB_list_new["lg_prob_factor"],
      Rcpp::_["UB2min"]         = UB_list_new["UB2min"]
    // ,  Rcpp::_["thetabar_const_base"] = UB_list_new["thetabar_const_base"],
    //   Rcpp::_["New_LL_Slope"]= UB_list_new["New_LL_Slope"],
    //   Rcpp::_["shape3_face"]= gamma_list_new["shape3_face"]
  );
}

if (disp_grid_type == 2) {
  Env3 = EnvelopeSort(
      Rcpp::_["l1"]      = l1,
      Rcpp::_["l2"]      = l2,
      Rcpp::_["GIndex"]  = Env3_raw["GridIndex"],
      Rcpp::_["G3"]      = Env3_raw["thetabars"],
      Rcpp::_["cbars"]   = cbars,
      Rcpp::_["logU"]    = Env3_raw["logU"],
      Rcpp::_["logrt"]   = Env3_raw["logrt"],
      Rcpp::_["loglt"]   = Env3_raw["loglt"],
      Rcpp::_["logP"]    = logP_mat,   // <-- correct shape, untouched values
      Rcpp::_["LLconst"] = Env3_raw["LLconst"],
      Rcpp::_["PLSD"]    = Env3_raw["PLSD"],
      Rcpp::_["a1"]      = Env3_raw["a1"],
      Rcpp::_["E_draws"] = Env3_raw["E_draws"],
      Rcpp::_["lg_prob_factor"] = UB_list_new["lg_prob_factor"],
      Rcpp::_["UB2min"]         = UB_list_new["UB2min"]
    // ,Rcpp::_["thetabar_const_base"] = UB_list_new["thetabar_const_base"],
    //   Rcpp::_["New_LL_Slope"]= UB_list_new["New_LL_Slope"],
    //   Rcpp::_["shape3_face"]= gamma_list_new["shape3_face"]
  );
}

if (Env3.containsElementNamed("sort_ok") && !Rcpp::as<bool>(Env3["sort_ok"])) {
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeSort] Using unsorted envelope (memory fallback).\n";
  }
  Env3 = Env3_raw;
  Env3["lg_prob_factor"] = UB_list_new["lg_prob_factor"];
  Env3["UB2min"]         = UB_list_new["UB2min"];
  if (disp_grid_type == 1) {
    Env3["thetabar_const_base"] = UB_list_new["thetabar_const_base"];
    Env3["New_LL_Slope"]        = UB_list_new["New_LL_Slope"];
    Env3["shape3_face"]         = gamma_list_new["shape3_face"];
  }
}
if (verbose) {
  Rcpp::Rcout << "[EnvelopeSort] Exiting: "
              << glmbayes::progress::timestamp_cpp()
              << "\n";
}


// --- Step 4: Update UB_list_new with reordered values from Env3 ---

UB_list_new["lg_prob_factor"] = Env3["lg_prob_factor"];
UB_list_new["UB2min"]         = Env3["UB2min"];

if(disp_grid_type==1){
  UB_list_new["thetabar_const_base"] = Env3["thetabar_const_base"];
  UB_list_new["New_LL_Slope"]        = Env3["New_LL_Slope"];
  gamma_list_new["shape3_face"]      = Env3["shape3_face"];
}



// Extract gamma_list, diagnostics, low, upp from disp_env_out
Rcpp::List diagnostics    = disp_env_out["diagnostics"];
double low = gamma_list_new["disp_lower"];
double upp = gamma_list_new["disp_upper"];

// --- Final return structure (matches R orchestrator exactly) ---
return Rcpp::List::create(
  Rcpp::_["Env"]        = Env3,
  Rcpp::_["gamma_list"] = gamma_list_new,
  Rcpp::_["UB_list"]    = UB_list_new,
  Rcpp::_["diagnostics"] = diagnostics,
  Rcpp::_["low"]        = low,
  Rcpp::_["upp"]        = upp
);

}

}
}