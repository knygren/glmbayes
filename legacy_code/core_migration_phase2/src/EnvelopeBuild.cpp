// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include <RcppParallel.h>
#include "openclPort.h"
#include "progress_utils.h"
#include "R_interface.h"


using namespace Rcpp;
using namespace openclPort;
using namespace glmbayes::fam;
using namespace glmbayes::env;
using namespace glmbayes::progress;



/*
 EnvelopeBuild_c — envelope grid construction for models with Gaussian priors
 and log-concave likelihoods
 
 Purpose
 Construct an axis-aligned envelope grid and the auxiliary objects required
 by EnvelopeSet_Grid and the envelope sampler. This routine assumes a multivariate
 normal prior (precision-like matrix A) and a log-concave likelihood so that
 the (transformed) posterior is unimodal and well-behaved near the mode.
 The routine is family-agnostic otherwise: envelope placement uses only the
 posterior center (bStar) and curvature-like information (diag(A)).
 
 Lint construction (per-dimension cutpoints)
 - Lint is a 2 × l1 matrix; column i defines the central interval cutpoints
 for dimension i:
 Lint[0,i] := ℓ_{i,1} = θ*_i − 0.5 · ω_i
 Lint[1,i] := ℓ_{i,2} = θ*_i + 0.5 · ω_i
 where θ* = bStar (the posterior mode on the parameter scale used for the
 envelope) and ω is a spread parameter derived from the i-th diagonal of A.
 - Implementation details:
 a_2   = diag(A)            // curvature/precision per dimension
 omega = (sqrt(2) - exp(-1.20491 - 0.7321 * sqrt(0.5 + a_2))) / sqrt(1 + a_2)
 yy_1  = [1, 1]; yy_2 = [-0.5, 0.5]
 Lint  = yy_1 %*% transpose(bStar) + yy_2 %*% transpose(omega)
 The constants in the omega formula are empirical calibrations chosen to
 produce robust interval widths across typical precision regimes.
 
 Why Gaussian-prior + log-concave-likelihood matters
 - Log-concavity of the likelihood combined with a Gaussian prior implies the
 (unnormalized) log-posterior is concave around the mode; a single-mode
 assumption justifies axis-aligned intervals centered at the posterior mode.
 - Using diag(A) to set ω leverages local curvature information: larger
 precision (larger a_i) yields narrower ω_i, smaller a_i yields wider ω_i.
 - These choices ensure the central interval ℓ_{i,1}..ℓ_{i,2} captures the
 high-density mass along each axis while keeping tails handled by the
 left/right intervals.
 
 Gridpoint generation and Gridtype
 - G1 contains candidate per-dimension points {θ* - ω, θ*, θ* + ω}.
 - Gridtype controls whether a dimension uses the three-point set or only the
 mode:
 Gridtype 1: static threshold test using (1 + a_i) ≤ 2/√π → single-point
 Gridtype 2: dynamic selection via EnvelopeOpt(a_2, n) (cost-based)
 Gridtype 3: always three-point
 Gridtype 4: always single-point
 - Rationale: trade-off between build cost and sampling cost; with larger n
 or available parallelism, richer grids can be worthwhile.
 
 Full-grid expansion and cbars
 - Per-dimension G2 lists are combined with expand.grid to form the full set
 of grid locations (G3) and corresponding region codes (GIndex).
 - G4 = transpose(G3); cbars[j,i] = G4[j,i] − bStar[i] is the j-th component's
 offset from the mode for dimension i.
 - In EnvelopeSet_Grid these cbars shift Lint per-row to produce Down[j,i] and Up[j,i],
 the final bounds for truncated-normal evaluations.
 
 Numerical and modelling notes
 - The envelope construction uses only bStar and diag(A); any family-specific
 transforms (e.g., link functions, canonical transforms) should be applied
 upstream so bStar and A are expressed on the scale used by the envelope.
 - Log-concavity guarantees the mode is meaningful for centering intervals;
 if the likelihood is not log-concave the assumptions behind axis-aligned
 modal envelopes break down and the envelope builder may need alternative
 strategies.
 - Keep the Lint construction here as the single source of truth; if ω
 calibration constants are refined, update them only in this function.
 
 OpenCL hinting and verbosity
 - If compiled with OpenCL and use_opencl == true, the function may scale n
 by the detected core count to exploit parallel envelope optimization.
 - If OpenCL is unavailable at compile-time, the function logs a diagnostic
 (when verbose) and disables use_opencl, falling back to CPU.
 
 Outputs used downstream
 - G2: per-dimension gridpoint lists
 - G3/G4: full grid and transpose
 - GIndex: integer region codes per grid component and dimension
 - Lint: two unshifted cutpoints per dimension (lower, upper)
 - cbars: per-component offsets from bStar used to shift Lint in EnvelopeSet_Grid
 - These objects enable EnvelopeSet_Grid to compute per-dimension truncated-normal
 probabilities U_{j,i} and their log-sums logP[j], suitable for envelope
 weighting and acceptance tests under the Gaussian-prior + log-concave
 likelihood assumption.
 */









namespace glmbayes {

namespace env {
List EnvelopeBuild(NumericVector bStar,
                       NumericMatrix A, /// Diagonal Precision Matrix for Adjusted Likelihood Function
                       NumericVector y, 
                       NumericMatrix x,
                       NumericMatrix mu,
                       NumericMatrix P, /// Part of the prior precision matrix that is shifted to the likelihood
                       NumericVector alpha,
                       NumericVector wt,
                       std::string family,
                       std::string link,
                       int Gridtype, 
                       int n,
                       int n_envopt,
                       bool sortgrid,
                       bool use_opencl ,        // Enables OpenCL acceleration during envelope construction
                       bool verbose             // Enables diagnostic output
                       
){
  
  if (n_envopt < 0) {
    n_envopt = n;
  }  
  
  // Handle OpenCL availability at compile/runtime
  // (fall back to CPU if requested but not supported)
  
  
  
#ifdef USE_OPENCL
  // OpenCL support detected at compile time — proceed as requested
#else
  if (use_opencl) {
    if (verbose) {
      Rcpp::Rcout << "[NOTE] OpenCL support was not detected during configuration.\n";
      Rcpp::Rcout << "       Disabling use_opencl and falling back to CPU implementation.\n";
      Rcpp::Rcout << "       To enable OpenCL, install an OpenCL SDK and ensure CL/cl.h is discoverable.\n";
      Rcpp::Rcout << "       You may need to set OPENCL_HOME or add the SDK to your system PATH.\n";
    }
    use_opencl = false;
  }
#endif
  
  // Print call parameters if verbose
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] Arguments:\n"
                << "  [Gridtype]   " << Gridtype   << "\n"
                << "  [n]          " << n          << "\n"
                << "  [use_opencl] " << use_opencl << "\n"
                << "  [sortgrid]   " << sortgrid   << "\n";
  }
  
  
  // Basic setup: dimensions, Armadillo views, and working vectors
  // l1 = number of parameters, k = number of predictors
  
  
  // int progbar=0;
  
  int l1 = A.nrow(), k = A.ncol();
  arma::mat A2(A.begin(), l1, k, false);
  arma::colvec bStar_2(bStar.begin(), bStar.size(), false);
  
  
  NumericVector a_1(l1);
  arma::vec a_2(a_1.begin(), a_1.size(), false);
  
  NumericVector xx_1(3, 1.0);
  NumericVector xx_2=NumericVector::create(-1.0,0.0,1.0);
  NumericVector yy_1(2, 1.0);
  NumericVector yy_2=NumericVector::create(-0.5,0.5);
  NumericMatrix G1(3,l1);
  NumericMatrix Lint1(2,l1);
  arma::mat G1b(G1.begin(), 3, l1, false);
  arma::mat Lint(Lint1.begin(), 2, l1, false);
  
  arma::colvec xx_1b(xx_1.begin(), xx_1.size(), false);
  arma::colvec xx_2b(xx_2.begin(), xx_2.size(), false);
  arma::colvec yy_1b(yy_1.begin(), yy_1.size(), false);
  arma::colvec yy_2b(yy_2.begin(), yy_2.size(), false);
  // List G2(a_1.size());
  // List GIndex1(a_1.size());
  Rcpp::Function EnvelopeOpt = glmbayes_R::r_envelope_opt();
  Rcpp::Function expGrid("expand.grid");
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function EnvSort = glmbayes_R::r_envelope_sort();
  
  // int i;  
  
  
  // Construct per-dimension tangent points (G1) and linear intercepts (Lint)
  // using diagonal precisions and offsets (omega)
  
  a_2=arma::diagvec(A2);
  arma::vec omega=(sqrt(2)-arma::exp(-1.20491-0.7321*sqrt(0.5+a_2)))/arma::sqrt(1+a_2);
  G1b=xx_1b*arma::trans(bStar_2)+xx_2b*arma::trans(omega);
  Lint=yy_1b*arma::trans(bStar_2)+yy_2b*arma::trans(omega);
  
  
  // Call EnvelopeSize to determine grid structure and expected draws
  
  Rcpp::List size_info = EnvelopeSize(a_2, G1,
                                      Gridtype,
                                      n,
                                      n_envopt,
                                      use_opencl,
                                      verbose);
  
  
  // Unpack results
  Rcpp::List G2        = size_info["G2"];
  Rcpp::List GIndex1   = size_info["GIndex1"];
  double E_draws       = Rcpp::as<double>(size_info["E_draws"]);
  Rcpp::NumericVector gridindex = size_info["gridindex"];
  
  // 
  // Expand grid indices and candidate points (GIndex, G3, G4)
  // l2 = total number of grid combinations
  
  
  NumericMatrix GIndex=asMat(expGrid(GIndex1));
  int l2=GIndex.nrow();
  NumericMatrix G3=asMat(expGrid(G2));
  NumericMatrix G4(G3.ncol(),G3.nrow());
  

  if (verbose) {
  Rcpp::Rcout << "[EnvelopeBuild] Grid size (l2) = "
              << glmbayes::progress::format_int_with_commas(l2)
              << "\n";
}
  
  

  
  arma::mat G3b(G3.begin(), G3.nrow(), G3.ncol(), false);
  arma::mat G4b(G4.begin(), G4.nrow(), G4.ncol(), false);
  
  G4b=trans(G3b);
  
  // Allocate containers for evaluation results (cbars, NegLL, logP, etc.)
  
  NumericMatrix cbars(l2,l1);
  NumericMatrix Up(l2,l1);
  NumericMatrix Down(l2,l1);
  NumericMatrix logP(l2,2);
  NumericMatrix logU(l2,l1);
  NumericMatrix loglt(l2,l1);
  NumericMatrix logrt(l2,l1);
  NumericMatrix logct(l2,l1);
  
  NumericMatrix LLconst(l2,1);
  NumericVector NegLL(l2);    
  
  NumericVector NegLL_Alt(l2);    /// Temporary
  
  
  arma::mat cbars2(cbars.begin(), l2, l1, false); 
  arma::mat cbars3(cbars.begin(), l2, l1, false); 
  
  // Note: NegLL_2 only added to allow for QC printing of results 
  
  arma::colvec NegLL_2(NegLL.begin(), NegLL.size(), false);
  
  // Call EnvelopeEval to compute negative log-likelihood and gradients
  // at each grid point

  if (verbose) {
    
    Rcpp::Rcout << "[EnvelopeBuild:EnvelopeEval] Entering: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                << glmbayes::progress::timestamp_cpp()
                << "\n";
  }
  
  
  Rcpp::List eval_info = EnvelopeEval(G4, y, x, mu, P, alpha, wt,
                                      family, link, use_opencl, verbose);
  
  
  if (verbose) {
    
    Rcpp::Rcout << "[EnvelopeBuild:EnvelopeEval] Exiting: "
                << glmbayes::progress::timestamp_cpp()
                << "\n";
  }
  
  
  // Copy results into cbars/NegLL structures used downstream
  
  NegLL = eval_info["NegLL"];
  cbars2 = Rcpp::as<arma::mat>(eval_info["cbars"]);
  
  // Do a temporary correction here cbars3 should point to correct memory
  // See if this sets cbars
  
  cbars3=cbars2;
  

  if (verbose) {
    
    Rcpp::Rcout << "[EnvelopeBuild:EnvelopeSet_Grid] Entering: "
                << glmbayes::progress::timestamp_cpp()
                << "\n";
  }
  
  // Set Grid
  
  EnvelopeSet_Grid_C2_pointwise(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);
  
  if (verbose) {
    
    Rcpp::Rcout << "[EnvelopeBuild:EnvelopeSet_Grid] Exiting: "
                   << glmbayes::progress::timestamp_cpp()
                   << "\n";
  }
  
  

  // Set LOG P
  
  if (verbose) {
    
    Rcpp::Rcout << "[EnvelopeBuild:EnvelopeSet_LogP] Entering: "
                << glmbayes::progress::timestamp_cpp()
                << "\n";
  }
  
  
  setlogP_C2(logP,NegLL,cbars,G3,LLconst);
  
  if (verbose) {
    
    Rcpp::Rcout << "[EnvelopeBuild:EnvelopeSet_LogP] Exiting: "
              << glmbayes::progress::timestamp_cpp()
              << "\n";
  }
  
  // Normalize probabilities (PLSD) for likelihood subgradient densities
  
  NumericMatrix::Column logP2 = logP( _, 1);
  double  maxlogP=max(logP2);
  
  NumericVector PLSD=exp(logP2-maxlogP);
  
  double sumP=sum(PLSD);
  
  PLSD=PLSD/sumP;
  
  
  bool bad = false;
  
  for (int i = 0; i < PLSD.size(); ++i) {
    double v = PLSD[i];
    if (!R_finite(v) || v <= 0.0) {
      Rcout << "[ERROR] Invalid PLSD at index " << i
            << "  value=" << v
            << "  logP2=" << logP2[i]
            << "  maxlogP=" << maxlogP
            << "\n";
      bad = true;
      break;
    }
  }
  
  if (bad) {
    // --- Design-matrix diagnostics (unweighted) ---
    arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
    arma::vec s = arma::svd(X);
    
    double kappa_X   = s.max() / s.min();                         // ~ kappa(X)
    double kappa_XtX = (s.max() * s.max()) / (s.min() * s.min()); // ~ kappa(crossprod(X))
    
    double maxlogP = max(logP2);
    double minlogP = min(logP2);
    double spread  = maxlogP - minlogP;
    
    Rcpp::stop(
      "[EnvelopeBuild] PLSD construction failed for this model.\n"
      "  - kappa(X)           ≈ " + std::to_string(kappa_X)   + "\n"
      "  - kappa(t(X) %*% X)  ≈ " + std::to_string(kappa_XtX) + "\n"
      "  - PLSD log-weights span " + std::to_string(spread) + " on the log scale.\n"
      "Interpretation:\n"
      "  The design matrix used in the standardized model is highly collinear\n"
      "  (very large kappa), so X'X has extremely uneven curvature. In combination\n"
      "  with the chosen prior, this creates posterior directions that are either\n"
      "  very sharp or nearly flat. The fixed PLSD partition from Nygren & Nygren\n"
      "  (2006) cannot remain tight under this curvature, so the envelope becomes\n"
      "  extremely loose and the PLSD weights blow up.\n"
    );
  }
  
  // Optionally sort grid for efficiency if sortgrid = TRUE
  
  if(sortgrid==true){
    
    if (verbose) {
      
      Rcpp::Rcout << "[EnvelopeBuild:EnvelopeSort] Entering: "
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
    }

    // C++ sort commented out (slower than R; R result used downstream)
    // if (verbose) {
    //   Rcpp::Rcout << "[EnvelopeBuild:EnvelopeSort_cpp] Entering: "
    //               << glmbayes::progress::timestamp_cpp()
    //               << "\n";
    // }
    // (void) EnvelopeSort_cpp(
    //   l1, l2,
    //   GIndex, G3, cbars,
    //   logU, logrt, loglt, logP,
    //   LLconst, PLSD, a_1, E_draws
    // );
    // if (verbose) {
    //   Rcpp::Rcout << "[EnvelopeBuild:EnvelopeSort_cpp] Exiting: "
    //               << glmbayes::progress::timestamp_cpp()
    //               << "\n";
    // }

    Rcpp::List outlist = EnvSort(l1, l2, GIndex, G3, cbars, logU, logrt, loglt, logP, LLconst, PLSD, a_1, E_draws);

    if (outlist.containsElementNamed("sort_ok") && !Rcpp::as<bool>(outlist["sort_ok"])) {
      if (verbose) {
        Rcpp::Rcout << "[EnvelopeBuild:EnvelopeSort] Using unsorted envelope (memory fallback).\n";
      }
      outlist = Rcpp::List::create(
        Rcpp::Named("GridIndex") = GIndex,
        Rcpp::Named("thetabars") = G3,
        Rcpp::Named("cbars")    = cbars,
        Rcpp::Named("logU")     = logU,
        Rcpp::Named("logrt")    = logrt,
        Rcpp::Named("loglt")    = loglt,
        Rcpp::Named("LLconst")  = LLconst,
        Rcpp::Named("logP")     = logP(_, 0),
        Rcpp::Named("PLSD")     = PLSD,
        Rcpp::Named("a1")       = a_1,
        Rcpp::Named("E_draws")  = E_draws
      );
    }
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild:EnvelopeSort] Exiting:"
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
    }
    return(outlist);
    
  }
  
  // Return assembled envelope components as a list
  
  
  return Rcpp::List::create(Rcpp::Named("GridIndex")=GIndex,
                            Rcpp::Named("thetabars")=G3,
                            Rcpp::Named("cbars")=cbars,
                            Rcpp::Named("logU")=logU,
                            Rcpp::Named("logrt")=logrt,
                            Rcpp::Named("loglt")=loglt,
                            Rcpp::Named("LLconst")=LLconst,
                            Rcpp::Named("logP")=logP(_,0),
                            Rcpp::Named("PLSD")=PLSD,
                            Rcpp::Named("a1")=a_1,
                            Rcpp::Named("E_draws")=E_draws
  );
  
}

}
}