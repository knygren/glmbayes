core <- c(
  "EnvelopeBuild", "EnvelopeCentering", "EnvelopeDispersionBuild",
  "EnvelopeEval", "EnvelopeOpt", "EnvelopeOrchestrator", "EnvelopeSetGrid",
  "EnvelopeSetLogP", "EnvelopeSize", "EnvelopeSort", "Prior_Check",
  "compute_gaussian_prior", "dBeta", "dGamma", "dIndependent_Normal_Gamma",
  "dNormal", "dNormal_Gamma", "glmb.wfit",
  "glmb_Standardize_Model", "glmbfamfunc", "multi_prior_setup",
  "pfamily", "pinvgamma_ct", "pnorm_ct", "qinvgamma_ct", "rBeta_prior",
  "rBeta_reg", "rGamma_Conjugate_prior", "rGamma_Conjugate_reg",
  "rGamma_prior", "rGamma_reg", "rIndepNormalGammaReg_std",
  "rIndependent_Normal_Gamma_prior", "rNormalGLM_std", "rNormalGamma_reg",
  "rNormal_Gamma_prior", "rNormal_prior", "rNormal_reg", "rNormal_reg.wfit",
  "rgamma_ct", "rglmb", "rindepNormalGamma_reg", "rinvgamma_ct", "rlmb",
  "rnorm_ct", "simfunction"
)
lines <- c(
  "## Type (2): re-exports from glmbayesCore (Phase 2 backend migration).",
  "## Archived implementations: legacy_code/core_migration_phase2/R/.",
  "",
  "NULL"
)
for (fn in core) {
  lines <- c(
    lines,
    sprintf(
      "#' @inherit glmbayesCore::%s return title description params details examples references format note",
      fn
    ),
    "#' @export",
    sprintf("%s <- glmbayesCore::%s", fn, fn),
    ""
  )
}
out <- "C:/Rpackages/glmbayes/R/reexports-core-engine.R"
writeLines(lines, out)
message("Wrote ", out, " (", length(lines), " lines)")
