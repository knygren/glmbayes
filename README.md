# glmbayes

![CRAN status](https://www.r-pkg.org/badges/version/glmbayes)
![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/glmbayes)
![Monthly downloads](https://cranlogs.r-pkg.org/badges/glmbayes)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/knygren/glmbayes?label=version)
![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/knygren/glmbayes/R-CMD-check.yaml?label=R%20CMD%20Check)

glmbayes provides independent and identically distributed (iid) samples for Bayesian Generalized Linear Models (GLMs).
Its primary interface, glmb(), serves as a Bayesian analogue to R's glm() function, supporting Gaussian, Poisson,
Binomial, and Gamma families under log-concave likelihoods. Sampling for most models is performed using accept-reject
methods based on likelihood subgradients (Nygren and Nygren, 2006). For Gaussian models, the package also includes
lmb(), a Bayesian counterpart to R's lm().

The package includes a rich set of supporting tools for prior specification, model diagnostics, and method functions
that mirror those for lm() and glm(). Most functions are extensively documented, and a comprehensive set of vignettes
are available to guide users through the package's capabilities.

The current **CRAN release is version 0.9.7**
([CRAN](https://CRAN.R-project.org/package=glmbayes)).
The [GitHub](https://github.com/knygren/glmbayes) repository holds the source; [R-Universe](https://knygren.r-universe.dev/glmbayes) builds binaries from it.
See [NEWS.md](https://github.com/knygren/glmbayes/blob/main/NEWS.md) for changes.

## Installation

**CRAN (release 0.9.7)**

```r
install.packages("glmbayes")
```

**GitHub / R-Universe** (install from both CRAN and R-Universe repositories if you want R-Universe binaries or faster mirrors):

```r
install.packages("glmbayes",
                 repos = c("https://cloud.r-project.org",
                           "https://knygren.r-universe.dev"))
```

Prebuilt binaries from CRAN (0.9.7) and R-Universe are built **without OpenCL GPU
support**. For the CRAN release, OpenCL requires installing **from source** on a
system with OpenCL development files available. To set up GPU acceleration, follow

**Chapter 16 — Large models: GPU acceleration using OpenCL**
https://knygren.r-universe.dev/articles/glmbayes/Chapter-16.html

## Minimal Working Example

    library(glmbayes)

    # Dobson (1990), p. 93: Randomized Controlled Trial
    counts <- c(18,17,15,20,10,20,25,13,12)
    outcome <- gl(3,1,9)
    treatment <- gl(3,3)
    print(d.AD <- data.frame(treatment, outcome, counts))

    ## Classical glm
    glm.D93 <- glm(counts ~ outcome + treatment,
                   family = poisson())

    ## Bayesian glmb
    # Step 1: Set up prior
    ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())
    mu <- ps$mu
    V  <- ps$Sigma

    # Step 2: Fit using glmb
    glmb.D93 <- glmb(counts ~ outcome + treatment,
                     family = poisson(),
                     pfamily = dNormal(mu = mu, Sigma = V))

    summary(glmb.D93)

## Supported families, links, and pfamilies

As with `glm()`, models are defined by a formula for the linear predictor and a `family()` describing the likelihood and 
link. In addition, `glmb()` requires a **pfamily** object specifying the prior.

### Priors on regression coefficients

The primary table below covers priors on the regression coefficients **β**. The standard prior for
all families is `dNormal`. The conjugate priors `dBeta` and `dGamma(Inv_Dispersion = FALSE)` provide
closed-form IID posterior draws for intercept-only models with an identity link.

| Likelihood family           | Link functions                         | Compatible pfamilies (coefficient priors)                            |
|-----------------------------|----------------------------------------|----------------------------------------------------------------------|
| Gaussian                    | identity                               | dNormal, dNormal_Gamma, dIndependent_Normal_Gamma                    |
| Poisson / Quasi-Poisson     | log                                    | dNormal                                                              |
| Poisson                     | identity *(intercept-only)*            | dGamma(Inv_Dispersion = FALSE) — conjugate Gamma–Poisson rate prior  |
| Binomial / Quasi-Binomial   | logit, probit, cloglog                 | dNormal                                                              |
| Binomial                    | identity *(intercept-only)*            | dBeta — conjugate Beta–Binomial probability prior                    |
| Gamma                       | log                                    | dNormal                                                              |
| Gamma                       | identity *(intercept-only)*            | dGamma(Inv_Dispersion = FALSE) — conjugate Gamma–Gamma rate prior    |

`dNormal_Gamma` and `dIndependent_Normal_Gamma` also model precision jointly with the coefficients;
see the precision/dispersion table below.

### Priors on precision / dispersion

`dGamma(Inv_Dispersion = TRUE)` (the default when `Inv_Dispersion` is omitted) places a Gamma prior
on the inverse dispersion **1/φ** with the regression coefficients **β** held fixed. This is the
precision prior used in Gibbs sampling steps for dispersion estimation.

| Likelihood family | Link     | Compatible pfamilies (precision prior)      |
|-------------------|----------|---------------------------------------------|
| Gaussian          | identity | dGamma — prior on 1/σ² (precision)          |
| Gamma             | log      | dGamma — prior on 1/φ (shape / dispersion)  |

`dNormal_Gamma` and `dIndependent_Normal_Gamma` model **β** and precision jointly in a single
conjugate step, avoiding the need for a separate Gibbs precision update.

### Prior_Setup

For a default, data‑aligned prior using the same formula and family as `glm()`, call `Prior_Setup(formula, family, data = ..., ...)`. 
The returned list includes default settings for the following:

- **mu**, **Sigma** — Zellner‑style normal prior components for use with most priors  
- **Additional Gaussian‑specific calibration components**:  
  - `dispersion` for use with the `dNormal()` prior (gaussian and Gamma families)
  - `Sigma_0`, `shape` and `rate` for use with the `dNormal_Gamma()` prior  
  - `shape_ING` and `rate` for use with `dIndependent_Normal_Gamma()` prior 
  - `shape`, `rate_gamma` and `coefficients` for use with the `dGamma()` precision prior  
- **Conjugate prior calibration components** (intercept-only models):  
  - `conj_beta` (`shape1`, `shape2`, `beta`) for use with `dBeta()` (Binomial/identity)  
  - `conj_poisson` (`shape`, `rate`, `beta`) for use with `dGamma(Inv_Dispersion = FALSE)` (Poisson/identity)  

Optional arguments adjust prior weight, centering, and related settings (see the function help and vignette Chapter 04).

### Typical Prior_Setup wiring

Assuming `ps <- Prior_Setup(...)`:

- **Non‑Gaussian families (log/logit/probit/cloglog links):**  
  Use `dNormal(mu = ps$mu, Sigma = ps$Sigma)`.  
  (For Gamma GLMs, also supply `dispersion` from the fitted GLM or from `ps`; see `example("glmb")`.)

- **Binomial — conjugate Beta prior (identity link, intercept-only):**  
  Use `dBeta(shape1 = ps$conj_beta$shape1, shape2 = ps$conj_beta$shape2, beta = ps$conj_beta$beta)`.

- **Poisson — conjugate Gamma rate prior (identity link, intercept-only):**  
  Use `dGamma(shape = ps$conj_poisson$shape, rate = ps$conj_poisson$rate, beta = ps$conj_poisson$beta, Inv_Dispersion = FALSE)`.

- **Gaussian — normal prior with known dispersion:**  
  Use `dNormal(mu = ps$mu, Sigma = ps$Sigma, dispersion = ps$dispersion)`.

- **Gaussian — conjugate Normal–Gamma:**  
  Use `dNormal_Gamma(mu = ps$mu, Sigma_0 = ps$Sigma_0, shape = ps$shape, rate = ps$rate)`.

- **Gaussian — independent Normal–Gamma:**  
  Use `dIndependent_Normal_Gamma(mu = ps$mu, Sigma = ps$Sigma, shape = ps$shape_ING, rate = ps$rate)`.

- **Gaussian / Gamma — precision prior (coefficients fixed, for Gibbs):**  
  With `rate_dg <- if (!is.null(ps$rate_gamma)) ps$rate_gamma else ps$rate`, use  
  `dGamma(shape = ps$shape, rate = rate_dg, beta = ps$coefficients)`.

The default priors have limiting behaviors that produce estimates resembling classical estimates as priors get weak 
(see documentation and vignettes for details).

All supported models have log‑concave likelihoods, enabling efficient iid sampling via enveloping functions
and subgradient‑based accept–reject algorithms, especially for models lacking standard iid samplers. 


## Examples and Demos

Use `example()` and `demo()` to explore built-in examples and demos for supported families and links:

    ## Bayesian linear regression
    example("lmb")

    ## Bayesian generalized linear models
    example("glmb")

    ## Beta-Binomial conjugacy: dBeta() prior; Bechdel test (requires bayesrules)
    ## See also: vignette("Chapter-02-S03", package = "glmbayes")
    demo("Ex_12_BetaBinomial")

    ## Gamma-Poisson conjugacy: dGamma(Inv_Dispersion=FALSE); bike counts + heart
    ## transplant mortality (requires bayesrules; Appendix A requires LearnBayes)
    ## See also: vignette("Chapter-02-S04", package = "glmbayes")
    demo("Ex_13_GammaPoisson")

    ## Predictions for fitted glmb objects (newdata, type, etc.)
    example("predict.glmb")

    ## Deviance residuals and simulate() for posterior predictive checks (menarche)
    example("residuals.glmb")

    ## Two-block Gibbs sampler compared with iid sampling (linear model)
    example("rlmb")

    ## Default prior specification using Prior_Setup
    example("Prior_Setup")

    ## Matrix-input GLM example with an informative prior
    example("rglmb")

    ## Two-step Boston example: estimates and summarizes models with unknown
    ## dispersion using dGamma priors via rGamma_reg, rglmb, rlmb, glmb, and lmb
    example("summary.rGamma_reg")

    ## High-dimensional Gaussian model (14 predictors) with GPU acceleration (requires OpenCL)
    example("Boston_centered")

    ## High-dimensional binomial model (14 predictors) with GPU acceleration (requires OpenCL)
    example("Cleveland")

    ## Hierarchical linear model (Rubin/Gelman 8-schools) via rlmb
    demo("Ex_07_Schools")

    ## Hierarchical generalized linear model (Poisson BikeSharing) via rglmb
    demo("Ex_09_BikeSharingPoisson")

    ## Detailed simulation pipeline for rNormalGLM models (JASA 2006; Vignette Chapter A05)
    example("rNormalGLM_std")

    ## Detailed simulation pipeline for rIndepNormalGammaReg models (Vignette Chapter A07)
    example("rIndepNormalGammaReg_std")

## Methodology

For generalized linear models where well known sampling methods are unavailable, sampling follows the
framework from Nygren and Nygren (2006), using likelihood subgradients to construct enveloping functions for
the posterior distribution. When the posterior is approximately normal, the expected number of draws per
acceptance is bounded as per that paper and as discussed in our vignettes.
Dispersion can be sampled via `rGamma_reg()` (standalone) or jointly with coefficients via
`rNormalGamma_reg()` and `rindepNormalGamma_reg()`.

## GPU Acceleration Using OpenCL

The implemented algorithms tend to have acceptable performance on CPUs up to around 10-14 dimensions.
For larger models, the envelope construction is embarrassingly parallel. To accelerate envelope construction
in such cases, the package provides optional GPU acceleration using OpenCL. This requires that users have
GPU enabled machines and an OpenCL installation. These features are discussed in more detail in two of
our vignettes.

## Vignettes

The glmbayes package includes a comprehensive set of vignettes organized into five major parts.
These vignettes guide users from introductory material through applied modeling, advanced topics,
and the underlying simulation methods that support the package.

### Part 1: An Introduction

Overview of the package, its design philosophy, single-parameter conjugate models, and the basic workflow for fitting Bayesian linear and generalized linear models.

- **Chapter 00 - Introduction**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-00.html

- **Chapter 01 - Getting Started with glmbayes**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-01.html

- **Chapter 02 — Conjugate inference for single parameters** (S01–S05)  
  Start with [Chapter 02-S01](https://knygren.r-universe.dev/articles/glmbayes/Chapter-02-S01.html); then S02 (Normal–Normal), S03 (Beta–Binomial), S04 (Gamma–Poisson), S05 (Gamma–Gamma).

### Part 2: Bayesian regression models

These chapters focus on Bayesian **linear** regression (Gaussian family). Topics include **`lmb()`** fitting, **`Prior_Setup()`**, posterior predictive checks (**bayesplot**), deviance residuals and model summaries, **bayestestR**-style summaries, and the bridge to Bayesian GLMs in Part 3.

- **Chapter 03 — Estimating Bayesian linear models**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-03.html

- **Chapter 04 — Tailoring priors — leveraging the Prior_Setup function**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-04.html

- **Chapter 05 — Model predictions and posterior predictive checks (+ bayesplot `ppc_*`)**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-05.html

- **Chapter 06 — Deviance residuals, model statistics and posterior inference (+ bayestestR)**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-06.html

### Part 3: Generalized Linear Models
This part presents Bayesian GLMs across the major likelihood families, including binomial,
quasi-binomial, Poisson, quasi-Poisson, and Gamma models. It covers model specification,
link functions, log-concavity, diagnostics, interpretation of posterior results, and tooling
(**bayesplot**, **bayestestR**) for visualization and summaries.

- **Chapter 07 — Foundations of GLMs — families, links, and log-concave likelihoods**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-07.html

- **Chapter 08 — Estimating Bayesian generalized linear models**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-08.html

- **Chapter 09 — Models for the Binomial family**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-09.html

- **Chapter 10 — Models for the Poisson family**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-10.html

- **Chapter 11 — Models for the Gamma family**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-11.html

- **Chapter 12 — Visualizing posteriors with bayesplot**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-12.html

- **Chapter 13 — Bayesian inference and decision making with bayestestR**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-13.html

### Part 4: Advanced Topics
These chapters explore more complex modeling scenarios and computational strategies, such as
informative priors, two-block Gibbs sampling, linear and generalized linear mixed-effects models,
models with unknown dispersion parameters, and large-scale model fitting using GPU acceleration
using OpenCL.

- **Chapter 14 — Informative priors — centering and differential prior weights**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-14.html

- **Chapter 15 — Estimating models with unknown dispersion parameters**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-15.html

- **Chapter 16 — Large models: GPU acceleration using OpenCL**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-16.html

- **Chapter 17 — Linear mixed-effects models**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-17.html

- **Chapter 18 — Generalized linear mixed-effects models**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-18.html

### Part 5: Simulation Methods and Technical Implementation
This part documents the mathematical and algorithmic foundations of the package. Topics include
estimation procedures, likelihood subgradient densities, envelope construction, accept-reject
sampling, and technical reports on sampler design including implementation aspects for GPU acceleration using
OpenCL.

- **Chapter A01 - A detailed overview of the glmbayes package**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A01.html

- **Chapter A02 - Overview of Estimation Procedures**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A02.html

- **Chapter A03 - Methods Available in glmbayes**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A03.html

- **Chapter A04 - Directional Tail Diagnostics for Prior-Posterior Disagreement**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A04.html

- **Chapter A05 - Simulation Methods - Likelihood Subgradient Densities**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A05.html

- **Chapter A06 - Accept-Reject Sampling for Dispersion in Gamma Regression**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A06.html

- **Chapter A07 - Accept-Reject Sampling for gaussian Regression models with independent normal-gamma priors**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A07.html

- **Chapter A08 - Overview of Envelope Related Functions**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A08.html

- **Chapter A09 - Parallel Sampling Implementation using RcppParallel**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A09.html

- **Chapter A10 - Accelerated EnvelopeBuild Implementation using OpenCL**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A10.html

- **Chapter A11 - Implementation Companion for Independent Normal-Gamma**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A11.html

- **Chapter A12 - Technical Derivations for Priors Returned by `Prior_Setup()`**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A12.html


Together, these vignettes form a comprehensive reference that supports users at all levels, 
from first-time Bayesian GLM users to researchers interested in the mathematical and computational
details behind the samplers.

## Feature Highlights

- S3 interface mirroring the structure of base glm()
- Posterior predictive checks via `pp_check()` from the 'bayesplot' package for fitted `glmb` objects
- Accept-reject sampling for log-concave likelihoods
- Samplers for both fixed and variable dispersion
- Extensive vignettes to guide users through the package's capabilities
- Modular prior setup function

## Limitations

- Non-log-concave likelihoods are not currently supported

## Future Plans

The main near-term goal is to migrate much of the sampling back end — the
**`src/*.cpp`** code and related R-layer simulation/envelope functions — into
[**glmbayesCore**](https://github.com/knygren/glmbayesCore), leaving **glmbayes**
focused on the user-facing **`glmb()`** / **`lmb()`** interface, priors,
diagnostics, and documentation. **glmbayesCore** already hosts shared CPU/OpenCL
nmath and several core samplers; the migration will extend that split.

Further performance and algorithm work:

- Poisson speed (OpenCL and simulation): Precompute the log-factorial term `log(y!)`
  once per observation and reuse it in both OpenCL envelope construction and
  accept-reject simulation, since it depends only on the response, to reduce
  redundant `lgamma` evaluation and improve performance for large Poisson models.
- Grid selection (simulation): Precompute cumulative PLSD and use inverse CDF
  sampling (e.g. binary search) to select the grid component per candidate
  instead of scanning PLSD, improving the simulation loop when many candidates
  are evaluated.