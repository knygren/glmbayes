# glmbayes

<<<<<<< HEAD
![GitHub release (latest by date)](https://img.shields.io/github/v/release/knygren/glmbayes?label=version)
![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/knygren/glmbayes/R-CMD-check.yaml?label=R%20CMD%20Check)

**glmbayes** provides independent and identically distributed (iid) samples for Bayesian Generalized Linear Models (GLMs), serving as a Bayesian counterpart to R’s classical `glm()` function. It supports Gaussian, Poisson, Binomial, and Gamma families using log-concave likelihoods and leverages accept-reject sampling via likelihood subgradients (Nygren & Nygren, 2006).

This package is currently in beta (`v0.1.0-beta`). For details on recent updates and planned enhancements, see the [`NEWS.md`](https://github.com/knygren/glmbayes/blob/main/NEWS.md).

## 📦 Installation

```r
# Install the beta release from GitHub
devtools::install_github("knygren/glmbayes@v0.1.0-beta")
```

## 🧪 Minimal Working Example

```r
library(glmbayes)

# Setup prior
ps <- Prior_Setup(counts ~ outcome + treatment)

=======
**glmbayes** provides independent and identically distributed (iid) samples for Bayesian Generalized Linear Models (GLMs), serving as a Bayesian counterpart to R’s classical `glm()` function. It supports Gaussian, Poisson, Binomial, and Gamma families using log-concave likelihoods and leverages accept-reject sampling via likelihood subgradients (Nygren & Nygren, 2006).

## 📦 Installation

```r
# Install the beta release from GitHub
devtools::install_github("knygren/glmbayes@v0.1.0-beta")
```

## 🧪 Minimal Working Example

```r
library(glmbayes)

# Setup prior
ps <- Prior_Setup(counts ~ outcome + treatment)

>>>>>>> f98b84a9ca45edd79584a165efa5b3000c8fabc1
# Fit Bayesian GLM
fit <- glmb(counts ~ outcome + treatment,
            family = poisson(),
            pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma))

# Summarize results
summary(fit)
```

## ✨ Feature Highlights

- S3 interface mirroring the structure of base `glm()`
- Accept-reject sampling for log-concave likelihoods
- Samplers for both fixed and variable dispersion
- Vignette-based comparisons with classical GLM estimates
- Modular prior setup and checking tools

## 📘 Supported Families & Links

- **Gaussian** (identity)
- **Poisson / Quasi-Poisson** (log)
- **Binomial / Quasi-Binomial** (logit, probit, cloglog)
- **Gamma** (log)

All supported models feature log-concave likelihoods, enabling efficient iid sampling via enveloping functions and subgradient-based accept-reject algorithms.

## 📚 Vignettes & Demos

Use `demo()` to explore built-in examples for supported families and links:

```r
demo("Gaussian_identity")
demo("Poisson_log")
demo("Binomial_logit")
```

To view vignette documentation:

```r
browseVignettes("glmbayes")
```

Topics include comparisons with `glm()` outputs, two-block Gibbs sampling strategies, and handling of dispersion parameters.

## 🧠 Methodology

Sampling follows the framework from Nygren & Nygren (2006), using likelihood subgradients to construct enveloping functions for the posterior distribution. When the posterior is approximately normal, the expected number of draws per acceptance is bounded. Dispersion parameters can be estimated separately using `rglmbdisp()` for Gaussian and Gamma families.

## 🚧 Limitations

- Non-log-concave likelihoods are not currently supported
- Dispersion estimation requires a second sampler (`rglmbdisp()`)
- Hierarchical modeling and GPU support remain experimental

## 📈 Future Plans

- Unified two-block Gibbs sampling with joint dispersion updates
- GPU acceleration for high-dimensional models
- Full CRAN submission and expanded vignette documentation