## Non-Zellner prior: Independent Normal-Gamma with a strongly shrunk diagonal
## covariance (0.001 * diag(diag(ps$Sigma))) instead of full Prior_Setup Sigma.
## Formerly in inst/examples/Ex_lmb.R as a "temporary non zellner test".

test_that("lmb: Independent Normal-Gamma with scaled diagonal Sigma (non-Zellner)", {
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)

  ps <- Prior_Setup(weight ~ group, gaussian())
  Sigma_non_zellner <- 0.001 * diag(diag(ps$Sigma))

  fit <- lmb(
    weight ~ group,
    dIndependent_Normal_Gamma(
      ps$mu,
      Sigma_non_zellner,
      shape = ps$shape_ING,
      rate  = ps$rate
    ),
    n = 500L,
    verbose = FALSE,
    use_parallel = FALSE
  )

  expect_s3_class(fit, "lmb")
  expect_equal(nrow(fit$coefficients), 500L)
  expect_equal(ncol(fit$coefficients), nrow(ps$mu))
  expect_true(all(is.finite(fit$coef.means)))
})

test_that("lmb: strongly anisotropic Independent Normal-Gamma prior does not trigger UB2 sign violations", {
  ## Deliberately anisotropic (high condition number) coefficient prior with
  ## a small sample size, mirroring the scenarios that used to trigger
  ## "Sign violation: UB2 < 0" errors under the old endpoint-only UB2_Min
  ## calculation in src/EnvelopeDispersionBuild.cpp::bound_ub2_over_dispersion
  ## (the endpoint-only shortcut is only exact when the coefficient prior is
  ## isotropic in the standardized K = Q^{-1/2} P Q^{-1/2} sense, e.g. a
  ## Zellner g-prior; see glmbayesCore's
  ## data-raw/README_ub2_rootfinding_fix.md and
  ## vignettes/Chapter-A07.Rmd Remark 5.5.7 / Claim 7 for the underlying
  ## theory and the exact root-finding fix ported here).
  set.seed(42)
  n <- 14L
  x1 <- rnorm(n)
  beta_true <- c(0.5, -1.2)
  y <- as.numeric(cbind(1, x1) %*% beta_true + rnorm(n, sd = 1.5))
  dat <- data.frame(y = y, x1 = x1)

  mu <- c(0, 0)
  Sigma_aniso <- diag(c(2000, 0.05))  # condition number 40000 in Sigma space

  for (rep_i in 1:10) {
    set.seed(100 + rep_i)
    fit <- lmb(
      y ~ x1,
      dIndependent_Normal_Gamma(mu, Sigma_aniso, shape = 3, rate = 2),
      data = dat,
      n = 300L,
      verbose = FALSE,
      use_parallel = FALSE
    )
    expect_s3_class(fit, "lmb")
    expect_equal(nrow(fit$coefficients), 300L)
    expect_true(all(is.finite(fit$coef.means)))
  }
})

test_that("Prior_Setup Gaussian calibration: shape from n_prior, E[sigma^2|y] = dispersion", {
  ctl <- c(4.17, 5.58)
  trt <- c(4.81, 4.17)
  group <- gl(2, 2, 4)
  weight <- c(ctl, trt)
  ps <- Prior_Setup(
    weight ~ group,
    gaussian(),
    pwt = 0.01
  )
  p <- ncol(ps$x)
  n_w <- ps$PriorSettings$n_effective
  n_prior <- ps$PriorSettings$n_prior
  S_marg <- ps$dispersion * (n_w - p)
  expect_equal(ps$shape, (n_prior + 1) / 2)
  post_mean_sigma2 <- (ps$rate + S_marg / 2) / (ps$shape + n_w / 2 - 1)
  expect_equal(post_mean_sigma2, ps$dispersion)
})
