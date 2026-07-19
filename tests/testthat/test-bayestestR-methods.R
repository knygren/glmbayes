## Tests for bayestestR-methods.R: simulate_prior.glmb, check_prior.glmb,
## describe_prior.glmb
## -----------------------------------------------------------------------
## Covers:
##   1. simulate_prior.glmb() shape/colnames match get_parameters.glmb()
##   2. simulate_prior.glmb() draws from the true joint N(mu, Sigma), i.e.
##      recovers off-diagonal correlation, not just per-parameter marginals
##   3. simulate_prior.glmb() dispersion column absent when fixed, present
##      (and inverse-gamma-distributed) when sampled
##   4. simulate_prior.glmb() exact (non-normal-surrogate) draws for dBeta
##   5. check_prior.glmb() "gelman"/"lakeland" methods run and return the
##      expected shape; a near-flat prior is "uninformative", a
##      near-point-mass prior is "informative"
##   6. describe_prior.glmb() returns pfamily(model), identical to
##      get_priors.glmb() and pfamily() itself
##   8. simulate_prior.glmb() delegates to pfamily$pfun (same draws as direct call)
##
## bayestestR is a hard dependency (Imports) and its generics are re-exported
## by glmbayes (R/reexports.R), so these are called unqualified below rather
## than as bayestestR::simulate_prior()/etc.

test_that("simulate_prior.glmb() shape/colnames match get_parameters.glmb()", {
  counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  outcome <- gl(3, 1, 9)
  treatment <- gl(3, 3)
  ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())

  fit <- glmb(
    n = 200,
    counts ~ outcome + treatment,
    family = poisson(),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    verbose = FALSE,
    use_parallel = FALSE
  )

  prior_draws <- simulate_prior(fit, n = 500)
  expect_s3_class(prior_draws, "data.frame")
  expect_equal(nrow(prior_draws), 500L)
  expect_identical(colnames(prior_draws), colnames(fit$coefficients))
  expect_false("dispersion" %in% colnames(prior_draws))
})

test_that("simulate_prior.glmb() recovers the true joint covariance, not just marginals", {
  set.seed(2026)
  n <- 40L
  x1 <- rnorm(n)
  x2 <- 0.5 * x1 + rnorm(n, sd = 0.5)
  y <- as.numeric(cbind(1, x1, x2) %*% c(0.5, -1.2, 0.8) + rnorm(n, sd = 1.5))
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  mu <- c(0, 0, 0)
  Sigma <- matrix(c(10, 0, 0,
                     0, 10, 4,
                     0,  4, 10), nrow = 3)

  fit <- lmb(
    y ~ x1 + x2,
    dNormal(mu = mu, Sigma = Sigma, dispersion = 1),
    data = dat, n = 200L, verbose = FALSE, use_parallel = FALSE
  )

  set.seed(1)
  prior_draws <- simulate_prior(fit, n = 20000)

  ## Marginal variances match diag(Sigma) ...
  expect_equal(stats::var(prior_draws[["(Intercept)"]]), 10, tolerance = 0.15)
  expect_equal(stats::var(prior_draws$x1), 10, tolerance = 0.15)
  expect_equal(stats::var(prior_draws$x2), 10, tolerance = 0.15)

  ## ... but crucially the off-diagonal correlation between x1 and x2 (which a
  ## marginal-only reconstruction from insight's usual flattened Location/
  ## Scale table could never recover) matches
  ## Sigma[x1, x2] / sqrt(Sigma[x1,x1] * Sigma[x2,x2]).
  true_cor <- Sigma[2, 3] / sqrt(Sigma[2, 2] * Sigma[3, 3])
  expect_equal(stats::cor(prior_draws$x1, prior_draws$x2), true_cor, tolerance = 0.05)
  ## (Intercept) is uncorrelated with x1/x2 in Sigma above.
  expect_equal(stats::cor(prior_draws[["(Intercept)"]], prior_draws$x1), 0, tolerance = 0.05)
})

test_that("simulate_prior.glmb() dispersion column matches get_parameters.glmb()'s convention", {
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  ps <- Prior_Setup(weight ~ group, gaussian())

  ## Fixed dispersion: no dispersion column.
  fit_fixed <- lmb(
    weight ~ group,
    dNormal(mu = ps$mu, Sigma = ps$Sigma, dispersion = 1),
    n = 200L, verbose = FALSE, use_parallel = FALSE
  )
  expect_false("dispersion" %in% colnames(simulate_prior(fit_fixed, n = 100)))

  ## Sampled dispersion: dispersion ~ InverseGamma(shape, rate), matching
  ## pfamily(fit_sampled)$prior_list$shape/rate.
  fit_sampled <- lmb(
    weight ~ group,
    dNormal_Gamma(ps$mu, Sigma_0 = ps$Sigma_0, shape = ps$shape, rate = ps$rate),
    n = 200L, verbose = FALSE, use_parallel = FALSE
  )
  set.seed(1)
  prior_draws <- simulate_prior(fit_sampled, n = 20000)
  expect_true("dispersion" %in% colnames(prior_draws))
  expect_true(all(prior_draws$dispersion > 0))

  ## InverseGamma(shape, rate) has mean rate / (shape - 1) for shape > 1.
  if (ps$shape > 1) {
    expect_equal(mean(prior_draws$dispersion), ps$rate / (ps$shape - 1), tolerance = 0.1)
  }
})

test_that("simulate_prior.glmb() draws exactly from Beta(shape1, shape2) for dBeta", {
  b_init <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list(NULL, "(Intercept)"))
  pf <- dBeta(shape1 = 2, shape2 = 7, beta = b_init)

  set.seed(101)
  y_dat <- rbinom(30, size = 1, prob = 0.3)
  fit <- glmb(
    n = 200,
    y_dat ~ 1,
    data = data.frame(y_dat = y_dat),
    weights = rep(1L, 30),
    family = binomial(link = "identity"),
    pfamily = pf,
    verbose = FALSE,
    use_parallel = FALSE
  )

  set.seed(1)
  prior_draws <- simulate_prior(fit, n = 20000)
  expect_true(all(prior_draws[["(Intercept)"]] > 0 & prior_draws[["(Intercept)"]] < 1))
  ## Beta(2, 7): mean = 2/9, variance = 2*7 / (9^2 * 10).
  expect_equal(mean(prior_draws[["(Intercept)"]]), 2 / 9, tolerance = 0.02)
  expect_equal(stats::var(prior_draws[["(Intercept)"]]), 2 * 7 / (81 * 10), tolerance = 0.01)
})

test_that("simulate_prior.glmb() truncates the dIndependent_Normal_Gamma() dispersion prior", {
  ## dIndependent_Normal_Gamma()'s dispersion prior is *always* two-sided
  ## truncated in practice: rglmb()/rlmb()/glmb()/lmb() write the envelope's
  ## actual disp_lower/disp_upper bounds back into the fitted pfamily even
  ## when the user leaves them at their default NULL (see
  ## rindepNormalGamma_reg() / R/simfunction.R). simulate_prior() must sample
  ## that truncated Inverse-Gamma, not the plain (wider) untruncated one.
  data(mtcars)
  mt <- mtcars
  mt$c_wt  <- as.numeric(scale(mtcars$wt, center = TRUE, scale = FALSE))
  mt$c_cyl <- as.numeric(scale(mtcars$cyl, center = TRUE, scale = FALSE))
  ps <- Prior_Setup(mpg ~ c_wt + c_cyl, gaussian(), data = mt, n_prior = 5)

  fit <- lmb(
    mpg ~ c_wt + c_cyl,
    dIndependent_Normal_Gamma(mu = ps$mu, Sigma = ps$Sigma, shape = ps$shape_ING, rate = ps$rate),
    data = mt, n = 200L, verbose = FALSE, use_parallel = FALSE
  )

  disp_lower <- fit$pfamily$prior_list$disp_lower
  disp_upper <- fit$pfamily$prior_list$disp_upper
  expect_true(is.numeric(disp_lower) && is.numeric(disp_upper))

  set.seed(1)
  prior_draws <- simulate_prior(fit, n = 20000)
  expect_true("dispersion" %in% colnames(prior_draws))
  expect_true(all(prior_draws$dispersion >= disp_lower - 1e-8))
  expect_true(all(prior_draws$dispersion <= disp_upper + 1e-8))

  ## Confirm the bounds are actually binding, so this test would fail on an
  ## un-truncated Inverse-Gamma(shape, rate) draw.
  shape <- fit$pfamily$prior_list$shape
  rate  <- fit$pfamily$prior_list$rate
  set.seed(1)
  untrunc <- 1 / rgamma(20000, shape = shape, rate = rate)
  expect_true(mean(untrunc < disp_lower | untrunc > disp_upper) > 0.05)
})

test_that("check_prior.glmb() distinguishes an uninformative prior from an informative one", {
  set.seed(2026)
  n <- 60L
  x1 <- rnorm(n)
  y <- as.numeric(2 + 3 * x1 + rnorm(n, sd = 0.5))
  dat <- data.frame(y = y, x1 = x1)

  ## Vague/wide prior relative to the data -> "uninformative".
  fit_vague <- lmb(
    y ~ x1,
    dNormal(mu = c(0, 0), Sigma = diag(c(1000, 1000)), dispersion = 1),
    data = dat, n = 500L, verbose = FALSE, use_parallel = FALSE
  )
  cp_vague <- check_prior(fit_vague)
  expect_s3_class(cp_vague, "data.frame")
  expect_identical(colnames(cp_vague), c("Parameter", "Prior_Quality"))
  expect_identical(cp_vague$Prior_Quality[cp_vague$Parameter == "x1"], "uninformative")

  ## Tight/informative prior relative to the data -> "informative".
  fit_tight <- lmb(
    y ~ x1,
    dNormal(mu = c(0, 0), Sigma = diag(c(1e-6, 1e-6)), dispersion = 1),
    data = dat, n = 500L, verbose = FALSE, use_parallel = FALSE
  )
  cp_tight <- check_prior(fit_tight)
  expect_identical(cp_tight$Prior_Quality[cp_tight$Parameter == "x1"], "informative")

  ## "lakeland" method runs and returns the same shape.
  cp_lakeland <- check_prior(fit_vague, method = "lakeland")
  expect_identical(colnames(cp_lakeland), c("Parameter", "Prior_Quality"))
  expect_true(all(cp_lakeland$Prior_Quality %in% c("informative", "misinformative", "not determinable")))
})

test_that("describe_prior.glmb() returns pfamily(model), not a flattened table", {
  ## Correlated multivariate normal prior: describe_prior() must preserve the
  ## full off-diagonal Sigma, unlike insight's usual marginal-only
  ## Location/Scale table.
  set.seed(2026)
  n <- 40L
  x1 <- rnorm(n)
  x2 <- 0.5 * x1 + rnorm(n, sd = 0.5)
  y <- as.numeric(cbind(1, x1, x2) %*% c(0.5, -1.2, 0.8) + rnorm(n, sd = 1.5))
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  mu <- c(0, 0, 0)
  Sigma <- matrix(c(10, 0, 0,
                     0, 10, 4,
                     0,  4, 10), nrow = 3)

  fit <- lmb(
    y ~ x1 + x2,
    dNormal(mu = mu, Sigma = Sigma, dispersion = 1),
    data = dat, n = 200L, verbose = FALSE, use_parallel = FALSE
  )

  dp <- describe_prior(fit)
  expect_identical(dp, pfamily(fit))
  expect_s3_class(dp, "pfamily")
  ## The full covariance matrix, including the off-diagonal correlation
  ## between x1 and x2, survives untouched.
  expect_equal(dp$prior_list$Sigma, Sigma, ignore_attr = TRUE)
  ## get_priors.glmb() (R/insight-methods.R) is the same pfamily(model) call
  ## under insight's naming, so the two generics agree exactly.
  expect_identical(dp, get_priors(fit))
})

test_that("simulate_prior.glmb()/check_prior.glmb() report only dispersion for dGamma(Inv_Dispersion = TRUE)", {
  ## dGamma(Inv_Dispersion = TRUE) ("rGamma_reg" fits): the prior is on the
  ## inverse dispersion (precision) *only* -- beta is a fixed, known input,
  ## not itself a modeled/estimated quantity. simulate_prior()/check_prior()
  ## must report only dispersion and must not fabricate a "prior" for the
  ## fixed coefficients.
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  ps <- Prior_Setup(weight ~ group, family = gaussian())
  rate_dg <- if (!is.null(ps$rate_gamma)) ps$rate_gamma else ps$rate

  fit <- glmb(
    weight ~ group, family = gaussian(),
    pfamily = dGamma(shape = ps$shape, rate = rate_dg, beta = ps$coefficients),
    n = 200, verbose = FALSE, use_parallel = FALSE
  )
  expect_s3_class(fit, "rGamma_reg")

  set.seed(1)
  prior_draws <- simulate_prior(fit, n = 5000)
  expect_identical(colnames(prior_draws), "dispersion")
  expect_equal(mean(1 / prior_draws$dispersion), ps$shape / rate_dg, tolerance = 0.05)

  cp <- check_prior(fit)
  expect_identical(cp$Parameter, "dispersion")
  expect_true(cp$Prior_Quality %in% c("informative", "uninformative"))
})

test_that("simulate_prior.glmb() delegates to pfamily$pfun", {
  set.seed(2026)
  n <- 40L
  x1 <- rnorm(n)
  x2 <- 0.5 * x1 + rnorm(n, sd = 0.5)
  y <- as.numeric(cbind(1, x1, x2) %*% c(0.5, -1.2, 0.8) + rnorm(n, sd = 1.5))
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  mu <- c(0, 0, 0)
  Sigma <- matrix(c(10, 0, 0,
                     0, 10, 4,
                     0,  4, 10), nrow = 3)

  fit <- lmb(
    y ~ x1 + x2,
    dNormal(mu = mu, Sigma = Sigma, dispersion = 1),
    data = dat, n = 200L, verbose = FALSE, use_parallel = FALSE
  )

  pf <- pfamily(fit)
  expect_identical(pf$pfun, rNormal_prior)

  set.seed(1)
  via_method <- simulate_prior(fit, n = 100)
  set.seed(1)
  via_pfun <- pf$pfun(n = 100, prior_list = pf$prior_list, params = colnames(fit$coefficients))
  expect_identical(via_method, via_pfun)
})
