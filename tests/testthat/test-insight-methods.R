## Tests for insight-methods.R: model_info.glmb, get_parameters.glmb,
## find_parameters.glmb, find_algorithm.glmb, get_data.glmb, get_priors.glmb
## -----------------------------------------------------------------------
## Covers:
##   1. Smoke test: model_info(fit)$is_bayesian is TRUE (glmb + lmb)
##   2. model_info.glmb family/link-specific flags for Poisson
##   3. get_parameters.glmb shape/colnames match fit$coefficients
##   4. dispersion column absent when dispersion is fixed, present when sampled
##   5. find_parameters.glmb component split (conditional/dispersion)
##   6. find_algorithm.glmb reports iid rejection sampling, not "ML"/"OLS"
##   7. get_data.glmb returns the model data with no warning, with or without
##      an explicit data= argument to glmb()/lmb()
##   8. get_priors.glmb returns pfamily(model) directly (identical(), full
##      covariance preserved), across dNormal, dNormal_Gamma, and dBeta
##
## insight is a hard dependency (Imports) and its generics are re-exported by
## glmbayes (R/reexports.R), so these are called unqualified below rather
## than as insight::model_info()/insight::get_parameters()/etc.

test_that("model_info() reports is_bayesian = TRUE for glmb and lmb fits", {
  counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  outcome <- gl(3, 1, 9)
  treatment <- gl(3, 3)

  ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())
  fit_glmb <- glmb(
    n = 200,
    counts ~ outcome + treatment,
    family = poisson(),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    verbose = FALSE,
    use_parallel = FALSE
  )

  ## Smoke test: run this first, before any other insight-related assertion.
  expect_true(model_info(fit_glmb)$is_bayesian)

  set.seed(2026)
  n <- 20L
  x1 <- rnorm(n)
  y <- as.numeric(cbind(1, x1) %*% c(0.5, -1.2) + rnorm(n, sd = 1.5))
  dat <- data.frame(y = y, x1 = x1)

  fit_lmb <- lmb(
    y ~ x1,
    dNormal(mu = c(0, 0), Sigma = diag(c(10, 10)), dispersion = 1),
    data = dat, n = 200L, verbose = FALSE, use_parallel = FALSE
  )

  expect_true(model_info(fit_lmb)$is_bayesian)
})

test_that("model_info.glmb() derives family/link-specific flags correctly", {
  counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  outcome <- gl(3, 1, 9)
  treatment <- gl(3, 3)
  ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())

  fit_glmb <- glmb(
    n = 200,
    counts ~ outcome + treatment,
    family = poisson(),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    verbose = FALSE,
    use_parallel = FALSE
  )

  mi <- model_info(fit_glmb)
  expect_true(mi$is_count)
  expect_true(mi$is_poisson)
  expect_false(mi$is_linear)
  expect_identical(mi$family, "poisson")
  expect_identical(mi$link_function, "log")
  expect_identical(mi$n_obs, length(counts))

  set.seed(2026)
  n <- 20L
  x1 <- rnorm(n)
  y <- as.numeric(cbind(1, x1) %*% c(0.5, -1.2) + rnorm(n, sd = 1.5))
  dat <- data.frame(y = y, x1 = x1)
  fit_lmb <- lmb(
    y ~ x1,
    dNormal(mu = c(0, 0), Sigma = diag(c(10, 10)), dispersion = 1),
    data = dat, n = 200L, verbose = FALSE, use_parallel = FALSE
  )

  ## lmb objects have no $family; model_info.glmb() must default to gaussian().
  mi_lmb <- model_info(fit_lmb)
  expect_true(mi_lmb$is_linear)
  expect_identical(mi_lmb$family, "gaussian")
  expect_identical(mi_lmb$link_function, "identity")
})

test_that("get_parameters.glmb() shape/colnames match fit$coefficients", {
  counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  outcome <- gl(3, 1, 9)
  treatment <- gl(3, 3)
  ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())

  fit_glmb <- glmb(
    n = 200,
    counts ~ outcome + treatment,
    family = poisson(),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    verbose = FALSE,
    use_parallel = FALSE
  )

  draws <- get_parameters(fit_glmb)
  expect_s3_class(draws, "data.frame")
  expect_equal(dim(draws), dim(fit_glmb$coefficients))
  expect_identical(colnames(draws), colnames(fit_glmb$coefficients))
  expect_equal(as.matrix(draws), fit_glmb$coefficients, ignore_attr = TRUE)

  ## Fixed dispersion (dNormal, dispersion = 1): no dispersion column.
  expect_false("dispersion" %in% colnames(draws))
})

test_that("get_parameters.glmb() appends a dispersion column only when sampled", {
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)

  ps <- Prior_Setup(weight ~ group, gaussian())

  ## Fixed dispersion: dNormal() with an explicit scalar dispersion.
  fit_fixed <- lmb(
    weight ~ group,
    dNormal(mu = ps$mu, Sigma = ps$Sigma, dispersion = 1),
    n = 200L, verbose = FALSE, use_parallel = FALSE
  )
  draws_fixed <- get_parameters(fit_fixed)
  expect_false("dispersion" %in% colnames(draws_fixed))

  ## Sampled dispersion: dNormal_Gamma() samples dispersion jointly with coefficients.
  fit_sampled <- lmb(
    weight ~ group,
    dNormal_Gamma(ps$mu, Sigma_0 = ps$Sigma_0, shape = ps$shape, rate = ps$rate),
    n = 200L, verbose = FALSE, use_parallel = FALSE
  )
  draws_sampled <- get_parameters(fit_sampled)
  expect_true("dispersion" %in% colnames(draws_sampled))
  expect_equal(draws_sampled$dispersion, fit_sampled$dispersion)
  expect_true(length(unique(draws_sampled$dispersion)) > 1L)
})

test_that("find_parameters.glmb() splits conditional/dispersion components", {
  counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  outcome <- gl(3, 1, 9)
  treatment <- gl(3, 3)
  ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())

  fit_glmb <- glmb(
    n = 200,
    counts ~ outcome + treatment,
    family = poisson(),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    verbose = FALSE,
    use_parallel = FALSE
  )

  fp <- find_parameters(fit_glmb)
  expect_identical(fp, list(conditional = colnames(fit_glmb$coefficients)))

  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  ps2 <- Prior_Setup(weight ~ group, gaussian())

  fit_sampled <- lmb(
    weight ~ group,
    dNormal_Gamma(ps2$mu, Sigma_0 = ps2$Sigma_0, shape = ps2$shape, rate = ps2$rate),
    n = 200L, verbose = FALSE, use_parallel = FALSE
  )
  fp_sampled <- find_parameters(fit_sampled)
  expect_identical(
    fp_sampled,
    list(conditional = colnames(fit_sampled$coefficients), dispersion = "dispersion")
  )
})

test_that("find_algorithm.glmb() reports iid rejection sampling, not ML/OLS", {
  counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  outcome <- gl(3, 1, 9)
  treatment <- gl(3, 3)
  ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())

  fit_glmb <- glmb(
    n = 200,
    counts ~ outcome + treatment,
    family = poisson(),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    verbose = FALSE,
    use_parallel = FALSE
  )

  fa <- find_algorithm(fit_glmb)
  expect_identical(fa$algorithm, "iid rejection sampling")
  expect_identical(fa$iterations, nrow(fit_glmb$coefficients))
  ## Guard against silently falling through to insight's inherited glm/lm
  ## methods (which report "ML"/"OLS" -- wrong for a Bayesian iid sampler).
  expect_false(fa$algorithm %in% c("ML", "OLS"))
})

test_that("get_data.glmb() returns model data without warning, with or without data=", {
  counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  outcome <- gl(3, 1, 9)
  treatment <- gl(3, 3)
  ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())

  ## No explicit data= argument: fit_glmb$data is the calling environment,
  ## not a data frame (as for stats::glm()); get_data() must fall back to
  ## the model frame (fit_glmb$model) instead.
  fit_glmb <- glmb(
    n = 200,
    counts ~ outcome + treatment,
    family = poisson(),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    verbose = FALSE,
    use_parallel = FALSE
  )
  expect_false(is.data.frame(fit_glmb$data))

  d <- expect_no_warning(get_data(fit_glmb))
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), length(counts))

  ## Explicit data= argument: lmb() stores the model frame directly as $data.
  set.seed(2026)
  n <- 20L
  x1 <- rnorm(n)
  y <- as.numeric(cbind(1, x1) %*% c(0.5, -1.2) + rnorm(n, sd = 1.5))
  dat <- data.frame(y = y, x1 = x1)
  fit_lmb <- lmb(
    y ~ x1,
    dNormal(mu = c(0, 0), Sigma = diag(c(10, 10)), dispersion = 1),
    data = dat, n = 200L, verbose = FALSE, use_parallel = FALSE
  )
  expect_true(is.data.frame(fit_lmb$data))

  d_lmb <- expect_no_warning(get_data(fit_lmb))
  expect_s3_class(d_lmb, "data.frame")
  expect_equal(nrow(d_lmb), n)
})

test_that("get_priors.glmb() returns pfamily(model), not a flattened table", {
  ## get_priors() must not lose information insight's usual Parameter/
  ## Distribution/Location/Scale shape cannot represent (e.g. off-diagonal
  ## covariance for a genuinely multivariate normal coefficient prior), so
  ## it returns the full pfamily object directly instead of that shape.
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  ps <- Prior_Setup(weight ~ group, gaussian())

  fit_fixed <- lmb(
    weight ~ group,
    dNormal(mu = ps$mu, Sigma = ps$Sigma, dispersion = 1),
    n = 200L, verbose = FALSE, use_parallel = FALSE
  )
  gp_fixed <- get_priors(fit_fixed)
  expect_identical(gp_fixed, pfamily(fit_fixed))
  expect_s3_class(gp_fixed, "pfamily")
  expect_equal(gp_fixed$prior_list$mu, ps$mu, ignore_attr = TRUE)
  expect_equal(gp_fixed$prior_list$Sigma, ps$Sigma, ignore_attr = TRUE)

  fit_sampled <- lmb(
    weight ~ group,
    dNormal_Gamma(ps$mu, Sigma_0 = ps$Sigma_0, shape = ps$shape, rate = ps$rate),
    n = 200L, verbose = FALSE, use_parallel = FALSE
  )
  gp_sampled <- get_priors(fit_sampled)
  expect_identical(gp_sampled, pfamily(fit_sampled))
  expect_equal(gp_sampled$prior_list$shape, ps$shape)
  expect_equal(gp_sampled$prior_list$rate, ps$rate)
})

test_that("get_priors.glmb() preserves the full covariance matrix for a correlated prior", {
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

  gp <- get_priors(fit)
  expect_equal(gp$prior_list$Sigma, Sigma, ignore_attr = TRUE)
  ## The off-diagonal correlation (which a flattened Location/Scale table
  ## could never represent) survives untouched.
  expect_equal(gp$prior_list$Sigma[2, 3], 4)
})

test_that("get_priors.glmb() returns pfamily(model) for a dBeta prior too", {
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

  gp <- get_priors(fit)
  expect_identical(gp, pfamily(fit))
  expect_identical(gp$prior_list$shape1, 2)
  expect_identical(gp$prior_list$shape2, 7)
})

test_that("get_parameters.glmb()/find_parameters.glmb() omit coefficients for dGamma(Inv_Dispersion = TRUE)", {
  ## dGamma(Inv_Dispersion = TRUE) ("rGamma_reg" fits): the prior is on the
  ## inverse dispersion (precision) *only* -- beta is a fixed, known input
  ## (typically held fixed in a Gibbs step), not itself a modeled/estimated
  ## quantity. get_parameters()/find_parameters() must report *only*
  ## dispersion and omit the coefficients entirely, rather than treating the
  ## fixed beta as if it were a (degenerate) posterior draw (regression test
  ## for a design issue found while checking all five pfamily constructors
  ## against these methods).
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
  expect_identical(nrow(fit$coefficients), 1L)
  expect_gt(length(unique(fit$dispersion)), 1L)

  gpars <- get_parameters(fit)
  expect_identical(colnames(gpars), "dispersion")
  expect_identical(nrow(gpars), length(fit$dispersion))
  expect_identical(gpars$dispersion, fit$dispersion)

  fp <- find_parameters(fit)
  expect_identical(fp, list(dispersion = "dispersion"))
})

test_that("get_priors.glmb() returns pfamily(model) for dIndependent_Normal_Gamma too", {
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  ps <- Prior_Setup(weight ~ group, gaussian())

  fit <- lmb(
    weight ~ group,
    dIndependent_Normal_Gamma(ps$mu, ps$Sigma, shape = ps$shape_ING, rate = ps$rate),
    n = 200L, verbose = FALSE, use_parallel = FALSE
  )
  gp <- get_priors(fit)
  expect_identical(gp, pfamily(fit))
  expect_equal(gp$prior_list$Sigma, ps$Sigma, ignore_attr = TRUE)
})
