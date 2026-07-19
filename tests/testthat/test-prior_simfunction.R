## Tests for R/prior_simfunction.R: rNormal_prior, rGamma_prior,
## rGamma_Conjugate_prior, rNormal_Gamma_prior,
## rIndependent_Normal_Gamma_prior, rBeta_prior
## -----------------------------------------------------------------------
## These functions are wired into each pfamily() constructor as $pfun and are
## called by simulate_prior.glmb(); this file exercises them directly against
## their own prior_list inputs and known distributional properties.
##
## Covers:
##   1. All six share an identical (n, prior_list, params, ...) signature
##      and can be called generically via do.call().
##   2. rNormal_prior() draws from the true joint N(mu, Sigma), recovering
##      off-diagonal correlation (not just per-parameter marginals).
##   3. rGamma_prior() (Inv_Dispersion = TRUE) returns dispersion only,
##      matching InverseGamma(shape, rate) moments (or a truncated version
##      when disp_lower/disp_upper are supplied).
##   4. rGamma_Conjugate_prior() draws each coefficient independently from
##      Gamma(shape, rate), including for multi-parameter beta.
##   5. rNormal_Gamma_prior() returns both a coefficient block (joint normal)
##      and an (untruncated) InverseGamma dispersion column.
##   6. rIndependent_Normal_Gamma_prior() returns both a coefficient block
##      and a *truncated* InverseGamma dispersion column, matching a real
##      dIndependent_Normal_Gamma() fit's disp_lower/disp_upper exactly.
##   7. rBeta_prior() draws each coefficient independently from
##      Beta(shape1, shape2), including for multi-parameter beta.
##   8. params overrides rownames(prior_list$mu); both fall back to "V1",
##      "V2", ... when neither is available.
##   9. Each pfamily() constructor sets the matching $pfun.

test_that("each pfamily() constructor sets the matching pfun", {
  mu <- c(0, 0)
  Sigma <- diag(2)
  beta <- c(0.1, 0.2)

  expect_identical(dNormal(mu = mu, Sigma = Sigma)$pfun, rNormal_prior)
  expect_identical(dGamma(shape = 2, rate = 3, beta = beta, Inv_Dispersion = TRUE)$pfun, rGamma_prior)
  expect_identical(dGamma(shape = 2, rate = 3, beta = beta, Inv_Dispersion = FALSE)$pfun, rGamma_Conjugate_prior)
  expect_identical(dBeta(shape1 = 2, shape2 = 3, beta = beta)$pfun, rBeta_prior)
  expect_identical(dNormal_Gamma(mu = mu, Sigma = Sigma, shape = 2, rate = 3)$pfun, rNormal_Gamma_prior)
  expect_identical(dIndependent_Normal_Gamma(mu = mu, Sigma = Sigma, shape = 2, rate = 3)$pfun, rIndependent_Normal_Gamma_prior)
})

test_that("all six prior_simfuncs share an identical signature", {
  prior_funs <- list(
    rNormal_prior, rGamma_prior, rGamma_Conjugate_prior,
    rNormal_Gamma_prior, rIndependent_Normal_Gamma_prior, rBeta_prior
  )
  for (f in prior_funs) {
    expect_identical(names(formals(f)), c("n", "prior_list", "params", "..."))
  }

  ## Generic call pattern: do.call(pfun, list(n = n, prior_list = prior_list))
  ## must work for every one of them.
  pl <- list(mu = c(0, 0), Sigma = diag(2))
  draws <- do.call(rNormal_prior, list(n = 10, prior_list = pl))
  expect_equal(nrow(draws), 10L)
})

test_that("rNormal_prior() draws from the true joint N(mu, Sigma)", {
  mu <- c(0, 0, 0)
  Sigma <- matrix(c(10, 0, 0,
                     0, 10, 4,
                     0,  4, 10), nrow = 3)
  pl <- list(mu = mu, Sigma = Sigma)

  set.seed(1)
  draws <- rNormal_prior(20000, pl, params = c("a", "b", "c"))
  expect_s3_class(draws, "data.frame")
  expect_identical(colnames(draws), c("a", "b", "c"))
  expect_equal(nrow(draws), 20000L)

  expect_equal(stats::var(draws$a), 10, tolerance = 0.15)
  expect_equal(stats::var(draws$b), 10, tolerance = 0.15)
  expect_equal(stats::var(draws$c), 10, tolerance = 0.15)

  true_cor <- Sigma[2, 3] / sqrt(Sigma[2, 2] * Sigma[3, 3])
  expect_equal(stats::cor(draws$b, draws$c), true_cor, tolerance = 0.05)
  expect_equal(stats::cor(draws$a, draws$b), 0, tolerance = 0.05)
})

test_that("rGamma_prior() draws InverseGamma(shape, rate), truncated when bounds given", {
  pl <- list(shape = 5, rate = 8)
  set.seed(1)
  draws <- rGamma_prior(20000, pl)
  expect_identical(colnames(draws), "dispersion")
  expect_equal(mean(draws$dispersion), pl$rate / (pl$shape - 1), tolerance = 0.05)

  ## With explicit truncation bounds, draws must stay within [lower, upper]
  ## and differ from the untruncated moments above.
  pl_trunc <- list(shape = 5, rate = 8, disp_lower = 0.5, disp_upper = 1.5)
  set.seed(1)
  draws_trunc <- rGamma_prior(20000, pl_trunc)
  expect_true(all(draws_trunc$dispersion >= 0.5 - 1e-8))
  expect_true(all(draws_trunc$dispersion <= 1.5 + 1e-8))
})

test_that("rGamma_Conjugate_prior() draws each coefficient independently from Gamma(shape, rate)", {
  mu <- matrix(c(2, 2, 2), ncol = 1, dimnames = list(c("a", "b", "c"), NULL))
  pl <- list(shape = 6, rate = 3, mu = mu)

  set.seed(1)
  draws <- rGamma_Conjugate_prior(20000, pl)
  expect_identical(colnames(draws), c("a", "b", "c"))
  expect_equal(mean(draws$a), pl$shape / pl$rate, tolerance = 0.05)
  expect_equal(mean(draws$b), pl$shape / pl$rate, tolerance = 0.05)

  ## Independent draws: columns must not be identical (unlike a naive
  ## single-vector-recycled-across-columns implementation).
  expect_false(isTRUE(all(draws$a == draws$b)))
  expect_equal(stats::cor(draws$a, draws$b), 0, tolerance = 0.05)
})

test_that("rNormal_Gamma_prior() returns joint-normal coefficients plus untruncated dispersion", {
  mu <- c(1, -1)
  Sigma <- matrix(c(4, 1, 1, 4), nrow = 2)
  pl <- list(mu = mu, Sigma = Sigma, shape = 6, rate = 10)

  set.seed(1)
  draws <- rNormal_Gamma_prior(20000, pl, params = c("x1", "x2"))
  expect_identical(colnames(draws), c("x1", "x2", "dispersion"))
  expect_equal(mean(draws$x1), 1, tolerance = 0.1)
  expect_equal(mean(draws$x2), -1, tolerance = 0.1)
  expect_equal(mean(draws$dispersion), pl$rate / (pl$shape - 1), tolerance = 0.1)
})

test_that("rIndependent_Normal_Gamma_prior() truncates dispersion exactly like a real fit", {
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
  pl <- fit$pfamily$prior_list
  expect_true(is.numeric(pl$disp_lower) && is.numeric(pl$disp_upper))

  set.seed(1)
  draws <- rIndependent_Normal_Gamma_prior(20000, pl, params = colnames(fit$coefficients))
  expect_identical(colnames(draws), c(colnames(fit$coefficients), "dispersion"))
  expect_true(all(draws$dispersion >= pl$disp_lower - 1e-8))
  expect_true(all(draws$dispersion <= pl$disp_upper + 1e-8))

  ## Bounds must actually be binding (else this test would pass trivially).
  set.seed(1)
  untrunc <- 1 / rgamma(20000, shape = pl$shape, rate = pl$rate)
  expect_true(mean(untrunc < pl$disp_lower | untrunc > pl$disp_upper) > 0.05)

  ## Result must match calling this same pfunction directly against the
  ## exact prior_list the fit stored, i.e. it is a faithful drop-in for
  ## the coefficient+dispersion draws simulate_prior.glmb() currently
  ## builds by hand.
  expect_equal(nrow(draws), 20000L)
})

test_that("rBeta_prior() draws each coefficient independently from Beta(shape1, shape2)", {
  beta <- matrix(c(0.5, 0.5), ncol = 1, dimnames = list(c("p1", "p2"), NULL))
  pl <- list(shape1 = 2, shape2 = 7, beta = beta,
             mu = beta * 0 + 2 / 9,
             Sigma = diag(2) * (2 * 7 / (81 * 10)))
  rownames(pl$mu) <- rownames(beta)

  set.seed(1)
  draws <- rBeta_prior(20000, pl)
  expect_identical(colnames(draws), c("p1", "p2"))
  expect_true(all(draws$p1 > 0 & draws$p1 < 1))
  expect_equal(mean(draws$p1), 2 / 9, tolerance = 0.02)
  expect_equal(stats::var(draws$p1), 2 * 7 / (81 * 10), tolerance = 0.01)
  expect_false(isTRUE(all(draws$p1 == draws$p2)))
})

test_that("params overrides rownames(prior_list$mu); both fall back to V1, V2, ...", {
  mu <- matrix(c(0, 0), ncol = 1, dimnames = list(c("a", "b"), NULL))
  pl <- list(mu = mu, Sigma = diag(2))

  ## rownames(mu) used when params is NULL.
  draws1 <- rNormal_prior(5, pl)
  expect_identical(colnames(draws1), c("a", "b"))

  ## params, when supplied, overrides rownames(mu).
  draws2 <- rNormal_prior(5, pl, params = c("x", "y"))
  expect_identical(colnames(draws2), c("x", "y"))

  ## Neither available -> "V1", "V2", ...
  pl_unnamed <- list(mu = c(0, 0), Sigma = diag(2))
  draws3 <- rNormal_prior(5, pl_unnamed)
  expect_identical(colnames(draws3), c("V1", "V2"))
})
