## During CRAN checks, run examples sequentially.
use_parallel <- identical(Sys.getenv("NOT_CRAN"), "true")

## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
print(d.AD <- data.frame(treatment, outcome, counts))

## Set up Prior for Poisson Model
ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())
ps

## Normal prior for glmb
glmb.D93 <- glmb(
  counts ~ outcome + treatment,
  family = poisson(),
  pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
  use_parallel = use_parallel
)

## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)

## Set up prior for gaussian model
ps2 <- Prior_Setup(weight ~ group, family = gaussian())
ps2

## Conjugate Normal Prior (fixed dispersion)
lmb.D9 <- lmb(
  weight ~ group,
  pfamily = dNormal(mu = ps2$mu, ps2$Sigma, dispersion = ps2$dispersion),
  use_parallel = use_parallel
)

## Conjugate Normal_Gamma Prior
lmb.D9_v2 <- lmb(
  weight ~ group,
  pfamily = dNormal_Gamma(
    ps2$mu,
    Sigma_0 = ps2$Sigma_0,
    shape = ps2$shape,
    rate = ps2$rate
  ),
  use_parallel = use_parallel
)

## Independent_Normal_Gamma_Prior
lmb.D9_v3 <- lmb(
  weight ~ group,
  dIndependent_Normal_Gamma(
    ps2$mu,
    ps2$Sigma,
    shape = ps2$shape_ING,
    rate = ps2$rate
  ),
  use_parallel = use_parallel
)

## -------------------------------------------------------------------------
## Matrix-input bridge: use Prior_Setup outputs with rglmb() and rlmb()
## -------------------------------------------------------------------------
y <- ps2$y
x <- as.matrix(ps2$x)
wt <- rep(1, length(y))

rglmb.D9 <- rglmb(
  n = 1000,
  y = y,
  x = x,
  pfamily = dIndependent_Normal_Gamma(
    ps2$mu,
    ps2$Sigma,
    shape = ps2$shape_ING,
    rate = ps2$rate
  ),
  weights = wt,
  family = gaussian(),
  use_parallel = use_parallel
)

rlmb.D9 <- rlmb(
  n = 1000,
  y = y,
  x = x,
  pfamily = dIndependent_Normal_Gamma(
    ps2$mu,
    ps2$Sigma,
    shape = ps2$shape_ING,
    rate = ps2$rate
  ),
  weights = wt,
  use_parallel = use_parallel
)

## -------------------------------------------------------------------------
## Prior-list templates for lower-level samplers
## -------------------------------------------------------------------------
prior_list_rNormalGamma <- list(
  mu = ps2$mu,
  Sigma = ps2$Sigma_0,
  shape = ps2$shape,
  rate = ps2$rate
)

prior_list_rindepNormalGamma <- list(
  mu = ps2$mu,
  Sigma = ps2$Sigma,
  dispersion = ps2$dispersion,
  shape = ps2$shape_ING,
  rate = ps2$rate,
  Precision = solve(ps2$Sigma),
  max_disp_perc = 0.99
)

rate_dg <- if (!is.null(ps2$rate_gamma)) ps2$rate_gamma else ps2$rate
prior_list_rGamma <- list(
  beta = ps2$coefficients,
  shape = ps2$shape,
  rate = rate_dg
)

## Note: for a full dGamma run across rGamma_reg/rglmb/rlmb/glmb/lmb, see:
## example("summary.rGamma_reg")

## -------------------------------------------------------------------------
## dGamma prior illustration: Prior_Setup(shape, rate_gamma or rate) + fixed beta
## -------------------------------------------------------------------------
## For gaussian models, Prior_Setup() provides Gamma hyperparameters for the
## precision: tau = 1/dispersion ~ Gamma(shape, rate). The lmb() / rGamma_reg()
## constructors consume these via dGamma(pfamily) and rGamma_reg(prior_list),
## respectively.

## dGamma via pfamily for lmb()
lmb.D9_dGamma <- lmb(
  weight ~ group,
  pfamily = dGamma(shape = ps2$shape, rate = rate_dg, beta = ps2$coefficients),
  use_parallel = use_parallel
)

## dGamma via prior_list for rGamma_reg()
out.rGamma_reg <- rGamma_reg(
  n = 1000,
  y = y,
  x = x,
  prior_list = prior_list_rGamma,
  offset = rep(0, length(y)),
  weights = wt,
  family = gaussian(),
  use_parallel = use_parallel
)

## -------------------------------------------------------------------------
## Poisson(link = "identity"), intercept-only: `conj_poisson` + dGamma(Inv_Dispersion=FALSE)
## -------------------------------------------------------------------------
y_p <- c(rep(1L, 3L), rep(0L, 6L))
df_p <- data.frame(y = y_p)
ps_p <- Prior_Setup(
  y ~ 1,
  family = poisson(link = "identity"),
  data = df_p,
  pwt = 0.4
)
if (!is.null(ps_p$conj_poisson)) {
  cp <- ps_p$conj_poisson
  pf_conj <- dGamma(shape = cp$shape, rate = cp$rate, beta = cp$beta, Inv_Dispersion = FALSE)
  ## glmb(n = 500, formula = y ~ 1, data = df_p,
  ##      family = poisson(link = "identity"), pfamily = pf_conj)
}
