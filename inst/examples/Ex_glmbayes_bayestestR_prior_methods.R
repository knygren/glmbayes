## bayestestR prior-checking methods for glmb/lmb fits (simulate_prior,
## check_prior, describe_prior). 'bayestestR' is a hard dependency (Imports)
## and these generics are re-exported by glmbayes, so no separate
## library(bayestestR) call is needed.

## ----setup: mtcars, a multivariate normal coefficient prior with real data----
## wt and cyl are strongly correlated in mtcars (heavier cars tend to have
## more cylinders), so Prior_Setup()'s data-driven Sigma -- based on
## (X'WX)^-1 -- has a substantial off-diagonal entry between c_wt and c_cyl:
## simulate_prior() draws from this exact joint N(mu, Sigma), not just the
## per-parameter marginals insight's usual get_priors() shape would report
## for other model classes.
##
## sd = c(10, 10, 0.5) gives the intercept and c_wt a deliberately vague
## prior (SD = 10) but c_cyl a tight, informative one (SD = 0.5) centered at
## 0 -- in tension with c_cyl's actual (negative) effect on mpg -- so
## check_prior() below reports a genuine mix of "informative"/"uninformative"
## (gelman) and "informative"/"misinformative" (lakeland) verdicts, rather
## than the same verdict for every coefficient. (A per-coefficient `pwt`
## vector, unlike `sd`, does not survive Prior_Setup()'s Gaussian
## dispersion-calibration step, which recomputes Sigma from (X'WX)^-1 times
## a single scalar dispersion; `sd` is passed straight through instead.)
data(mtcars)
mt <- mtcars
mt$c_wt  <- as.numeric(scale(mtcars$wt, center = TRUE, scale = FALSE))
mt$c_cyl <- as.numeric(scale(mtcars$cyl, center = TRUE, scale = FALSE))
form <- mpg ~ c_wt + c_cyl

ps <- Prior_Setup(form, gaussian(), data = mt, sd = c(10, 10, 0.5), n_prior = 3)

fit <- lmb(
  form,
  dNormal(mu = ps$mu, Sigma = ps$Sigma, dispersion = ps$dispersion),
  data = mt, n = 2000L, verbose = FALSE, use_parallel = FALSE
)

## ----simulate_prior-----------------------------------------------------------
prior_draws <- simulate_prior(fit, n = 5000)
dim(prior_draws)                          ## 5000 draws x 3 coefficients
cor(prior_draws$c_wt, prior_draws$c_cyl)  ## recovers Sigma's off-diagonal correlation

## ----check_prior---------------------------------------------------------------
check_prior(fit)                       ## "gelman" rule (posterior SD vs. prior SD)
check_prior(fit, method = "lakeland")  ## posterior mass inside the prior's 95% HDI

## ----describe_prior-------------------------------------------------------------
## describe_prior() and get_priors() both return pfamily(fit) directly -- the
## same full report (including the complete Sigma, not just its diagonal)
## shown by pfamily(fit) itself.
describe_prior(fit)
identical(describe_prior(fit), pfamily(fit))  ## TRUE
identical(get_priors(fit), pfamily(fit))      ## TRUE

## ----dIndependent_Normal_Gamma: priors on *both* coefficients and dispersion---
## dNormal_Gamma() and dIndependent_Normal_Gamma() both place a prior on the
## coefficients (joint normal) *and* a separate inverse-gamma prior on the
## dispersion, so simulate_prior() below returns both a coefficient block and
## a dispersion column. Unlike dNormal_Gamma() (fully conjugate, no
## truncation), dIndependent_Normal_Gamma()'s independent coefficient/
## precision priors are not jointly conjugate, so its accept-reject sampler
## needs a *bounded* dispersion domain -- the dispersion prior is therefore
## two-sided truncated to [disp_lower, disp_upper] (derived from
## max_disp_perc; here left at its default). simulate_prior() draws from this
## exact *truncated* Inverse-Gamma, not the wider untruncated one.
ps_ing <- Prior_Setup(form, gaussian(), data = mt, n_prior = 5)

fit_ing <- lmb(
  form,
  dIndependent_Normal_Gamma(mu = ps_ing$mu, Sigma = ps_ing$Sigma,
                             shape = ps_ing$shape_ING, rate = ps_ing$rate),
  data = mt, n = 2000L, verbose = FALSE, use_parallel = FALSE
)

## The fit writes the envelope's actual bounds back into pfamily(fit_ing),
## even though disp_lower/disp_upper were left at their default NULL above.
fit_ing$pfamily$prior_list[c("disp_lower", "disp_upper")]

prior_draws_ing <- simulate_prior(fit_ing, n = 5000)
colnames(prior_draws_ing)                 ## coefficients *and* dispersion
range(prior_draws_ing$dispersion)         ## stays within [disp_lower, disp_upper]

check_prior(fit_ing)  ## "dispersion" is checked alongside the coefficients
