## insight accessor methods for glmb/lmb fits, covering all 16 of insight's
## documented "core" functions (see ?glmbayes_insight_methods and insight's
## own JOSS paper for this list): model_info, get_parameters, find_parameters,
## get_priors, find_algorithm, get_data, get_response, get_predictors,
## get_random, get_variance, find_formula, find_variables, find_terms,
## find_predictors, find_random, find_response.
##
## 6 of these (model_info, get_parameters, find_parameters, find_algorithm,
## get_data, get_priors) have custom glmb methods (R/insight-methods.R) and
## are re-exported by glmbayes, so they can be called unqualified below. The
## other 10 already work correctly via glmb/lmb's inherited "glm"/"lm" class
## entries or insight's own defaults, are NOT re-exported by glmbayes, and
## are called below with an explicit insight:: prefix.

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)

ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())
glmb.D93 <- glmb(
  n = 1000,
  counts ~ outcome + treatment,
  family = poisson(),
  pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma)
)

## ----custom glmb methods: model_info, get_parameters--------------------------
model_info(glmb.D93)$is_bayesian  ## TRUE
model_info(glmb.D93)$is_poisson   ## TRUE

draws <- get_parameters(glmb.D93)
dim(draws)   ## one row per posterior draw, one column per coefficient
head(draws)

## ----custom glmb methods: find_parameters, find_algorithm---------------------
find_parameters(glmb.D93)   ## list(conditional = <coefficient names>)
find_algorithm(glmb.D93)    ## "iid rejection sampling", not "ML"

## ----custom glmb methods: get_data, get_priors---------------------------------
get_data(glmb.D93)          ## model data frame, no environment-recovery warning

## get_priors() returns pfamily(glmb.D93) directly -- the same Call / Prior
## Family / Prior List report pfamily() itself prints, including the
## complete covariance matrix for this multivariate normal coefficient
## prior (not just its diagonal, which insight's usual per-parameter
## Location/Scale table for other model classes could not represent).
get_priors(glmb.D93)
identical(get_priors(glmb.D93), pfamily(glmb.D93))  ## TRUE

## ----already work via inherited glm methods: response/predictor discovery-----
insight::get_response(glmb.D93)     ## the observed counts
insight::get_predictors(glmb.D93)   ## data frame of predictor columns (no response)
insight::find_response(glmb.D93)    ## "counts"
insight::find_predictors(glmb.D93)  ## list(conditional = c("outcome", "treatment"))
insight::find_formula(glmb.D93)     ## list(conditional = counts ~ outcome + treatment)
insight::find_variables(glmb.D93)   ## response + conditional predictor names
insight::find_terms(glmb.D93)       ## response/conditional terms, incl. "1" for intercept

## ----already correctly NULL/empty: no random-effects structure----------------
## glmbayes fits have no random effects, so these correctly report "none".
## get_variance()'s default fallback also warns "not supported" (accurate,
## since it targets variance-decomposition models); verbose = FALSE silences
## that expected warning for this clean demo run.
insight::get_random(glmb.D93)     ## NULL
insight::get_variance(glmb.D93, verbose = FALSE)   ## NULL (no variance components)
insight::find_random(glmb.D93)    ## NULL
