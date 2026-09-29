# ---- Prior_Setup roxygen (inactive); live help from @inherit glmbayesCore::Prior_Setup in R/reexports.R ----
# #' Setup Prior Objects
# #'
# #' Helper function to facilitate the Setup of Prior Distributions for glm models.
# #' @name Prior_Setup
# #' @param na.action how \code{NAs} are treated. The default is first, any \code{\link{na.action}} attribute of
# #' data, second a \code{na.action} setting of \link{options}, and third \code{na.fail} if that is unset.
# #' The \code{factory-fresh} default is \code{na.omit}. Another possible value is \code{NULL}.
# #' @param family a description of the error distribution and link function to be used in the model.
# #' @param pwt Weight on the prior relative to the likelihood function at the maximum likelihood
# #' estimate. If supplied, this value is used directly (scalar or one value per coefficient).
# #' If \code{n_prior} is provided and \code{pwt} is still a **scalar** and \code{sd} was **not**
# #' supplied, \code{pwt} is set to \code{n_prior / (n_prior + n_effective)}. If \code{length(pwt) > 1}
# #' (including from \code{sd}) or \code{sd} was supplied, \code{n_prior} does **not** overwrite
# #' \code{pwt}; it is used only as a scalar for Gamma / S_marg steps. If \code{sd} is provided,
# #' \code{pwt} is computed from the prior standard deviations. If none of these are supplied,
# #' \code{pwt} defaults to \code{pwt_default_low} for models with fewer than 14 coefficients, and
# #' \code{pwt_default_high} otherwise.
# #' @param pwt_default_low Default prior weight used when \code{pwt} is not supplied and the model
# #' dimension is below 14. Defaults to 0.01.
# #' @param pwt_default_high Default prior weight used when \code{pwt} is not supplied and the model
# #' dimension is 14 or greater. Defaults to 0.05.
# #' @param n_prior Optional scalar effective prior sample size (on the \code{n_effective} scale).
# #' If provided with scalar \code{pwt} and without \code{sd}, \code{pwt} is recomputed from
# #' \code{n_prior}. With vector \code{pwt} or with \code{sd}, \code{pwt} is left unchanged and
# #' \code{n_prior} is used for the Gamma prior on precision and related Gaussian calibration only.
# #' If missing and \code{pwt} is scalar, \code{n_prior = (pwt/(1-pwt))*n_effective}.
# #' @param sd Optional vector argument with the prior standard deviations for the coefficients
# #' @param dispersion Optional scalar dispersion override (default \code{NULL}).
# #' For now, this is documented as an optional argument used to scale the
# #' \code{Sigma} (variance-covariance) matrix; see Details for additional context.
# #' @param k Scalar (default \code{1}), non-negative (\eqn{k \geq 0}), with \eqn{k + p \geq 2}
# #'   where \eqn{p} is the number of coefficients (columns of the model matrix). \code{k}
# #'   controls the tail behavior and effective degrees of freedom of the variance prior. It does
# #'   not change the posterior mean of \eqn{\sigma^2} or the covariance of \eqn{\beta}, but larger
# #'   \code{k} makes the prior and posterior for \eqn{\sigma^2} more concentrated and less
# #'   heavy-tailed. Not yet used in calibration; passed through to \code{\link{compute_gaussian_prior}}
# #'   for future use.
# #' @param intercept_source Specifies the method through which the prior mean for the intercept term is set. Options are based on the null intercept only model (null_model) or full_models. The default is the null model which is safer if variables are not centered.
# #' @param effects_source Specifies the method through which the prior means for the effects terms are set. Options are null_effects (prior means set to zero) or full_model (effect means set to match maximum likelihood estimates).
# #' @param mu Optional vector argument with the prior means for the coefficients
# #' @param x An object of class \code{"PriorSetup"}
# #' @inheritParams stats::glm
# #' @inheritParams stats::model.frame
# #' 
# #' @details
# #' \strong{Inputs to the function}
# #'
# #' The inputs to `Prior_Setup()` fall into three conceptual categories:
# #'
# #' \emph{1. Model specification}
# #' * \code{formula}: structure of the GLM (response and predictors).
# #' * \code{family}: error distribution and link.
# #' * \code{data}, \code{weights}, \code{subset}, \code{na.action},
# #'   \code{offset}, \code{contrasts}, \code{control}, \code{...}: as in
# #'   \code{\link[stats]{glm}}.
# #'
# #' \emph{2. Prior variance–covariance specification}
# #' * \code{pwt}: prior weight relative to the likelihood. If scalar, used to
# #'   construct a Zellner-type g-prior. If vector, applied elementwise.
# #' * \code{n_prior}: optional scalar effective prior sample size. Replaces
# #'   scalar \code{pwt} only when \code{pwt} is scalar and \code{sd} is not
# #'   used; otherwise supplies precision-prior / calibration only.
# #' * \code{sd}: optional vector of prior standard deviations. If provided,
# #'   used to compute \code{pwt} from the diagonal of \code{vcov(glm_full)}.
# #' * \code{pwt_default_low}, \code{pwt_default_high}: defaults for \code{pwt}
# #'   when not supplied.
# #'
# #' \emph{3. Prior mean specification}
# #' * \code{intercept_source}: method for setting the prior mean of the intercept
# #'   (\code{"null_model"} or \code{"full_model"}).
# #' * \code{effects_source}: method for setting the prior mean of the effects
# #'   (\code{"null_effects"} or \code{"full_model"}).
# #' * \code{mu}: optional user-specified prior mean vector; overrides other
# #'   centering logic if provided.
# #'
# #' \strong{Prior covariance and Zellner scaling}
# #'
# #' Let \eqn{V_0 = \mathrm{vcov}(\hat\beta)} be the covariance matrix of the
# #' full-model GLM coefficients. For non-Gaussian families, the prior covariance
# #' is:
# #' \deqn{
# #'   \Sigma =
# #'   \begin{cases}
# #'     \dfrac{1 - \mathrm{pwt}}{\mathrm{pwt}} V_0, & \text{scalar pwt},\\[4pt]
# #'     V_0 \circ \left[\sqrt{\dfrac{1 - \mathrm{pwt}_i}{\mathrm{pwt}_i}}
# #'                    \sqrt{\dfrac{1 - \mathrm{pwt}_j}{\mathrm{pwt}_j}}\right],
# #'     & \text{vector pwt},
# #'   \end{cases}
# #' }
# #' where \eqn{\circ} denotes elementwise multiplication.
# #'
# #' \strong{Intercept-only Poisson(link = "identity") conjugate prior on the rate}
# #'
# #' When the design is a single column (intercept only), \code{family = poisson()},
# #' \code{link = "identity"}, scalar \code{pwt}, and offsets are zero, the effective
# #' conjugate prior observation count \eqn{n_{\mathrm{prior}}} already satisfies
# #' \eqn{n_{\mathrm{prior}} / (n_{\mathrm{prior}} + n_{\mathrm{effective}}) = \mathrm{pwt}}
# #' (otherwise \code{conj_poisson} remains \code{NULL} with a warning).  Writing the
# #' weighted mean \eqn{\bar{y}_w = \sum_i w_i y_i / \sum_i w_i}, the output list
# #' component \code{conj_poisson} stores
# #' \eqn{\texttt{shape} = n_{\mathrm{prior}} \bar{y}_w} and \eqn{\texttt{rate} =
# #' n_{\mathrm{prior}}}, so the prior mean for the rate matches \eqn{\bar{y}_w}.
# #' Omitting optional \code{mu} resets the surrogate Normal summaries \code{mu} and the
# #' sole diagonal element of \code{Sigma} to these Gamma moments (\code{Sigma_{11}} =
# #' \eqn{\bar{y}_w / n_{\mathrm{prior}}}).
# #'
# #' For Gaussian families, `Prior_Setup()` also constructs the dispersion-free
# #' covariance
# #' \deqn{
# #'   \Sigma_0 = \Sigma / \texttt{dispersion},
# #' }
# #' which under scalar \code{pwt} and the default calibration reduces to
# #' \deqn{
# #'   \Sigma_0 = \frac{1 - \mathrm{pwt}}{\mathrm{pwt}} (X^\top W X)^{-1}.
# #' }
# #'
# #' \strong{Gaussian Normal–Gamma calibration and \eqn{S_{\mathrm{marg}}}}
# #'
# #' For \code{family = gaussian()}, the function performs the Normal–Gamma
# #' calibration described in \insertCite{glmbayesChapterA12}{glmbayes}. Let:
# #' * \eqn{p = \texttt{ncol}(x)},
# #' * \eqn{n_{\mathrm{effective}} = \sum_i w_i},
# #' * \eqn{\hat\beta} the weighted least-squares estimator,
# #' * \eqn{\Sigma_0} the dispersion-free prior covariance.
# #'
# #' The marginal quadratic term is
# #' \deqn{
# #'   S_{\mathrm{marg}}
# #'     = \mathrm{RSS}_w
# #'       + (\hat\beta - \mu)^\top
# #'         \left(\Sigma_0 + (X^\top W X)^{-1}\right)^{-1}
# #'         (\hat\beta - \mu),
# #' }
# #' where \eqn{\mathrm{RSS}_w} is the weighted residual sum of squares at
# #' \eqn{\hat\beta}. Under the default scalar-\code{pwt} Zellner mapping
# #' \eqn{\Sigma_0 = \frac{1 - \mathrm{pwt}}{\mathrm{pwt}}(X^\top W X)^{-1}},
# #' this simplifies to
# #' \deqn{
# #'   S_{\mathrm{marg}}
# #'     = \mathrm{RSS}_w
# #'       + \mathrm{pwt}\,
# #'         (\hat\beta - \mu)^\top (X^\top W X)(\hat\beta - \mu),
# #' }
# #' which makes the limiting behavior as \eqn{\mathrm{pwt} \to 0} transparent.
# #'
# #' The calibrated dispersion is
# #' \deqn{
# #'   \texttt{dispersion}
# #'     = \frac{S_{\mathrm{marg}}}{n_{\mathrm{effective}} - p},
# #' }
# #' and the Normal–Gamma hyperparameters are
# #' \deqn{
# #'   \text{shape} = \frac{n_{\mathrm{prior}} + k}{2},\qquad
# #'   \text{rate}
# #'     = \frac{1}{2} S_{\mathrm{marg}}
# #'       \frac{n_{\mathrm{prior}} + k + p - 2}{n_{\mathrm{effective}} - p}.
# #' }
# #' The independent Normal–Gamma shape is
# #' \deqn{
# #'   \text{shape}_{ING} = \text{shape} + \frac{p}{2}.
# #' }
# #'
# #' \strong{Posterior summaries for the conjugate Normal–Gamma prior}
# #'
# #' Under the conjugate Normal–Gamma prior (used by \code{dNormal_Gamma()}),
# #' the posterior has:
# #' * Posterior mean
# #'   \deqn{
# #'     E[\beta \mid y]
# #'       = (1 - \mathrm{pwt})\,\hat\beta
# #'         + \mathrm{pwt}\,\mu.
# #'   }
# #' * Posterior expectation of \eqn{\sigma^2}
# #'   \deqn{
# #'     E[\sigma^2 \mid y]
# #'       = \frac{S_{\mathrm{marg}}}{n_{\mathrm{effective}} - p}.
# #'   }
# #' * Posterior covariance
# #'   \deqn{
# #'     \mathrm{Cov}(\beta \mid y)
# #'       = E[\sigma^2 \mid y]\,
# #'         \left(\Sigma_0^{-1} + X^\top W X\right)^{-1}.
# #'   }
# #'
# #' \strong{Weak-prior limits (Theorems 2 and 3)}
# #'
# #' As \eqn{n_{\mathrm{prior}} \to 0^+} (equivalently \eqn{\mathrm{pwt} \to 0}),
# #' \eqn{S_{\mathrm{marg}} \to \mathrm{RSS}_w}, and the conjugate Normal–Gamma
# #' posterior converges to the classical weighted least-squares limit:
# #' \deqn{
# #'   E[\beta \mid y] \to \hat\beta,\qquad
# #'   E[\sigma^2 \mid y] \to \frac{\mathrm{RSS}_w}{n_{\mathrm{effective}} - p},\qquad
# #'   \mathrm{Cov}(\beta \mid y) \to
# #'     \frac{\mathrm{RSS}_w}{n_{\mathrm{effective}} - p}
# #'     (X^\top W X)^{-1}.
# #' }
# #'
# #' For the independent Normal–Gamma prior used by
# #' \code{dIndependent_Normal_Gamma()}, neither the posterior mean nor the
# #' posterior covariance is available in closed form; the posterior must be
# #' obtained by numerical integration or sampling (e.g.,
# #' \code{rindepNormalGamma_reg()}). Theorem 3 in
# #' \insertCite{glmbayesChapterA12}{glmbayes} shows that the ING posterior has
# #' the same weak-prior limit as the conjugate Normal–Gamma posterior:
# #' \deqn{
# #'   E[\beta \mid y] \to \hat\beta,\qquad
# #'   \mathrm{Cov}(\beta \mid y) \to
# #'     \frac{\mathrm{RSS}_w}{n_{\mathrm{effective}} - p}
# #'     (X^\top W X)^{-1}.
# #' }
# #'
# #' @return
# #' A list of class \code{"PriorSetup"} with components:
# #' \item{mu}{Prior mean vector (length equal to the number of coefficients).}
# #' \item{Sigma}{Coefficient-scale prior variance–covariance matrix.}
# #' \item{Sigma_0}{For \code{family = gaussian()} only: dispersion-independent
# #'   prior covariance on the precision-weighted coefficient scale (the
# #'   \code{Sigma_0} passed to \code{\link{compute_gaussian_prior}}). Under
# #'   scalar \code{pwt},
# #'   \eqn{\Sigma_0^{-1} = \frac{\mathrm{pwt}}{1-\mathrm{pwt}} X^\top W X}.}
# #' \item{dispersion}{Calibrated dispersion (Gaussian models only), equal to
# #'   \eqn{S_{\mathrm{marg}}/(n_{\mathrm{effective}} - p)} under the default
# #'   calibration.}
# #' \item{shape}{Derived prior Gamma shape parameter for the Normal–Gamma prior
# #'   on precision (Gaussian only), \eqn{(n_{\mathrm{prior}} + k)/2}.}
# #' \item{shape_ING}{For \code{gaussian()} only when \code{shape} is available:
# #'   dedicated shape parameter for \code{dIndependent_Normal_Gamma()},
# #'   \eqn{\texttt{shape} + p/2}.}
# #' \item{rate}{Derived prior Gamma rate parameter (Gaussian only), using the
# #'   calibrated \eqn{S_{\mathrm{marg}}}.}
# #' \item{rate_gamma}{For \code{gaussian()} only, when Gaussian calibration runs:
# #'   prior Gamma rate for \code{dGamma()} / fixed-\eqn{\beta} use, based on
# #'   \eqn{\mathrm{RSS}_w(\beta_\star)} at the Zellner blend.}
# #' \item{coefficients}{Named numeric vector of returned coefficient values.
# #'   For \code{gaussian()} with scalar or vector \code{pwt}, this is the
# #'   closed-form posterior-mean blend
# #'   \eqn{(1-\mathrm{pwt})\hat\beta + \mathrm{pwt}\mu} when inputs are valid;
# #'   otherwise it falls back to the full-model GLM coefficients.}
# #' \item{model}{The model frame used to construct the design matrix (if
# #'   \code{model = TRUE}).}
# #' \item{x}{The model matrix used (if \code{x = TRUE}).}
# #' \item{y}{The response vector used (if \code{y = TRUE}).}
# #' \item{call}{The matched call to \code{Prior_Setup()}.}
# #' \item{PriorSettings}{A list containing prior configuration details, including
# #'   \code{pwt}, \code{n_prior}, \code{n_effective}, \code{n_likelihood},
# #'   \code{intercept_source}, and \code{effects_source}.}
# #' \item{conj_poisson}{\code{NULL} unless \code{family = poisson(link = "identity")} with a
# #'   single-column design (intercept-only), scalar \code{pwt}, finite positive scalar
# #'   \code{n_prior}, and negligible \code{offset}.  Then a named list with conjugate
# #'   \code{Gamma(shape, rate)} hyperparameters on the Poisson rate, with \eqn{\texttt{shape} =
# #'   n_{\mathrm{prior}}\,\bar{y}_w}, \eqn{\texttt{rate} = n_{\mathrm{prior}}}, and \code{beta}
# #'   centered at the weighted sample mean \eqn{\bar{y}_w}; pass to \code{\link{dGamma}} with \code{Inv_Dispersion = FALSE}.
# #'   See Details.}
# #'
# #' @family prior
# #' @seealso
# #' \code{\link{pfamily}} for prior-family objects and the constructors
# #' \code{\link{dNormal}}, \code{\link{dNormal_Gamma}}, \code{\link{dGamma}},
# #' and \code{\link{dIndependent_Normal_Gamma}}.
# #'
# #' \code{\link{glmb}}, \code{\link{lmb}} for formula-based fits with a
# #' \code{pfamily} built from \code{Prior_Setup()} output; \code{\link{rglmb}},
# #' \code{\link{rlmb}} for matrix-based sampling that consumes the same prior
# #' structure; \code{\link{simfuncs}} for functions that take a \code{prior_list}
# #' assembled from those components (including
# #' \code{\link{rindepNormalGamma_reg}} for
# #' \code{\link{dIndependent_Normal_Gamma}()}).
# #' \code{\link{multi_prior_setup}} for a matrix/cbind response with Gaussian;
# #' use with \code{\link{lmb}}
# #' \code{Prior_Setup} per column. 
# #'
# #' \insertCite{zellner1986gprior}{glmbayes};
# #' \insertCite{Raiffa1961}{glmbayes};
# #' \insertCite{Gelman2013}{glmbayes};
# #' \insertCite{McCullagh1989}{glmbayes};
# #' \insertCite{glmbayesChapter03}{glmbayes};
# #' \insertCite{glmbayesChapterA12}{glmbayes}.
# #'
# #' @references
# #' \insertAllCited{}
# #'
# #' @example inst/examples/Ex_Prior_Setup.R
# (Prior_Setup doc archive)
# NULL

## glmbayes implementation commented out; exported from glmbayesCore (R/reexports.R).
# ## Note arguments outside of first two are currently not used
#
# Prior_Setup <- function(
#     formula,
#     family      = gaussian(),
#     data=NULL,
#     weights=NULL,
#     subset=NULL,
#     na.action   = na.fail,
#     offset=NULL,
#     contrasts   = NULL,
#     pwt         = NULL,
#     pwt_default_low = 0.01,      # new: low-d default
#     pwt_default_high = 0.05,     # new: high-d default
#     n_prior     = NULL,
#     sd          = NULL,
#     dispersion  = NULL,
#     intercept_source = c("null_model", "full_model"),
#     effects_source   = c("null_effects",  "full_model"),
#     mu          = NULL,
#     k           = 1,
#     ...
#   ) 
#   
#   {
# 
#   ## ---------------------------------------------------------------------------
#   ## Step 1: Parse and normalize top-level arguments.
#   ## ---------------------------------------------------------------------------
#   call <- match.call()  
#   intercept_source <- match.arg(intercept_source)
#   effects_source <- match.arg(effects_source)
#   if (!is.null(dispersion)) {
#     if (!is.numeric(dispersion) || length(dispersion) != 1L ||
#         !is.finite(dispersion) || dispersion <= 0) {
#       stop("dispersion must be NULL or a single positive finite numeric value.")
#     }
#   }
#   dispersion_input <- dispersion
#   
#   
#   
#   #mf<-model.frame(formula,data,subset=subset,na.action=na.action,
#   #                drop.unused.levels=drop.unused.levels,xlev=xlev)
#   
#   
#   if (is.character(family)) 
#     family <- get(family, mode = "function", envir = parent.frame())
#   if (is.function(family)) 
#     family <- family()
#   if (is.null(family$family)) {
#     print(family)
#     stop("'family' not recognized")
#   }
#   
#   if (missing(data))   data <- environment(formula)
#   
#   ## ---------------------------------------------------------------------------
#   ## Step 2: Build model frame / response / design matrix.
#   ## ---------------------------------------------------------------------------
#   mf <- match.call(expand.dots = FALSE)
#   m <- match(c("formula", "data", "subset", "weights", "na.action", 
#                "etastart", "mustart", "offset"), names(mf), 0L)
#   mf <- mf[c(1L, m)]
#   mf$drop.unused.levels <- TRUE
#   mf[[1L]] <- quote(stats::model.frame)
#   mf <- eval(mf, parent.frame())
#   
##  mf<-model.frame(formula,data)
#   
#   mt <- attr(mf, "terms")
#   Y <- model.response(mf, "any")
#   
#   if (family$family == "binomial"
#       && is.numeric(Y) && is.vector(Y)
#       && all(Y >= 0 & Y <= 1)    # use <= instead of <=
#       && is.null(weights)) {
#     warning(
#       "You supplied a proportion response (0 <= y <= 1) to a binomial family\n",
#       "without `weights`. Each case will be treated as a single trial (n=1).\n",
#       "If you meant to model counts, either use `cbind(success, failure)`\n",
#       "or supply `weights =` the number of trials."
#     )
#   }
#   
#   
##  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(, NROW(Y), 0L)
#   X <- if (!is.empty.model(mt)) model.matrix(mt, mf, ...) else matrix(, NROW(Y), 0L)
#   
#   ## ---------------------------------------------------------------------------
#   ## Step 3: Resolve weights and effective sample size.
#   ## ---------------------------------------------------------------------------
#   ## --- WEIGHT HANDLING -------------------------------------------------------
#   
#   # Extract raw weights from the model frame (may be NULL)
#   raw_wt <- model.weights(mf)
#   
#   # Number of observations (always correct)
#   n_obs <- nrow(X)
#   
#   # Case 1: User supplied scalar weight (e.g., weights = 4)
#   # model.frame() cannot accept scalar weights, so we expand them *after* mf is built
#   if (!is.null(raw_wt) && length(raw_wt) == 1L) {
#     weights <- rep(raw_wt, n_obs)
#   }
#   
#   # Case 2: User supplied a full-length weight vector
#   else if (!is.null(raw_wt)) {
#     
#     if (!is.numeric(raw_wt))
#       stop("'weights' must be numeric")
#     
#     if (any(raw_wt < 0))
#       stop("negative weights not allowed")
#     
#     if (length(raw_wt) != n_obs)
#       stop("weights must be either a scalar or have length equal to number of observations")
#     
#     weights <- raw_wt
#   }
#   
#   # Case 3: No weights supplied -> default depends on family
#   else {
#     
#     if (family$family %in% c("gaussian", "Gamma")) {
#       # Gaussian/Gamma treat weights as replication weights -> default = 1
#       weights <- rep(1, n_obs)
#     } else {
#       # Binomial/Poisson: GLM semantics require weights = NULL
#       weights <- NULL
#     }
#   }
#   
#   # Compute effective sample size
#   if (is.null(weights)) {
#     n_effective <- n_obs
#   } else {
#     n_effective <- sum(weights)
#   }
#   
#   
#   #######################################33
#   
#   
#   
#   
#   
#   
#   offset <- as.vector(model.offset(mf))
#   if (!is.null(offset)) {
#     if (length(offset) != NROW(Y)) 
#       stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
#                     length(offset), NROW(Y)), domain = NA)
#   }
#   
#   mustart <- model.extract(mf, "mustart")
#   etastart <- model.extract(mf, "etastart")
#   
#   
#   x<-X  
##  x<-model.matrix(formula,mf)
#   
#   nvar=ncol(x)
#   if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 0) {
#     stop("k must be a single non-negative finite numeric value.", call. = FALSE)
#   }
#   if (k + nvar < 2) {
#     stop(
#       "Prior_Setup: require k + p >= 2, where p is the number of coefficients (ncol of model matrix). ",
#       "Got k = ", k, ", p = ", nvar, ".",
#       call. = FALSE
#     )
#   }
#   
#   ## ---------------------------------------------------------------------------
#   ## Step 4: Resolve prior-weight inputs (pwt, sd, n_prior).
#   ## Shared by Gaussian and non-Gaussian families.
#   ## ---------------------------------------------------------------------------
#   if (is.null(pwt)) {
#     pwt <- if (nvar < 14) pwt_default_low else pwt_default_high
#     ## n_prior (later) recomputes pwt; avoid implying the default applies.
#     if (is.null(n_prior)) {
#       message("Using default pwt = ", pwt,
#               " (", if (nvar < 14) "low-d" else "high-d", " default).")
#     }
#   }
#   
#   ## Make sure the *columns* of x are named correctly:
#   
#   
#   ## validate pwt  
#   if (!is.numeric(pwt) || any(is.na(pwt))) {  
#     stop("pwt must be numeric and non NA, either length 1 or length ", nvar)  
#   }  
#   if (! (length(pwt) %in% c(1, nvar)) ) {  
#     stop("pwt must have length 1 or length(coef) = ", nvar,  
#          "; you supplied length ", length(pwt))  
#   }  
#   if (any(pwt <= 0 | pwt >= 1)) {  
#     stop("All elements of pwt must lie strictly between 0 and 1; you supplied:",  
#          paste0(round(pwt, 3), collapse = ", "))  
#   }  
#   
#   
#   
#   var_names <- colnames(x)
#   colnames(x) <- var_names
# 
#   mu_internal <- matrix(0, nrow = nvar, ncol = 1, dimnames = list(var_names, "mu"))
#  
#   
#   ## ---------------------------------------------------------------------------
#   ## Step 5: Fit full GLM and extract baseline covariance (V0).
#   ## ---------------------------------------------------------------------------
#   
# 
#   glm_full <- glm.fit(
#     x       = X,
#     y       = Y,
#     weights = weights,
#     offset  = offset,
#     family  = family
#     ,control = glm.control(...)
#   )
#   
# 
#   glm_full$call      <- call
#   glm_full$formula   <- formula
#   glm_full$terms     <- mt
#   glm_full$data      <- mf
#   glm_full$offset    <- offset
#   glm_full$contrasts <- attr(X, "contrasts")
#   glm_full$xlevels   <- .getXlevels(mt, mf)
#   class(glm_full)    <- c("glm", "lm")
#   
# 
#   V0 <- vcov(glm_full)
#   ## Saturated Gaussian blocks (n = p, RSS = 0): vcov(glm) is NA; use dispersion * (X'WX)^{-1}.
#   if (anyNA(V0) && family$family %in% c("gaussian", "poisson")) {
#     d_v0 <- if (!is.null(dispersion_input)) {
#       dispersion_input
#     } else if (identical(family$family, "poisson")) {
#       1
#     } else {
#       NULL
#     }
#     if (!is.null(d_v0)) {
#     w_fit <- if (is.null(weights)) rep(1, n_obs) else as.numeric(weights)
#     XtW <- sweep(x, 1, w_fit, `*`)
#     Gm <- crossprod(XtW, x)
#     Ginv <- tryCatch(
#       solve(Gm),
#       error = function(e) {
#         stop(
#           "vcov(glm_full) is NA and (X'WX) is singular; the design is rank-deficient.",
#           call. = FALSE
#         )
#       }
#     )
#     V0 <- d_v0 * Ginv
#     dimnames(V0) <- list(var_names, var_names)
#     }
#   }
# 
#   glm_summary=summary(glm_full)
#   
#   ##n_likelihood <- glm_summary$df.residual + glm_summary$df[1]  # residual df + model rank
#   n_likelihood <- n_effective
#   
#   # If sd is provided, use it to compute pwt
# if (!is.null(sd)) {
#   if (!is.numeric(sd) || any(is.na(sd))) {
#     stop("sd must be a numeric vector with no missing values.")
#   }
#   if (length(sd) != nvar) {
#     stop("Length of sd must match number of coefficients (", nvar, ").")
#   }
# 
#   # Compute pwt from sd and V0
#   V0_diag <- diag(V0)
#   if (any(V0_diag <= 0)) {
#     stop("Diagonal entries of V0 must be positive to compute pwt from sd.")
#   }
# 
#   pwt <- V0_diag / (V0_diag + sd^2)
#   message("Computed pwt from user-specified prior standard deviations (sd).")
# }
#   
#     
#   ## n_prior may imply pwt only when pwt is still a single scalar not from `sd`.
#   ## If length(pwt) > 1 (vector pwt) or `sd` was supplied, do not overwrite pwt;
#   ## n_prior is then used only downstream as a scalar for Gamma / S_marg remap.
#   if (!is.null(n_prior)) {
#     if (!is.numeric(n_prior) || length(n_prior) != 1 || n_prior <= 0) {
#       stop("n_prior must be a single positive numeric value")
#     }
#     if (length(pwt) == 1L && is.null(sd)) {
##    pwt <- n_prior / (n_prior + n_likelihood)
#       pwt <- n_prior / (n_prior + n_effective)
#       message("Computed pwt = ", round(pwt, 4),
#               " from n_prior = ", n_prior,
#               " and n_effective = ", n_effective)
#     }
#   }
#   
#   # Compute n_prior if not supplied and pwt is scalar
#   if (is.null(n_prior) && length(pwt) == 1L) {
#     n_prior <- (pwt/(1-pwt)) * n_effective
#   #  message("Computed n_prior = ", round(n_prior, 4),
#   #          " from pwt = ", round(pwt, 4),
#   #          " and n_likelihood = ", n_likelihood)
#   }
#   if (identical(family$family, "gaussian") && is.null(n_prior)) {
#     stop(
#       "For Gaussian models, a scalar effective prior sample size `n_prior` is required. ",
#       "Use scalar `pwt` (it implies `n_prior`), or supply `n_prior` explicitly. ",
#       "Per-coefficient `sd` implies vector `pwt`; in that case you must pass `n_prior`.",
#       call. = FALSE
#     )
#   }
# 
#   ## ---------------------------------------------------------------------------
#   ## Step 6: Family-specific dispersion baseline.
#   ## Gaussian: weighted RSS ratio; Gamma: MASS::gamma.dispersion; else NULL.
#   ## ---------------------------------------------------------------------------
#   ## --- CONDITIONAL DISPERSION (Gaussian): explicit ratio from glm.fit object -----
#   ## Uses stats::glm.fit (not glm()). Baseline dispersion = RSS_w / (n_effective - p),
#   ## p = ncol(x), matching weighted residual df and \code{\link{compute_gaussian_prior}}.
#   ## # Old MLE-style ratio (retained for reference, not used):
#   ## # dispersion <- rss_weighted / n_effective
#   ## With rate = dispersion * shape (shape = (n_prior + k) / 2), posterior summaries of tau = 1/d
#   ## depend on pwt unless additional structure holds (see Details).
#   rss_weighted_stored <- NA_real_
#   dispersion_classical <- NA_real_
#   if (family$family == "gaussian") {
#     res <- residuals(glm_full, type = "response")
#     w   <- glm_full$prior.weights
#     rss_weighted <- sum(w * res^2)
#     rss_weighted_stored <- rss_weighted
#     if (!is.null(dispersion_input)) {
#       dispersion <- dispersion_input
#       dispersion_classical <- dispersion_input
#     } else {
#       if (!is.finite(rss_weighted) || rss_weighted <= 0) {
#         stop("Weighted RSS must be strictly positive for Gaussian dispersion priors.")
#       }
#       if (!is.finite(n_effective) || n_effective <= 0) {
#         stop("n_effective must be strictly positive to compute Gaussian dispersion.")
#       }
#       if (n_effective <= nvar) {
#         stop(
#           "Gaussian dispersion requires n_effective > p (number of coefficients); ",
#           "use denominator n_effective - p, or pass a positive scalar `dispersion`. ",
#           "Got n_effective = ", n_effective, ", p = ", nvar, ".",
#           call. = FALSE
#         )
#       }
#       dispersion <- rss_weighted / (n_effective - nvar)
#       if (!is.finite(dispersion) || dispersion <= 0) {
#         stop("Computed Gaussian dispersion must be strictly positive.")
#       }
#       dispersion_classical <- dispersion
#     }
# 
#   } else if (family$family == "Gamma") {
#     
#     # MASS::gamma.dispersion() already returns the correct quasi-likelihood
#     # dispersion estimate for Gamma GLMs.
#     dispersion <- MASS::gamma.dispersion(glm_full)
#     
#   } else {
#     
#     dispersion <- NULL
#   }
#   
# 
#     if (!is.matrix(V0) || nrow(V0) != ncol(V0)) {
#     stop("vcov(glm_full) (V0) must be a square matrix.")
#   }
#   if (anyNA(V0)) {
#     stop("vcov(glm_full) (V0) contains missing values.")
#   }
#   
#   # 2. symmetry (up to numerical tolerance)
#   if (!isSymmetric(V0, tol = sqrt(.Machine$double.eps))) {
#     stop("vcov(glm_full) (V0) is not symmetric.")
#   }
#   
#   # 3. positive-definiteness via Cholesky
#   pd_try <- try(chol(V0), silent = TRUE)
#   if (inherits(pd_try, "try-error")) {
#     stop(
#       "Variance-covariance matrix V0 is not positive-definite.\n",
#       "This usually means the classical GLM is rank-deficient."
#     )
#   }
# 
#   ## ---------------------------------------------------------------------------
#   ## Step 7: Construct prior mean vector mu.
#   ## Intercept and effects can be sourced from null/full model unless user sets mu.
#   ## ---------------------------------------------------------------------------
#   if (var_names[1] == "(Intercept)") {
#     # build 1-column design matrix for intercept only
#     X0 <- matrix(1, nrow = NROW(Y), ncol = 1,
#                  dimnames = list(NULL, "(Intercept)"))
#     
#     # fit intercept-only model via glm.fit()
#     fit0 <- glm.fit(
#       x       = X0,
#       y       = Y,
#       weights = weights,
#       offset  = offset,
#       family  = family,
#       control = glm.control(...)
#     )
#     
#     # pick the intercept from null or full model
#     chosen_int <- switch(
#       intercept_source,
#       null_model = fit0$coefficients[1],
#       full_model = glm_full$coefficients[1]
#     )
#     
#     mu_internal[1, 1] <- chosen_int
#   }
#   
# 
#   # 5) effects prior means
#   if (nvar > 1) {
#     effect_names <- var_names[-1]
#     if (effects_source == "full_model") {
#       coefs <- coef(glm_full)[effect_names]
#       mu_internal[effect_names, 1] <- coefs
#     }
#     # else null_effects leaves mu[...] as zero
#   }
#   
# 
#   # Validate user-supplied mu if provided
#   if (!is.null(mu)) {
#     if (!is.numeric(mu)) {
#       stop("mu must be numeric.")
#     }
#     if (is.vector(mu)) {
#       if (length(mu) != nvar) {
#         stop("Length of mu vector must match number of coefficients (", nvar, ").")
#       }
#       mu <- matrix(mu, ncol = 1, dimnames = list(var_names, "mu"))
#     } else if (is.matrix(mu)) {
#       if (!all(dim(mu) == c(nvar, 1))) {
#         stop("mu matrix must have dimensions [", nvar, ", 1].")
#       }
#       rownames(mu) <- var_names
#       colnames(mu) <- "mu"
#     } else {
#       stop("mu must be either a numeric vector or a matrix.")
#     }
#     message("Using user-specified prior mean vector (mu).")
#   } else {
#     mu <- mu_internal
#   }
#   
#     
#   ## ---------------------------------------------------------------------------
#   ## Step 8: Build prior covariance Sigma from pwt and V0.
#   ## Scalar pwt gives Zellner scaling; vector pwt applies element-wise scaling.
#   ## ---------------------------------------------------------------------------
#   
#  Sigma=as.matrix(diag(nvar))
#  
#  Sigma=(1-pwt)/pwt*V0
#    
#   ## build prior covariance  
#   if (length(pwt) == 1L) {  
#     ## full matrix prior  
#     Sigma <- ((1 - pwt) / pwt) * V0  
#   }
#   else {  
#     scale_vec <- sqrt((1 - pwt) / pwt)
#     scale_mat <- outer(scale_vec, scale_vec)
#     Sigma <- V0 * scale_mat
#   }  
# 
#   rownames(mu)=var_names
#   colnames(mu)=c("mu")
#   rownames(Sigma)=var_names
#   colnames(Sigma)=var_names
# 
#   rate_gamma <- NULL
#   shape_ING <- NULL
#   ## ---------------------------------------------------------------------------
#   ## Step 9: Build Gamma(shape, rate) hyperparameters when available.
#   ## For Gaussian this provides precision-prior terms used in calibration.
#   ## ---------------------------------------------------------------------------
#   ## Gamma on precision: shape = (n_prior + 1) / 2, rate = dispersion * (n_prior/2).
#   ## compute_gaussian_prior() calibrates shape/rate from n_prior and S_marg only.
#   ## Sigma may be rescaled below by Gaussian calibration; shape/rate use current dispersion.
#   dispersion_for_shape_rate <- dispersion
#   if (!is.null(n_prior) && length(n_prior) == 1L && !is.null(dispersion_for_shape_rate)) {
#     ## n_prior is interpreted as effective prior sample size, on the same scale as sum(weights).
#     shape <- (n_prior + 1L) / 2
#     if (!is.finite(shape) || shape <= 0) {
#       stop("Computed shape must be strictly positive.")
#     }
#     rate <- dispersion_for_shape_rate * (n_prior / 2)
#     if (!is.finite(rate) || rate <= 0) {
#       stop("Computed rate must be strictly positive.")
#     }
#   } else {
#     shape <- NULL
#     rate <- NULL
#   }
# 
#   ## ---------------------------------------------------------------------------
#   ## Step 10 (Gaussian): run Gaussian calibration pipeline and replace outputs.
#   ## compute_gaussian_prior() returns calibrated dispersion/shape/rate/Sigma.
#   ## ---------------------------------------------------------------------------
#   ## --- S_marg_new / S_marg_sigma0_vcov before Post_mean / Nelder (full Sigma_pre_nm: scalar or vector pwt)
#   ## Sigma_0 = Sigma_pre_nm / d with d = d_OLS or d_vcov = summary(glm)$dispersion (cancels vcov scale in V0).
#   ## S_marg keeps the same value for downstream b_0 / remap logic (identical to S_marg_new here).
#   Sigma_pre_nm <- Sigma
#   .gauss_helper_preview <- NULL
#   if (identical(family$family, "gaussian") &&
#       n_effective > nvar &&
#       is.finite(dispersion_classical) && dispersion_classical > 0 &&
#       !is.null(n_prior) && length(n_prior) == 1L && is.finite(n_prior) && n_prior > 0 &&
#       !is.null(mu) && length(as.numeric(mu)) == nvar && all(is.finite(as.numeric(mu)))) {
#     w_h <- if (is.null(weights)) rep(1, n_obs) else as.numeric(weights)
#     off_h <- if (is.null(offset)) rep(0, n_obs) else as.numeric(offset)
#     bhat_h <- coef(glm_full)
#     if (length(bhat_h) == nvar && all(is.finite(bhat_h))) {
#       Sigma_0_h <- Sigma_pre_nm / dispersion_classical
#       .gauss_helper_preview <- compute_gaussian_prior(
#         X = X,
#         Y = Y,
#         weights = w_h,
#         offset = off_h,
#         dispersion = dispersion_input,
#         n_effective = n_effective,
#         bhat = bhat_h,
#         mu = mu,
#         Sigma_0 = Sigma_0_h,
#         Sigma = if (!is.null(sd)) Sigma_pre_nm else NULL,
#         n_prior = n_prior,
#         k = k
#       )
#     }
#   }
#   coefficients_mle <- coef(glm_full)
#   coefficients <- coefficients_mle
#   ## Default returned coefficients: closed-form posterior mean blend when available.
#   if (identical(family$family, "gaussian") &&
#       length(coefficients_mle) == nvar &&
#       !is.null(mu) && length(mu) == nvar &&
#       all(is.finite(as.numeric(mu)))) {
#     mle_fp <- vapply(
#       var_names,
#       function(nm) {
#         if (!is.null(names(coefficients_mle)) && nm %in% names(coefficients_mle)) {
#           v <- unname(coefficients_mle[nm])
#           if (length(v) == 1L && is.finite(v)) v else NA_real_
#         } else {
#           NA_real_
#         }
#       },
#       NA_real_
#     )
#     mu_fp <- as.numeric(mu)
#     if (length(pwt) == 1L && is.finite(pwt)) {
#       coefficients <- (1 - pwt) * mle_fp + pwt * mu_fp
#       names(coefficients) <- var_names
#     } else if (length(pwt) == nvar && all(is.finite(pwt))) {
#       coefficients <- (1 - pwt) * mle_fp + pwt * mu_fp
#       names(coefficients) <- var_names
#     }
#   }
# 
#   ## When calibration ran, take Gaussian dispersion and Gamma hyperparameters from compute_gaussian_prior().
#   if (identical(family$family, "gaussian") &&
#       !is.null(.gauss_helper_preview)) {
#     dispersion <- .gauss_helper_preview$dispersion
#     shape <- .gauss_helper_preview$shape
#     shape_ING <- .gauss_helper_preview$shape_ING
#     rate <- .gauss_helper_preview$rate
#     rate_gamma <- .gauss_helper_preview$rate_gamma
#     Sigma <- .gauss_helper_preview$Sigma
#     rownames(Sigma) <- var_names
#     colnames(Sigma) <- var_names
#   }
#   if (identical(family$family, "gaussian") &&
#       !is.null(shape) && length(shape) == 1L && is.finite(shape)) {
#     if (is.null(shape_ING)) {
#       shape_ING <- shape + nvar / 2
#     }
#   }
# 
#   Sigma_0_out <- NULL
#   if (identical(family$family, "gaussian")) {
#     if (!is.null(.gauss_helper_preview)) {
#       Sigma_0_out <- .gauss_helper_preview$Sigma_0
#     } else if (is.finite(dispersion_classical) && dispersion_classical > 0) {
#       Sigma_0_out <- Sigma_pre_nm / dispersion_classical
#     }
#     if (!is.null(Sigma_0_out)) {
#       rownames(Sigma_0_out) <- var_names
#       colnames(Sigma_0_out) <- var_names
#     }
#   }
# 
#   ## ---------------------------------------------------------------------------
#   ## Step 10b: Gamma–Poisson conjugate λ prior (intercept-only Poisson(link = "identity"))
#   ##
#   ## Effective prior observation count links to scalar pwt via
#   ##   n_prior / (n_prior + n_eff) = pwt  <=>  n_prior = (pwt / (1 - pwt)) * n_eff
#   ## Weighted arithmetic mean λ̄_w = Σ w_i y_i / Σ w_i.  Choosing
#   ##   λ ~ Gamma(shape = n_prior λ̄_w, rate = n_prior),
#   ## gives E[λ] = λ̄_w so the conjugate Poisson posterior mean remains λ̄_w.
#   ## ---------------------------------------------------------------------------
# 
#   conj_poisson <- NULL
#   mu_arg_missing <- missing(mu)
# 
#   gamma_poisson_conj_ok <- (
#     identical(family$family, "poisson") &&
#       identical(family$link, "identity") &&
#       ncol(x) == 1L &&
#       length(pwt) == 1L
#   )
# 
#   if (gamma_poisson_conj_ok && !is.null(offset)) {
#     of <- suppressWarnings(as.numeric(offset))
#     if (length(of) != n_obs || any(!is.finite(of)) ||
#         max(abs(of), na.rm = TRUE) > sqrt(.Machine$double.eps)) {
#       gamma_poisson_conj_ok <- FALSE
#     }
#   }
# 
#   if (gamma_poisson_conj_ok) {
#     ww <- glm_full$prior.weights
#     if (is.null(ww) || length(ww) != n_obs) ww <- rep(1, n_obs)
# 
#     ## Nonnegative Poisson observations
#     if (any(as.numeric(Y) < 0, na.rm = TRUE)) {
#       stop(
#         "`Prior_Setup()` Poisson conjugate calibration requires nonnegative responses.", call. = FALSE
#       )
#     }
# 
#     ybar <- sum(as.numeric(Y) * as.numeric(ww)) / sum(as.numeric(ww))
#     if (!(is.finite(ybar) && ybar > 0)) {
#       stop(
#         "`Prior_Setup()` Poisson(link='identity'), intercept-only conjugate calibration requires ",
#         "positive weighted mean counts (Gamma prior requires a positive Poisson rate); weighted mean ",
#         "= ", paste0("`", prettyNum(ybar), "`."), call. = FALSE
#       )
#     }
# 
#     if (is.null(n_prior) || length(n_prior) != 1L || !is.finite(n_prior) || n_prior <= 0) {
#       warning(
#         "Poisson(link='identity') intercept-only: skipping `conj_poisson` Gamma rate prior ",
#         "calibration because scalar `n_prior` is unavailable or non-positive; supply finite positive ",
#         "`pwt` or `n_prior`.",
#         call. = FALSE
#       )
#     } else {
#       np <- as.numeric(n_prior)
#       ## Prior mean E[lambda] matches weighted data mean; conjugate Gamma(shape,rate) hyperparameters:
#       conj_shape <- np * ybar
#       conj_rate <- np
# 
#       bm <- matrix(as.numeric(ybar), nrow = 1L, ncol = 1L,
#                    dimnames = list(NULL, var_names))
# 
#       conj_poisson <- list(
#         shape = conj_shape,
#         rate = conj_rate,
#         beta = bm,
#         weighted_mean_rate = ybar,
#         n_prior_eff = np
#       )
# 
#       ## Moment-matched surrogate Normal summaries for `\u03bc`/`Sigma` (single intercept only):
#       if (mu_arg_missing && nvar == 1L && is.finite(ybar) && is.finite(np) && np > 0) {
#         mu[1, 1] <- as.numeric(ybar)
#         Sigma[1, 1] <- as.numeric(ybar / np)
#         rownames(mu) <- var_names
#         colnames(mu) <- "mu"
#         rownames(Sigma) <- colnames(Sigma) <- var_names
#       }
#     }
#   }
# 
#   ## ---------------------------------------------------------------------------
#   ## Step 10c: Beta–Binomial conjugate θ prior (intercept-only Binomial(link = "identity"))
#   ##
#   ## Effective prior observation count:
#   ##   n_prior / (n_prior + n_eff) = pwt  <=>  n_prior = (pwt / (1 - pwt)) * n_eff
#   ## Weighted proportion p̄_w = Σ w_i y_i / Σ w_i.  Choosing
#   ##   θ ~ Beta(shape1 = n_prior * p̄, shape2 = n_prior * (1 - p̄))
#   ## gives E[θ] = p̄ so the conjugate Binomial posterior mean starts at p̄.
#   ## ---------------------------------------------------------------------------
# 
#   conj_binomial <- NULL
# 
#   beta_binom_conj_ok <- (
#     identical(family$family, "binomial") &&
#       identical(family$link, "identity") &&
#       ncol(x) == 1L &&
#       length(pwt) == 1L
#   )
# 
#   if (beta_binom_conj_ok && !is.null(offset)) {
#     of <- suppressWarnings(as.numeric(offset))
#     if (length(of) != n_obs || any(!is.finite(of)) ||
#         max(abs(of), na.rm = TRUE) > sqrt(.Machine$double.eps)) {
#       beta_binom_conj_ok <- FALSE
#     }
#   }
# 
#   if (beta_binom_conj_ok) {
#     ww <- glm_full$prior.weights
#     if (is.null(ww) || length(ww) != n_obs) ww <- rep(1, n_obs)
# 
#     y_num <- as.numeric(Y)
#     if (any(y_num < 0, na.rm = TRUE) || any(y_num > 1, na.rm = TRUE)) {
#       stop(
#         "`Prior_Setup()` Binomial(identity) conjugate calibration requires response in [0, 1].",
#         call. = FALSE
#       )
#     }
# 
#     pbar <- sum(y_num * as.numeric(ww)) / sum(as.numeric(ww))
#     if (!(is.finite(pbar) && pbar > 0 && pbar < 1)) {
#       warning(
#         "Binomial(link='identity') intercept-only: skipping `conj_binomial` Beta prior ",
#         "calibration because the weighted proportion is not strictly in (0, 1); ",
#         "weighted proportion = ", paste0("`", prettyNum(pbar), "`."),
#         call. = FALSE
#       )
#     } else if (is.null(n_prior) || length(n_prior) != 1L || !is.finite(n_prior) || n_prior <= 0) {
#       warning(
#         "Binomial(link='identity') intercept-only: skipping `conj_binomial` Beta prior ",
#         "calibration because scalar `n_prior` is unavailable or non-positive; supply finite positive ",
#         "`pwt` or `n_prior`.",
#         call. = FALSE
#       )
#     } else {
#       np <- as.numeric(n_prior)
#       conj_shape1 <- np * pbar
#       conj_shape2 <- np * (1 - pbar)
# 
#       bm <- matrix(as.numeric(pbar), nrow = 1L, ncol = 1L,
#                    dimnames = list(NULL, var_names))
# 
#       conj_binomial <- list(
#         shape1             = conj_shape1,
#         shape2             = conj_shape2,
#         beta               = bm,
#         weighted_mean_prop = pbar,
#         n_prior_eff        = np
#       )
# 
#       ## Moment-matched surrogate Normal mu/Sigma for single intercept
#       if (mu_arg_missing && nvar == 1L && is.finite(pbar) && is.finite(np) && np > 0) {
#         mu[1, 1] <- as.numeric(pbar)
#         Sigma[1, 1] <- as.numeric(pbar * (1 - pbar) / np)
#         rownames(mu) <- var_names
#         colnames(mu) <- "mu"
#         rownames(Sigma) <- colnames(Sigma) <- var_names
#       }
#     }
#   }
# 
#   ## ---------------------------------------------------------------------------
#   ## Step 11: Assemble and return PriorSetup object.
#   ## ---------------------------------------------------------------------------
#   prior_list <- list(
#     mu = mu,
#     Sigma = Sigma,
#     Sigma_0 = Sigma_0_out,
#     dispersion = dispersion,
#     shape = shape,
#     shape_ING = shape_ING,
#     rate = rate,
#     rate_gamma = rate_gamma,
#     coefficients = coefficients,
#     conj_poisson  = conj_poisson,
#     conj_binomial = conj_binomial,
#     model = mf,
#     x = x,
#     y = Y,
#     call = call,
#     PriorSettings = list(
#       pwt = pwt,
#       n_prior = n_prior,
#       intercept_source = intercept_source,
#       effects_source = effects_source,
#       ## For now retain n_likelihood for backward compatibility
#       n_likelihood = n_likelihood,
#       n_effective = n_effective
#     )
#   )
#   
#   class(prior_list) <- "PriorSetup"
#   return(prior_list)
#   
# }

## print.* S3 method provided by glmbayesCore (Imports); local registration commented out.
# #' @export
# #' @method print PriorSetup
# #' @rdname Prior_Setup
# 
# print.PriorSetup <- function(x, ...) {
# 
#   cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
#     
#   settings <- x$PriorSettings
#   
#   if (!is.null(settings$pwt) && length(settings$pwt) == 1L) {
#     g <- (1 - settings$pwt) / settings$pwt
#     cat("Setting up a Zellner g-type prior: \n")
#     cat("  pwt =", round(settings$pwt, 4), "\n")
#     cat("  g   = (1 - pwt)/pwt =", round(g, 4), "\n\n")
#   }
#   
#   if (!is.null(settings$n_prior) && !is.null(settings$n_likelihood)) {
#     if (!is.null(settings$pwt) && length(settings$pwt) > 1L) {
#       cat("Note: n_prior was provided to Prior_Setup (scalar prior sample size for precision / calibration; pwt is per coefficient):\n")
#     } else {
#       n_eff_print <- settings$n_effective
#       if (is.null(n_eff_print)) {
#         n_eff_print <- settings$n_likelihood
#       }
#       pwt_s <- settings$pwt
#       n_prior_implies_pwt <- FALSE
#       if (length(pwt_s) == 1L && is.finite(pwt_s) && pwt_s > 0 && pwt_s < 1 &&
#           is.finite(n_eff_print) && n_eff_print > 0 &&
#           is.finite(settings$n_prior)) {
#         n_prior_from_pwt <- (pwt_s / (1 - pwt_s)) * n_eff_print
#         np_stored <- as.numeric(settings$n_prior)
#         scale_np <- max(abs(np_stored), abs(n_prior_from_pwt), 1e-12)
#         n_prior_implies_pwt <- abs(np_stored - n_prior_from_pwt) <= 1e-10 * scale_np
#       }
#       if (n_prior_implies_pwt) {
#         cat("Note: n_prior was computed as (pwt / (1 - pwt)) * n_likelihood: \n")
#       } else {
#         cat("Note: n_prior was provided to Prior_Setup (scalar prior sample size for precision / calibration):\n")
#       }
#     }
#     cat("  n_prior      =", round(settings$n_prior, 4), "\n")
#     cat("  n_likelihood =", round(settings$n_likelihood, 4), "\n\n")
#   }
# 
#   if (!is.null(settings$pwt) && length(settings$pwt) > 1L) {
#     cat("Note: Differential prior weights (pwt) were specified per coefficient.\n\n")
#   }
#   
#   
#   cat("Prior Setup Summary\n")
#   cat("====================\n\n")
#   
#   # Check for Zellner g-prior structure
#   Sigma <- x$Sigma
#   mu <- x$mu
#   var_names <- rownames(mu)
#   nvar <- length(var_names)
#   
#   # Extract diagonal SDs
#   prior_sd <- sqrt(diag(Sigma))
#   
#   # Always compute prior correlation matrix
#   prior_cor <- cov2cor(Sigma)
#   
#   # Extract pwt vector for display
#   if (!is.null(settings$pwt)) {
#     if (length(settings$pwt) == 1L) {
#       pwt_vec <- rep(settings$pwt, nvar)
#     } else if (length(settings$pwt) == nvar) {
#       pwt_vec <- settings$pwt
#     } else {
#       warning("Length of pwt does not match number of variables; skipping pwt column.")
#       pwt_vec <- rep(NA_real_, nvar)
#     }
#   } else {
#     pwt_vec <- rep(NA_real_, nvar)
#   }
#   
# 
#   # Compute 95% intervals
#   z <- qnorm(0.975)
#   lower <- mu[, 1] - z * prior_sd
#   upper <- mu[, 1] + z * prior_sd
#   
# 
#   # Build output table
# 
#   out <- data.frame(
#     Prior.Mean = round(mu[, 1], 6),
#     Prior.SD   = round(prior_sd, 6),
#     CI.Lower   = round(lower, 6),
#     CI.Upper   = round(upper, 6),
#     pwt        = round(pwt_vec, 6)
#   )
#   
#   # Print table
#   cat("Prior Estimates with 95% Confidence Intervals\n")
#   invisible(print(out))
#   
#   if (nvar <= 10) {
#     cat("\nPrior Correlation Matrix\n")
#     invisible(print(round(prior_cor, 4)))
#   }
#   
#   # Optional: print dispersion
#   if (!is.null(x$dispersion)) {
#     cat("\nConditional Dispersion (Gaussian family): ", round(x$dispersion, 4), "\n\n")
#   }
# 
#   if (!is.null(x$shape) && !is.null(x$rate)) {
#     cat("Gamma Prior on Residual Precision:\n")
#     for (nm in c("shape", "rate")) {
#       z <- x[[nm]]
#       zch <- if (length(z) != 1L || !is.numeric(z)) {
#         as.character(z)
#       } else if (!is.finite(z)) {
#         as.character(z)
#       } else {
#         az <- abs(z)
#         if (az == 0) {
#           "0"
#         } else if (az < 1e-4) {
#           format(z, scientific = TRUE, digits = 6L)
#         } else {
#           as.character(round(z, 4L))
#         }
#       }
#       cat("  ", nm, " = ", zch, "\n", sep = "")
#     }
#     cat(
#       "  Expected precision (inverse variance) =",
#       format(signif(x$shape / x$rate, 6), scientific = TRUE, digits = 7),
#       "which implies 1/Expected precision  =",
#       format(signif(x$rate / x$shape, 6), scientific = TRUE, digits = 7),
#       "\n\n"
#     )
#     cat("  Applicable to gaussian models with compound pfamilies (e.g., dNormal_Gamma, dIndependent_Normal_Gamma),\n")
#     cat("  as well as for Gamma regression, quasipoisson, and quasibinomial models.\n\n")
#   }
# 
#   if (!is.null(x$conj_poisson)) {
#     cp <- x$conj_poisson
#     cat("Poisson rate prior (Gamma conjugate; pass to dGamma(Inv_Dispersion=FALSE)):\n")
#     cat(
#       "  shape = ",
#       round(as.numeric(cp$shape), 6),
#       "; rate = ",
#       round(as.numeric(cp$rate), 6),
#       " (prior mean lambda = ",
#       round(as.numeric(cp$weighted_mean_rate), 6),
#       "; n_prior = ",
#       round(as.numeric(cp$n_prior_eff), 6),
#       ")\n\n",
#       sep = ""
#     )
#   }
# 
#   if (!is.null(x$conj_binomial)) {
#     cb <- x$conj_binomial
#     cat("Binomial probability prior (Beta conjugate; pass to dBeta()):\n")
#     cat(
#       "  shape1 = ",
#       round(as.numeric(cb$shape1), 6),
#       "; shape2 = ",
#       round(as.numeric(cb$shape2), 6),
#       " (prior mean theta = ",
#       round(as.numeric(cb$weighted_mean_prop), 6),
#       "; n_prior = ",
#       round(as.numeric(cb$n_prior_eff), 6),
#       ")\n\n",
#       sep = ""
#     )
#   }
# 
#   if (is.null(x$shape) && is.null(x$rate) && !is.null(x$dispersion)) {
#     cat("Note: Gaussian family detected, but shape/rate parameters were not computed.\n")
#     cat("This may occur if n_prior is not scalar.\n\n")
#   }
#   
#   invisible(x)
# }



#' Checks for Prior-data conflicts
#'
#' Checks if the credible intervals for the prior overlap with the implied confidence intervals
#' from the classical model (obtained via \code{\link[stats]{glm}}). The approach relates to
#' prior-data conflict checks \insertCite{EvansMoshonov2006}{glmbayes}.
#'
#' @param level the confidence level at which the Prior-data conflict should be checked.
#' @inheritParams glmb
#' @return A vector where each item provided the ratio of the absolue value for the difference between the 
#' prior and maximum likelihood estimate divided by the length of the sum of half of the two intervals 
#' (where normality is assumed)
#' @seealso \code{\link{Prior_Setup}}, \code{\link{glmb}}; see \insertCite{glmbayesChapter03}{glmbayes} for prior tailoring;
#' \insertCite{glmbayesChapterA12}{glmbayes} for full derivations.
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @family prior
#' @example inst/examples/Ex_Prior_Check.R
#' @export
#' @rdname Prior_Check
#' @order 1

Prior_Check<-function(formula,family,pfamily,level=0.95,data=NULL, weights, subset,na.action, 
                      start = NULL, etastart, mustart, offset ,control = list(...) , model = TRUE, 
                      method = "glm.fit",x = FALSE, y = TRUE, contrasts = NULL, ...){
  
  pf=pfamily
  prior_list=pfamily$prior_list
  
  ## For now, the below is really only correct for the dNormal pfamily
  
  mu=prior_list$mu
  Sigma=prior_list$Sigma
  
  
  object=glm(formula=formula,family=family,data=data)
  
  Like_est=object$coefficients
  Like_std=summary(object)$coefficients[,2]
  
  if(is.null(mu)){
    message("No Prior mean vector provided. Variables with needed Priors are:")
    message(paste(names(Like_est), collapse = ", "))
    names(Like_est)    
    
  }
  
  
  if(level<=0.5) stop("level must be greater than 0.5")
  
  Sigma=as.matrix(Sigma)
  Prior_std=sqrt(diag(Sigma))
  
  message("Variables in the Model Are: ", paste(names(Like_est), collapse = ", "))
  std_dev_sum=qnorm(level)*(Prior_std+Like_std)
  
  abs_ratio=matrix(rep(0,length(Like_est),nrow=length(Like_est),ncol=1))
  abs_diff=abs(mu-Like_est)
  abs_ratio[1:length(Like_est),]=abs_diff/std_dev_sum
  
  rownames(abs_ratio)=names(Like_est)
  colnames(abs_ratio)=c("abs_ratio")
  max_abs_ratio=max(abs_ratio)
  
  if(max_abs_ratio>1) {
    warning("At least one of the maximum likelihood estimates appears to be inconsistent with the prior")
  }
  
  else{
    message("The maximum likelihood estimates for all coefficients appear to be roughly consistent with the prior.")
  }
  return(abs_ratio)
  
}




