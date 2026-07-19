#' Prior Family Objects for Bayesian Models
#'
#' Prior family objects provide a convenient way to specify the details of the priors 
#' used by functions such as \code{\link{glmb}}. See the documentations for \code{\link{lmb}},
#' \code{\link{glmb}}, \code{\link{glmb}}, and \code{\link{rglmb}} for the details of how such model fitting 
#' takes place.
#' @name pfamily
#' @param object the function \code{pfamily} accesses the \code{pfamily} objects which
#' are stored within objects created by modelling functions (e.g., \code{glmb}).
#' @param mu a prior mean vector for the the modeling coefficients used in several pfamilies
#' @param Sigma a prior variance-covariance matrix for \code{dNormal()} and
#'   \code{dIndependent_Normal_Gamma()}.
#' @param Sigma_0 prior variance-covariance on the precision-weighted coefficient scale for
#'   \code{dNormal_Gamma()} only (Gaussian). Stored in \code{prior_list$Sigma} for compatibility
#'   with downstream samplers.
#' @param dispersion the dispersion to be assumed when it is not given a prior. Should be provided
#' when the Normal prior is for the \code{gaussian()}, \code{Gamma()}, \code{quasibinomial},
#' or \code{quasipoisson} families. The \code{binomial()} and \code{poisson()} families
#' do not have dispersion coefficients. Omitted or \code{NULL} uses the internal default
#' \code{1} and sets \code{ddef} in \code{prior_list} (see Details).
#' @param shape The prior shape parameter for the gamma piece (inverse dispersion / precision).
#'   When taking defaults from \code{\link{Prior_Setup}}, use \code{ps$shape} with
#'   \code{\link{dNormal_Gamma}()} and \code{\link{dGamma}()}, and \code{ps$shape_ING} with
#'   \code{\link{dIndependent_Normal_Gamma}()} on the Gaussian calibrated path (see Details).
#' @param rate The prior rate parameter paired with \code{shape}. With Gaussian
#'   \code{\link{Prior_Setup}}, \code{\link{dNormal_Gamma}()} and \code{\link{dIndependent_Normal_Gamma}()}
#'   use \code{ps$rate}; for \code{\link{dGamma}()} with fixed \code{beta}, prefer \code{ps$rate_gamma}
#'   when that field is non-\code{NULL} (see Details).
#' @param beta the regression coefficients to be assumed when it is not given a prior. 
#' Needs to be provided when the Gamma prior is used for the dispersion. This
#' specification is typically only used as part of Gibbs sampling where the beta and 
#' dispersion parameters are updated separately.
#' @param Inv_Dispersion Logical (default \code{TRUE}).  Controls which of the two Gamma prior
#'   roles \code{dGamma()} plays:
#'   \itemize{
#'     \item \code{TRUE} (default) — prior on the \emph{inverse dispersion} (precision / shape
#'       parameter \eqn{k = 1/\phi}).  This is the classical path used for dispersion estimation
#'       in Gaussian and \code{Gamma(log)} regression (\code{simfun = rGamma_reg}).
#'     \item \code{FALSE} — conjugate prior on the Gamma or Poisson \emph{rate} \eqn{\beta}
#'       directly (intercept-only, identity link).  The posterior is a closed-form Gamma draw
#'       (\code{simfun = rGamma_Conjugate_reg}).
#'   }
#' @param lik_shape Known shape parameter \eqn{k > 0} of the Gamma likelihood.  Only used when
#'   \code{Inv_Dispersion = FALSE} and \code{family = Gamma(link = "identity")}.  The intercept
#'   coefficient is then the Gamma \emph{rate} \eqn{\beta}, and the conjugate posterior is
#'   \eqn{\beta \mid y \sim \mathrm{Gamma}(\alpha_0 + n k,\; \beta_0 + \sum y_i)}.
#'   Defaults to \code{1} (exponential distribution). Ignored for Poisson families and whenever
#'   \code{Inv_Dispersion = TRUE}.
#' @param max_disp_perc Specifies the percentile used to truncate the posterior dispersion 
#' distribution when constructing the envelope for accept-reject sampling. This determines 
#' the lower and upper bounds for the dispersion (\eqn{\sigma^2}) used in the simulation. A value of 0.99
#' corresponds to using the central 98 percent of the posterior dispersion mass (i.e., excluding 
#' the outer 1 percent in each tail). Smaller values yield tighter bounds and may improve acceptance 
#' rates, while larger values allow broader dispersion support but may increase envelope complexity.
#' @param disp_lower lower bound truncation for dispersion 
#' @param disp_upper upper bound truncation for dispersion
#' @param x an object, a pfamily function that is to be printed
#' @param \ldots additional argument(s) for methods.
#' @details
#' \code{pfamily} is a generic with methods for fitted objects such as \code{\link{glmb}} and
#' \code{\link{lmb}}. The \code{dNormal()} prior is supported for all response families.
#' The \code{gaussian()} family additionally supports \code{dNormal_Gamma()},
#' \code{dIndependent_Normal_Gamma()}, and \code{dGamma()} (precision prior).
#' Intercept-only models with an identity link support two closed-form conjugate priors:
#' \code{dBeta()} for \code{binomial(link = "identity")} and
#' \code{dGamma(Inv_Dispersion = FALSE)} for \code{poisson(link = "identity")} and
#' \code{Gamma(link = "identity")}.
#'
#' A `pfamily` object represents a structured prior specification for use in Bayesian generalized linear modeling.
#' Each constructor function (e.g., `dNormal()`, `dGamma()`, `dNormal_Gamma()`, `dBeta()`) returns an object of
#' class `"pfamily"` containing the prior parameters, supported likelihood families, compatible link functions,
#' and a simulation function for posterior sampling.
#'
#' These priors are designed to integrate seamlessly with modeling functions such as `glmb()` and `rlmb()` in the 
#' \pkg{glmbayes} package, which consume the `pfamily` object to define the prior distribution over model parameters. 
#' The `pfamily()` generic retrieves the embedded prior from a fitted model object, while `print.pfamily()` displays its structure.
#'
#' **\code{prior_list} and \code{simfun}.** The named list \code{prior_list} holds the hyperparameters for the chosen
#' prior family. When a model function draws from the posterior, it passes \code{prior_list} into the element
#' \code{simfun} (e.g., \code{\link{rNormal_reg}}, \code{\link{rGamma_reg}}) so the low-level sampler receives
#' one consistent list structure regardless of which constructor built the \code{pfamily}.
#'
#' **\code{\link{Prior_Setup}} and default hyperparameters.** \code{Prior_Setup()} fits an auxiliary GLM and returns
#' default \code{mu}, \code{Sigma} / \code{Sigma_0}, \code{dispersion}, Gamma \code{shape} and \code{rate}, and
#' related fields aligned with the data and prior-weight (\code{pwt}) choices. Those values can be supplied as
#' arguments to the \code{pfamily} constructors when you want package-default priors on the same scale as the model
#' matrix. Recommended use of \code{shape} and \code{rate} is not identical across constructors: for
#' \code{\link{dIndependent_Normal_Gamma}()}, pass \code{shape = ps$shape_ING} from \code{Prior_Setup} (not the
#' scalar \code{ps$shape} used by \code{\link{dNormal_Gamma}()}). For \code{\link{dGamma}()} with fixed coefficients
#' (\code{beta}), pass \code{rate = ps$rate_gamma} when that field is present (otherwise \code{ps$rate}); see
#' \code{\link{Prior_Setup}} and \code{\link{compute_gaussian_prior}}.
#'
#' ## Prior Families
#'
#' - **`dNormal()`**: Specifies a multivariate normal prior over regression coefficients. It is conjugate for 
#'   Gaussian likelihoods with an identity link function, and serves as the primary implemented prior for all 
#'   other supported likelihood families in the current framework. This structure facilitates efficient posterior 
#'   sampling and analytical tractability. The returned \code{prior_list} includes \code{ddef}: \code{TRUE} when
#'   \code{dispersion} was omitted or \code{NULL} (so the default \code{1} was used), \code{FALSE} when
#'   \code{dispersion} was supplied explicitly (including \code{1}).
#'
#'   For models with log-concave likelihood functions-such as Poisson, Binomial, and Gamma families-
#'   posterior sampling under a `dNormal` prior is performed using a \insertCite{Nygren2006}{glmbayes} 
#'   likelihood subgradient approach. This method constructs tight enveloping functions around the posterior 
#'   using subgradients of the log-likelihood, enabling efficient accept-reject sampling even in high dimensions.
#'
#'   When the posterior distribution is approximately normal (typically the case for large sample sizes), the 
#'   area under the enveloping function is bounded above by a constant factor-approximately \eqn{2 / \sqrt{\pi} \approx 1.128} 
#'   in the univariate case, and \eqn{(2 / \sqrt{\pi})^k} in \eqn{k}-dimensional models. These bounds ensure that 
#'   the rejection rate remains manageable and that the sampler remains computationally efficient.
#'
#'   The concept of conjugate priors was first formalized by \insertCite{Raiffa1961}{glmbayes}, and further 
#'   developed for regression models using g-prior structures by \insertCite{zellner1986gprior}{glmbayes}.
#'
#' - **`dGamma()`**: A Gamma prior with two distinct roles controlled by \code{Inv_Dispersion}:
#'   \itemize{
#'     \item \code{Inv_Dispersion = TRUE} (default): prior on the inverse dispersion (precision
#'       \eqn{1/\phi} or shape \eqn{k}). Used for dispersion estimation in Gaussian and
#'       Gamma(log) models, typically in a Gibbs step with \code{beta} held fixed
#'       \insertCite{Gelman2013,Dobson1990,McCullagh1989}{glmbayes}.
#'       With Gaussian \code{\link{Prior_Setup}} output, prefer \code{rate_gamma} for \code{rate}
#'       (see Details above).
#'     \item \code{Inv_Dispersion = FALSE}: conjugate Gamma prior on the rate parameter
#'       \eqn{\beta} directly. Supports intercept-only models with an identity link:
#'       Poisson (Gamma–Poisson conjugacy) and Gamma (Gamma–Gamma conjugacy).
#'       Posterior draws are closed-form IID samples via \code{\link{rGamma_Conjugate_reg}}.
#'       The \code{lik_shape} argument specifies the known Gamma likelihood shape (default 1,
#'       i.e.\ exponential). \code{\link{Prior_Setup}} returns calibrated \code{conj_poisson}
#'       hyperparameters for this path.
#'   }
#'
#' - **`dBeta()`**: A Beta prior on the binomial probability \eqn{\theta} for intercept-only
#'   \code{binomial(link = "identity")} models. The posterior is a closed-form Beta draw
#'   (Beta–Binomial conjugacy) produced by \code{\link{rBeta_reg}}. Arguments \code{shape1}
#'   and \code{shape2} are the prior pseudo-success and pseudo-failure counts.
#'   \code{\link{Prior_Setup}} returns calibrated \code{conj_beta} hyperparameters for this path.
#'
#' - **`dNormal_Gamma()`**: Combines a multivariate normal prior on coefficients with a gamma prior on precision,
#'   forming a conjugate structure for Gaussian models with unknown variance. The second argument is \code{Sigma_0}
#'   (precision-weighted scale); it is aliased internally to \code{Sigma} in \code{prior_list}.
#'   This formulation parallels classical Normal-Gamma models and is compatible with hierarchical extensions
#'   \insertCite{Gelman2013,Raiffa1961}{glmbayes}.
#'
#' - **`dIndependent_Normal_Gamma()`**: Similar to `dNormal_Gamma()`, but assumes independence between the
#'   coefficient and precision priors. This structure is useful for models where prior independence is desired
#'   or analytically convenient. With \code{\link{Prior_Setup}} on a Gaussian model, pass \code{shape_ING} as
#'   the \code{shape} argument (see Details above).
#'
#' Each `pfamily` object includes:
#' - `pfamily`, `prior_list`, `okfamilies`, `plinks`, `simfun`, and `pfun` (see Value).
#'
#' @return An object of class \code{"pfamily"} (with a concise \code{print} method). A list with elements:
#' \item{pfamily}{Character string: the constructor name (\code{"dNormal"}, \code{"dGamma"},
#'   \code{"dNormal_Gamma"}, \code{"dIndependent_Normal_Gamma"}, or \code{"dBeta"}).}
#' \item{prior_list}{Named list of prior hyperparameters. It is passed into \code{simfun} when sampling so the
#'   relevant low-level routine receives the prior in a fixed list form. Contents depend on the constructor:
#'   \describe{
#'     \item{\code{dNormal}:}{\code{mu}, \code{Sigma}, \code{dispersion}, and logical \code{ddef}
#'       (\code{TRUE} if \code{dispersion} was omitted or \code{NULL}, so the default \code{1} was used;
#'       \code{FALSE} if set explicitly).}
#'     \item{\code{dGamma}:}{\code{shape}, \code{rate}, \code{beta}, \code{Inv_Dispersion},
#'       \code{max_disp_perc}, \code{disp_lower}, \code{disp_upper}. When \code{Inv_Dispersion = FALSE},
#'       also includes surrogate \code{mu} and \code{Sigma} (computed from the Gamma prior moments)
#'       and \code{lik_shape}.}
#'     \item{\code{dNormal_Gamma}:}{\code{mu}, \code{Sigma} (the \code{Sigma_0} precision-weighted input),
#'       \code{shape}, \code{rate}.}
#'     \item{\code{dIndependent_Normal_Gamma}:}{\code{mu}, \code{Sigma} (coefficient-scale covariance),
#'       \code{shape}, \code{rate}, \code{max_disp_perc}, \code{disp_lower}, \code{disp_upper}.}
#'     \item{\code{dBeta}:}{\code{shape1}, \code{shape2}, \code{beta}, and surrogate \code{mu} and
#'       \code{Sigma} computed from the Beta prior moments
#'       (\code{mu = shape1/(shape1+shape2)},
#'        \code{Sigma = shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1))}).}
#'   }
#' }
#' \item{okfamilies}{Character vector of implemented \code{\link[stats]{family}} names for which this
#'   \code{pfamily} may be used.}
#' \item{plinks}{Function of one \code{family} argument returning allowed link names for that family.}
#' \item{simfun}{Function used to generate posterior draws (e.g., \code{\link{rNormal_reg}},
#'   \code{\link{rGamma_reg}}, \code{\link{rGamma_Conjugate_reg}}, \code{\link{rNormalGamma_reg}}, \code{\link{rindepNormalGamma_reg}});
#'   for standard use these produce i.i.d.\ posterior samples for the implemented settings.}
#' \item{pfun}{Prior-simulation function paired with this constructor (e.g.,
#'   \code{\link{rNormal_prior}}, \code{\link{rGamma_prior}},
#'   \code{\link{rNormal_Gamma_prior}}); called by
#'   \code{\link[bayestestR]{simulate_prior}} methods on fitted objects via
#'   \code{pfun(n, prior_list, params)}.}
#' 
#' @author The design of the \code{pfamily} set of functions was developed by Kjell Nygren and was 
#' inspired by the family used by the \code{\link{glmb}} function to specify the likelihood 
#' function. That design in turn was inspired by S functions of the same names from
#' the statistical modeling literature.
#'
#' @seealso
#' \code{\link{glmb}}, \code{\link{rlmb}}, \code{\link{lmb}}, \code{\link{rglmb}} for modeling functions that consume \code{pfamily} objects.
#'
#' \code{\link{rNormal_reg}}, \code{\link{rNormalGamma_reg}}, \code{\link{rGamma_reg}}, \code{\link{rGamma_Conjugate_reg}}, \code{\link{rindepNormalGamma_reg}} for lower-level sampling functions used by \code{pfamily} constructors.
#'
#' \code{\link{Prior_Setup}}, \code{\link{Prior_Check}} for initializing and validating prior specifications.
#'
#' \code{\link{EnvelopeBuild}} for envelope construction methods used in likelihood subgradient sampling \insertCite{Nygren2006}{glmbayes}.
#'
#' See also \insertCite{Hastie1992}{glmbayes} for the original S modeling framework that inspired the design of \code{pfamily}.
#'
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#'
#' @example inst/examples/Ex_pfamily.R
#' @export 
#' @rdname pfamily
#' @order 1

pfamily <- function(object, ...) UseMethod("pfamily")

#' @export 
#' @method pfamily default

pfamily.default <- function(object, ...){

  if(is.null(object$pfamily)) stop("no pfamily object found")
  if (!inherits(object$pfamily, "pfamily"))  stop("Object named pfamily is not of class pfamily")
  
  return(object$pfamily)
}


#' @export
#' @method print pfamily
#' @rdname pfamily
#' @order 7

print.pfamily <- function(x, ...)
{
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Prior Family:", x$pfamily, "\n\n")
  cat("Prior List:\n\n")
  print(x$prior_list)
  
  invisible(x)
}

#' @export 
#' @rdname pfamily
#' @order 2

dNormal<-function(mu,Sigma,dispersion=NULL){
  
  ## Check that the inputs are numeric
  
  if(is.numeric(mu)==FALSE||is.numeric(Sigma)==FALSE) stop("non-numeric argument to numeric function")

  mu=as.matrix(mu,ncol=1) ## Force mu to matrix
  Sigma=as.matrix(Sigma)  ## Force Sigma to matrix 
  
  nvar=length(mu)
  nvar1=nrow(Sigma)
  nvar2=ncol(Sigma)
  
  if(!nvar==nvar1||!nvar==nvar2) stop("dimensions of mu and Sigma are not consistent")

  ## Check for symmetry and positive definiteness
  if(!isSymmetric(Sigma))stop("matrix Sigma must be symmetric")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(Sigma, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  thr <- -tol * abs(ev[1])   # = -1e-06 * 12.56941 ~= -1.256941e-05  
  
  
  
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  
  ddef <- missing(dispersion) || is.null(dispersion)
  if (ddef) dispersion <- 1
  if(!is.null(dispersion)){
    if(!is.numeric(dispersion)) stop("non-numeric argument to numeric function")
    if(!length(dispersion)==1) stop("dispersion has length>1")
    if(!length(dispersion)>0) stop("dispersion must be >0")
  }
    
  okfamilies <- c("gaussian","poisson","binomial","quasipoisson","quasibinomial","Gamma")

  plinks<-function(family){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-c("log")		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-c("logit","probit","cloglog")		
    if(family$family=="Gamma") oklinks<-c("log")	
    return(oklinks)
  }
  
  prior_list=list(mu=mu,Sigma=Sigma,dispersion=dispersion,ddef=ddef)
  attr(prior_list,"Prior Type")="dNormal"  

  outlist=list(pfamily="dNormal",prior_list=prior_list,okfamilies=okfamilies,
  plinks=plinks,             
  simfun=rNormal_reg,
  pfun=rNormal_prior)
  attr(outlist,"Prior Type")="dNormal"             
  class(outlist)="pfamily"
  outlist$call<-match.call()
  return(outlist)
  }

#' @export 
#' @rdname pfamily
#' @order 3

dGamma <- function(shape, rate, beta,
                   Inv_Dispersion = TRUE,
                   lik_shape      = 1,
                   max_disp_perc  = 0.99,
                   disp_lower     = NULL,
                   disp_upper     = NULL) {

  if (!is.numeric(shape) || !is.numeric(rate) || !is.numeric(beta))
    stop("non-numeric argument to numeric function")
  if (length(shape) > 1) stop("shape is not of length 1")
  if (length(rate)  > 1) stop("rate is not of length 1")
  if (shape <= 0) stop("shape must be > 0")
  if (rate  <= 0) stop("rate must be > 0")

  beta <- as.matrix(beta, ncol = 1L)

  ## -------------------------------------------------------------------------
  ## Inv_Dispersion = TRUE  →  prior on precision/shape (inverse dispersion).
  ## Supports Gaussian(identity) and Gamma(log); uses rGamma_reg sampler.
  ## -------------------------------------------------------------------------
  if (Inv_Dispersion) {

    okfamilies <- c("gaussian", "Gamma")

    plinks <- function(family) {
      if (family$family == "gaussian")                          oklinks <- c("identity")
      if (family$family %in% c("poisson", "quasipoisson"))     oklinks <- NULL
      if (family$family %in% c("binomial", "quasibinomial"))   oklinks <- NULL
      if (family$family == "Gamma")                            oklinks <- c("log")
      return(oklinks)
    }

    prior_list <- list(
      shape          = shape,
      rate           = rate,
      beta           = beta,
      Inv_Dispersion = TRUE,
      max_disp_perc  = max_disp_perc,
      disp_lower     = disp_lower,
      disp_upper     = disp_upper
    )
    attr(prior_list, "Prior Type") <- "dGamma"
    outlist <- list(pfamily    = "dGamma",
                    prior_list = prior_list,
                    okfamilies = okfamilies,
                    plinks     = plinks,
                    simfun     = rGamma_reg,
                    pfun       = rGamma_prior)
    attr(outlist, "Prior Type") <- "dGamma"

  ## -------------------------------------------------------------------------
  ## Inv_Dispersion = FALSE  →  conjugate prior on the rate β directly.
  ## Supports Poisson(identity) and Gamma(identity); uses rGamma_Conjugate_reg.
  ## mu / Sigma: Gamma(shape, rate) moments; mean = shape/rate, var = shape/rate^2.
  ## -------------------------------------------------------------------------
  } else {

    if (!is.numeric(lik_shape) || length(lik_shape) != 1L ||
        !is.finite(lik_shape) || lik_shape <= 0)
      stop("lik_shape must be a single positive finite number (the known Gamma likelihood shape parameter; default 1 for exponential)")

    sh <- as.numeric(shape)[[1L]]
    rt <- as.numeric(rate)[[1L]]
    mu <- beta * 0 + sh / rt
    p  <- nrow(mu)
    sigma_sq <- sh / (rt * rt)
    Sigma <- diag(rep.int(sigma_sq, times = p), nrow = p, ncol = p)
    coef_nm <- rownames(beta)
    if (is.null(coef_nm)) coef_nm <- colnames(beta)
    if (!is.null(coef_nm) && length(coef_nm) == p) {
      rownames(mu) <- coef_nm
      if (!is.null(colnames(beta))) colnames(mu) <- colnames(beta)
      dimnames(Sigma) <- list(coef_nm, coef_nm)
    }

    okfamilies <- c("poisson", "Gamma")

    plinks <- function(family) {
      oklinks <- NULL
      if (family$family %in% c("poisson", "quasipoisson")) oklinks <- c("identity")
      if (family$family == "Gamma")                        oklinks <- c("identity")
      oklinks
    }

    prior_list <- list(
      shape          = shape,
      rate           = rate,
      beta           = beta,
      lik_shape      = lik_shape,
      Inv_Dispersion = FALSE,
      mu             = mu,
      Sigma          = Sigma,
      max_disp_perc  = max_disp_perc,
      disp_lower     = disp_lower,
      disp_upper     = disp_upper
    )
    attr(prior_list, "Prior Type") <- "dGamma"
    outlist <- list(pfamily    = "dGamma",
                    prior_list = prior_list,
                    okfamilies = okfamilies,
                    plinks     = plinks,
                    simfun     = rGamma_Conjugate_reg,
                    pfun       = rGamma_Conjugate_prior)
    attr(outlist, "Prior Type") <- "dGamma"
  }

  class(outlist) <- "pfamily"
  outlist$call   <- match.call()
  return(outlist)
}


## dGamma_Conjugate() removed 2026-05-27 — functionality merged into dGamma(Inv_Dispersion = FALSE).
## Commented out rather than deleted to preserve the implementation history.
#
# #' @description
# #' \code{dGamma_Conjugate()} was a deprecated alias for \code{dGamma(..., Inv_Dispersion = FALSE)}.
# #' Use \code{dGamma(Inv_Dispersion = FALSE)} directly.
# #'
# #' @export
# #' @rdname pfamily
# #' @order 4
#
# dGamma_Conjugate <- function(shape, rate, beta, lik_shape = 1,
#                               max_disp_perc = 0.99,
#                               disp_lower    = NULL,
#                               disp_upper    = NULL) {
#   .Deprecated(
#     new = "dGamma",
#     msg = paste0(
#       "dGamma_Conjugate() is deprecated.\n",
#       "Use dGamma(..., Inv_Dispersion = FALSE) instead."
#     )
#   )
#   dGamma(shape = shape, rate = rate, beta = beta,
#          Inv_Dispersion = FALSE, lik_shape = lik_shape,
#          max_disp_perc  = max_disp_perc,
#          disp_lower     = disp_lower,
#          disp_upper     = disp_upper)
# }



#' Conjugate Beta prior family (\code{dBeta}: closed-form IID updates for intercept-only
#' Binomial models with an identity link).
#'
#' Under a Beta(\code{shape1}, \code{shape2}) prior on the binomial probability \eqn{\theta}
#' and a Binomial(\eqn{n_i}, \eqn{\theta}) likelihood with identity link (\eqn{\theta = \beta}
#' directly), the posterior is:
#' \deqn{\theta \mid y \sim \mathrm{Beta}(\texttt{shape1} + \sum n_i y_i,\;
#'   \texttt{shape2} + \sum n_i (1 - y_i)).}
#'
#' \code{mu} / \code{Sigma}: the surrogate Normal mean is \code{shape1/(shape1+shape2)} and
#' the surrogate variance is the Beta variance
#' \code{shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1))}.
#'
#' @param shape1 First shape parameter \eqn{\alpha > 0} of the Beta prior (prior successes + 1).
#' @param shape2 Second shape parameter \eqn{\beta > 0} of the Beta prior (prior failures + 1).
#' @param beta Initial coefficient matrix (1 \eqn{\times} 1); typically set to the prior mean
#'   \code{shape1/(shape1+shape2)}.
#' @export
#' @rdname pfamily
#' @order 5

dBeta <- function(shape1, shape2, beta) {

  if (!is.numeric(shape1) || !is.numeric(shape2) || !is.numeric(beta))
    stop("non-numeric argument to numeric function")
  if (length(shape1) != 1L) stop("shape1 must be a single positive number")
  if (length(shape2) != 1L) stop("shape2 must be a single positive number")
  if (!is.finite(shape1) || shape1 <= 0) stop("shape1 must be a finite positive number")
  if (!is.finite(shape2) || shape2 <= 0) stop("shape2 must be a finite positive number")

  beta <- as.matrix(beta, ncol = 1L)

  ## Normal-style surrogate for glmb() pre-simulation and downstream Prior$mean/Variance.
  ## Beta(shape1, shape2): mean = shape1/(shape1+shape2),
  ##   variance = shape1*shape2 / ((shape1+shape2)^2 * (shape1+shape2+1)).
  s1  <- as.numeric(shape1)[[1L]]
  s2  <- as.numeric(shape2)[[1L]]
  s12 <- s1 + s2
  prior_mean_val <- s1 / s12
  prior_var_val  <- s1 * s2 / (s12^2 * (s12 + 1))

  p      <- nrow(as.matrix(beta, ncol = 1L))
  mu     <- beta * 0 + prior_mean_val
  Sigma  <- diag(rep.int(prior_var_val, times = p), nrow = p, ncol = p)

  coef_nm <- rownames(beta)
  if (is.null(coef_nm)) coef_nm <- colnames(beta)
  if (!is.null(coef_nm) && length(coef_nm) == p) {
    rownames(mu) <- coef_nm
    if (!is.null(colnames(beta))) colnames(mu) <- colnames(beta)
    dimnames(Sigma) <- list(coef_nm, coef_nm)
  }

  okfamilies <- c("binomial", "quasibinomial")

  plinks <- function(family) {
    oklinks <- NULL
    if (family$family %in% c("binomial", "quasibinomial")) oklinks <- c("identity")
    oklinks
  }

  prior_list <- list(
    shape1 = shape1,
    shape2 = shape2,
    beta   = beta,
    mu     = mu,
    Sigma  = Sigma
  )
  attr(prior_list, "Prior Type") <- "dBeta"

  outlist <- list(
    pfamily    = "dBeta",
    prior_list = prior_list,
    okfamilies = okfamilies,
    plinks     = plinks,
    simfun     = rBeta_reg,
    pfun       = rBeta_prior
  )
  attr(outlist, "Prior Type") <- "dBeta"
  class(outlist) <- "pfamily"
  outlist$call   <- match.call()

  return(outlist)
}


#' @export
#' @rdname pfamily
#' @order 7

dNormal_Gamma <- function(mu, Sigma_0, shape, rate) {
  Sigma <- Sigma_0

  ############################################################  
  
  if(is.numeric(mu)==FALSE||is.numeric(Sigma)==FALSE) stop("non-numeric argument to numeric function")
  if(is.numeric(shape)==FALSE||is.numeric(rate)==FALSE) stop("non-numeric argument to numeric function")
  
  if(length(shape)>1) stop("shape is not of length 1")
  if(length(rate)>1) stop("rate is not of length 1")
  if(shape<=0) stop("shape must be>0")
  if(rate<=0) stop("rate must be>0")
  
  mu=as.matrix(mu,ncol=1) ## Force mu to matrix
  Sigma=as.matrix(Sigma)  ## Force Sigma to matrix 
    
  nvar=length(mu)
  nvar1=nrow(Sigma)
  nvar2=ncol(Sigma)
  
  if(!nvar==nvar1||!nvar==nvar2) stop("dimensions of mu and Sigma are not consistent")
  
  ## Check for symmetry and positive definiteness
  if(!isSymmetric(Sigma))stop("matrix Sigma must be symmetric")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(Sigma, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  
  
  ############################################################
  
  okfamilies <- c("gaussian") # Unclear if this could be used for Gamma  or quasi-families

  plinks<-function(family){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-NULL		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-NULL		
    if(family$family=="Gamma") oklinks<-NULL	
    return(oklinks)
  }
  
  prior_list=list(mu=mu,Sigma=Sigma,shape=shape,rate=rate)
  attr(prior_list,"Prior Type")="dNormal_Gamma"  
  outlist=list(pfamily="dNormal_Gamma",call=call,prior_list=prior_list,
    okfamilies=okfamilies,plinks=plinks,simfun=rNormalGamma_reg,
    pfun=rNormal_Gamma_prior)
  
  attr(outlist,"Prior Type")="dNormal_Gamma"             
  class(outlist)="pfamily"
  outlist$call<-match.call()
  
  return(outlist)
  }



#' @export 
#' @rdname pfamily
#' @order 8

dIndependent_Normal_Gamma <- function(mu, Sigma, shape, rate, max_disp_perc = 0.99,disp_lower=NULL,disp_upper=NULL) {

  ##############################################################
  
  if(is.numeric(mu)==FALSE||is.numeric(Sigma)==FALSE) stop("non-numeric argument to numeric function")
  if(is.numeric(shape)==FALSE||is.numeric(rate)==FALSE) stop("non-numeric argument to numeric function")
  
  if(length(shape)>1) stop("shape is not of length 1")
  if(length(rate)>1) stop("rate is not of length 1")
  if(shape<=0) stop("shape must be>0")
  if(rate<=0) stop("rate must be>0")
  if (!is.numeric(max_disp_perc) || length(max_disp_perc) != 1 || max_disp_perc <= 0.5 || max_disp_perc >= 1) {
    stop("max_disp_perc must be a single number between 0.5 and 1")
  }
  
  mu=as.matrix(mu,ncol=1) ## Force mu to matrix
  Sigma=as.matrix(Sigma)  ## Force Sigma to matrix 
  

  nvar=length(mu)
  nvar1=nrow(Sigma)
  nvar2=ncol(Sigma)
  
  if(!nvar==nvar1||!nvar==nvar2) stop("dimensions of mu and Sigma are not consistent")
  
  ## Check for symmetry and positive definiteness
  if(!isSymmetric(Sigma))stop("matrix Sigma must be symmetric")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(Sigma, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  
  
  ##############################################################
  
  okfamilies <- c("gaussian") # Unclear if this could be used for Gamma or quasi-families
  
  plinks<-function(family){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-NULL		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-NULL
    if(family$family=="Gamma") oklinks<-NULL	
    return(oklinks)
  }
  
  
  prior_list <- list(mu = mu, Sigma = Sigma, shape = shape, rate = rate, max_disp_perc = max_disp_perc,
                     disp_lower=disp_lower,disp_upper=disp_upper)
  attr(prior_list,"Prior Type")="dIndependent_Normal_Gamma"  
  outlist=list(pfamily="dIndependent_Normal_Gamma",prior_list=prior_list,
               okfamilies=okfamilies,plinks=plinks,simfun=rindepNormalGamma_reg,
               pfun=rIndependent_Normal_Gamma_prior)
  
  attr(outlist,"Prior Type")="dIndependent_Normal_Gamma"             
  class(outlist)="pfamily"
  outlist$call<-match.call()
  
  return(outlist)
  
}


