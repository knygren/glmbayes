#' The Bayesian Generalized Linear Model Distribution
#'
#' \code{rglmb} is used to generate iid samples for Bayesian Generalized Linear Models.
#' The model is specified by providing a data vector, a design matrix, 
#' the family (determining the likelihood function) and the pfamily (determining the 
#' prior distribution).
#' @name rglmb
#' @param y a vector of observations of length \code{m}.
#' @param x for \code{rglmb} a design matrix of dimension \code{m * p} and for \code{print.rglmb} the object to be printed. 
#' @inheritParams glmb
#' @param n_envopt Effective sample size passed to EnvelopeOpt for grid
#'   construction. Defaults to match `n`. Larger values encourage tighter
#'   envelopes.
#' @param use_parallel Logical. Whether to use parallel processing during simulation.
#' @param use_opencl Logical. Whether to use OpenCL acceleration during Envelope construction.
#' @param verbose Logical. Whether to print progress messages.
#' @return \code{rglmb} returns a object of class \code{"rglmb"}.  The function \code{summary} 
#' (i.e., \code{\link{summary.rglmb}}) can be used to obtain or print a summary of the results.
#' The generic accessor functions \code{\link{coefficients}}, \code{\link{fitted.values}},
#' \code{\link{residuals}}, and \code{\link{extractAIC}} can be used to extract
#' various useful features of the value returned by \code{\link{rglmb}}.
#' An object of class \code{"rglmb"} is a list containing at least the following components:
#' \item{coefficients}{a matrix of dimension \code{n} by \code{length(mu)} with one sample in each row}
#' \item{coef.mode}{a vector of \code{length(mu)} with the estimated posterior mode coefficients}
#' \item{dispersion}{Either a constant provided as part of the call, or a vector of length \code{n} with one sample in each row.}
#' \item{Prior}{A list with the priors specified for the model in question. Items in the
#' list may vary based on the type of prior}
#' \item{prior.weights}{a vector of weights specified or implied by the model} 
#' \item{y}{a vector with the dependent variable} 
#' \item{x}{a matrix with the implied design matrix for the model} 
#' \item{famfunc}{Family functions used during estimation process}
#' \item{iters}{an \code{n} by \code{1} matrix giving the number of candidates generated before acceptance for each sample.}
#' \item{Envelope}{the envelope that was used during sampling}
#' 
#' Objects of class \code{"rglmb"} are normally of class \code{c("rglmb","glmb","glm","lm")},
#' meaning they inherit from \code{glmb}, \code{glm}, and \code{lm}. This allows methods defined
#' for these upstream classes to be applied to \code{"rglmb"} objects when appropriate, while
#' supporting extensions for regularized Bayesian GLMs with structured priors.
#' 
#' @details
#' The function \code{rglmb} is a minimalistic engine for Bayesian generalized linear model simulation. 
#' It is designed to generate independent draws from the posterior distribution of a GLM, given a design matrix, 
#' response vector, likelihood family, and prior specification. Unlike \code{\link{glmb}}, which wraps formula parsing, 
#' model setup, and method dispatch, \code{rglmb} operates directly on numeric inputs and is optimized for speed, 
#' transparency, and integration into simulation workflows.
#'
#' The original R implementation of \code{glm} was written by Simon Davies (under Ross Ihaka at the University of Auckland) 
#' and has since been extensively rewritten by members of the R Core Team; its design was inspired by the S function 
#' described in \insertCite{Hastie1992}{glmbayes}, which in turn relies on the formula framework described in 
#' \insertCite{WilkinsonRogers1973}{glmbayes}.
#'
#' The design of the \code{pfamily} family of functions was created by Kjell Nygren and is modeled on how 
#' \code{glm} uses \code{family} to specify the likelihood. For any implemented combination of family, link, and 
#' \code{pfamily}, \code{rglmb} generates independent draws from the posterior density-no MCMC chains are required.
#'
#' A helper, \code{\link{Prior_Setup}}, assists users in choosing prior parameters. It ships with sensible defaults 
#' but also allows full customization. In particular, the default for \code{dNormal} is a reparameterization of 
#' Zellner's g-prior \insertCite{zellner1986gprior}{glmbayes}.
#'
#' Currently supported response families are \code{gaussian} (identity link), \code{poisson} and \code{quasipoisson} 
#' (log link), \code{gamma} (log link), and \code{binomial} and \code{quasibinomial} (logit, probit, cloglog). 
#' All families support a \code{dNormal} prior; the Gaussian family also offers \code{dNormalGamma} and 
#' \code{dIndependent_Normal_Gamma}.
#'
#' For the Gaussian family, draws under \code{dNormal} and \code{dNormalGamma} come from posterior distributions 
#' resulting from conjugate prior distributions \insertCite{Raiffa1961}{glmbayes}. For all other priors or response families, 
#' we use an accept-reject sampler built on the likelihood-subgradient envelope method 
#' \insertCite{Nygren2006}{glmbayes}. The \code{Gridtype} argument controls how many tangent points are used 
#' in the envelope-trading off envelope tightness against construction cost-and \code{iters} reports candidate 
#' counts before acceptance.
#'
#' By default, \code{rglmb} draws \code{n = 1} sample, uses parallel CPU simulation, and-if \code{use_opencl = TRUE}-
#' GPU-accelerated envelope building. The function returns a list containing posterior samples, prior specifications, 
#' dispersion estimates, and the envelope used during sampling. It does not return a full model object, and does not 
#' support formula-based modeling or method dispatch. Instead, it is called internally by \code{\link{glmb}} and 
#'and may be useful for Gibbs sampling implementations or other workflows where full model 
#' reconstruction is unnecessary.
#'  
#' @author The \R implementation of \code{rglmb} has been written by Kjell Nygren and
#' was built to be a Bayesian version of the \code{glm} function but with a more minimalistic interface 
#' than the \code{glmb} function. It also borrows some of its structure from other random generating function 
#' like \code{\link{rnorm}} and hence the \code{r} prefix. 
#' 
#' @family modelfuns
#' @seealso \code{\link{glmb}} for the formula interface; \code{\link[stats]{lm}} and
#' \code{\link[stats]{glm}} for classical modeling functions.
#'
#' \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeSize}}, \code{\link{EnvelopeEval}}
#' for envelope construction and grid evaluation used in non-conjugate sampling.
#' 
#' \code{\link{family}} for documentation of family functions used to specify priors.
#' \code{\link{pfamily}} for documentation of pfamily functions used to specify priors.
#' 
#' \code{\link{Prior_Setup}}, \code{\link{Prior_Check}} for functions used to initialize and to check priors,  
#'
#' Further reading: \insertCite{Nygren2006}{glmbayes};
#' \insertCite{glmbayesSimmethods,glmbayesChapterA08}{glmbayes};
#' OpenCL/GPU: \insertCite{glmbayesChapter12,glmbayesChapterA10}{glmbayes}.
#'
#' \code{\link{summary.glmb}}, \code{\link{predict.glmb}}, \code{\link{residuals.glmb}}, \code{\link{simulate.glmb}}, 
#' \code{\link{extractAIC.glmb}}, \code{\link{dummy.coef.glmb}} and methods(class="glmb") for \code{glmb} 
#' and the methods and generic functions for classes \code{glm} and \code{lm} from which class \code{glmb} inherits.
#' 
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_rglmb.R
#' @export 
#' @rdname rglmb
#' @order 1



rglmb<-function(n=1,y,x,family=gaussian(),pfamily,offset=NULL,
                weights=1,
                Gridtype=2,
                n_envopt = NULL,          # NEW
                use_parallel = TRUE, 
                use_opencl = FALSE, 
                verbose = FALSE){
  

  
  ## normalize n_envopt
  if (is.null(n_envopt)) n_envopt <- n
  n_envopt <- as.integer(n_envopt)
  
  ## Resolve family before pfamily checks (needed for Poisson + dGamma coercion)
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  ## family = Poisson + prior from dGamma(Inv_Dispersion=TRUE) -> auto-coerce to conjugate rate prior
  if (inherits(pfamily, "pfamily") && identical(pfamily$pfamily, "dGamma") &&
      identical(family$family, "poisson")) {
    pl <- pfamily$prior_list
    pfamily <- dGamma(
      shape          = pl$shape,
      rate           = pl$rate,
      beta           = pl$beta,
      Inv_Dispersion = FALSE,
      max_disp_perc  = pl$max_disp_perc,
      disp_lower     = pl$disp_lower,
      disp_upper     = pl$disp_upper
    )
  }
  
  
  
  ## Pull in information from the pfamily  
  pf=pfamily$pfamily
  okfamilies=pfamily$okfamilies  
  plinks=pfamily$plinks
  prior_list=pfamily$prior_list 
  simfun=pfamily$simfun
  
  ## dGamma(Inv_Dispersion=FALSE) + Gamma / Poisson: scalar conjugate requires intercept-only designs
  if (isFALSE(prior_list$Inv_Dispersion) &&
      family$family %in% c("Gamma", "poisson", "quasipoisson")) {
    x_ck <- if (is.matrix(x)) x else as.matrix(x)
    n_ck <- nrow(x_ck)
    if (length(y) != n_ck) {
      stop("`rglmb()`: length(y) differs from nrow(x); cannot validate intercept-only conjugate layout.", call. = FALSE)
    }
    alp <- offset
    if (is.null(alp)) {
      alp <- rep(0, n_ck)
    } else if (!is.numeric(alp)) {
      stop("`rglmb()`: offset must be numeric for dGamma conjugate models.", call. = FALSE)
    } else if (length(alp) == 1L) {
      alp <- rep(alp, n_ck)
    } else if (length(alp) != n_ck) {
      stop("`rglmb()`: offset must be scalar or length nrow(x) when using dGamma conjugate.", call. = FALSE)
    }
    wg <- weights
    if (!is.numeric(wg))
      stop("`rglmb()`: weights must be numeric.", call. = FALSE)
    if (length(wg) == 1L) wg <- rep(wg, n_ck)
    if (length(wg) != n_ck) {
      stop("`rglmb()`: weights must be scalar or length nrow(x) when using dGamma conjugate.", call. = FALSE)
    }
    .check_gamma_conjugate_scalar_design(
      x_ck,
      prior_list$beta,
      alp,
      wg,
      family_label = family$family,
      ctx = "`rglmb()`"
    )
  }

  ## Pull in information on families  
  # if (is.character(family)) 
  #   family <- get(family, mode = "function", envir = parent.frame())
  # if (is.function(family)) 
  #   family <- family()
  # if (is.null(family$family)) {
  #   print(family)
  #   stop("'family' not recognized")
  # }
  
  ## Check that the family is implemented for the pfamily
  ##
  ## dGamma(Inv_Dispersion=FALSE) keeps a narrow `okfamilies` while `plinks()` still
  ## advertises conjugate-compatible links for Poisson / quasi-Poisson; allow those families here.
  family_ok <- family$family %in% okfamilies
  if (isFALSE(prior_list$Inv_Dispersion) &&
      family$family %in% c("poisson", "quasipoisson")) {
    family_ok <- TRUE
  }

  if(family_ok){
    oklinks=plinks(family)
    if(!family$link %in% oklinks){      
      stop(gettextf("link \"%s\" not available for selected pfamily/family combination; available links are %s", 
                    family$link , paste(sQuote(oklinks), collapse = ", ")), domain = NA)
    }
  }
  else{
    stop(gettextf("family \"%s\" not available for current pfamily; available families are %s", 
                  family$family , paste(sQuote(okfamilies), collapse = ", ")), 
         domain = NA)
    
  }
  
  
  simfun_args <- list(
    n = n,
    y = y,
    x = x,
    prior_list = prior_list,
    offset = offset,
    weights = weights,
    family = family,
    Gridtype = Gridtype,
    n_envopt = n_envopt,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = verbose
  )
  
  ## Call relevant simulation function (for now without control2 list)
  outlist = simfun(n = n, y = y, x = x, prior_list = prior_list,offset = offset, weights = weights, family = family, 
                   Gridtype=Gridtype,n_envopt = n_envopt,   # pass through
                   use_parallel = use_parallel, use_opencl = use_opencl, verbose = verbose)

  if (pfamily$pfamily == "dIndependent_Normal_Gamma") {
    if (!is.null(outlist$sim_bounds)) {
      pfamily$prior_list$disp_lower=outlist$sim_bounds$low
      pfamily$prior_list$disp_upper=outlist$sim_bounds$upp
      
    } else {
      cat("No simbounds returned in outlist.\n")
    }
  }  
  
  
  
  outlist$simfun_call <- outlist$call 

  outlist$call <- match.call()  # overwrite with the rglmb call
  outlist$pfamily=pfamily
##  outlist$simfun_call <- simfun_call         # simulation call
  outlist$simfun_args <- simfun_args         # simulation arguments
  
  
  if (is.null(colnames(outlist$coefficients))) {
    colnames(outlist$coefficients) <- colnames(outlist$x)
  }
  
  if (!is.null(outlist$coef.mode) &&
      !is.null(outlist$x) &&
      !is.null(colnames(outlist$x))) {
    
    names(outlist$coef.mode) <- colnames(outlist$x)
  }
  
  return(outlist)
  
}


## print.* S3 method provided by glmbayesCore (Imports); local registration commented out.
# #' @rdname rglmb
# #' @order 2
# #' @method print rglmb
# #' @export
# 
# print.rglmb<-function (x, digits = max(3, getOption("digits") - 3), ...) 
# {
#   
#   cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
#       "\n\n", sep = "")
#   if (length(coef(x))) {
#     cat("Simulated Coefficients")
#     cat(":\n")
#     print.default(format(x$coefficients, digits = digits), 
#                   print.gap = 2, quote = FALSE)
#   }
#   else cat("No coefficients\n\n")
# }





