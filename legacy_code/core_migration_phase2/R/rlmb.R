#' The Bayesian Linear Model Distribution
#'
#' \code{rlmb} is used to generate iid samples from Bayesian Linear Models with multivariate normal priors. 
#' The model is specified by providing a data vector, a design matrix, and a pfamily (determining the 
#' prior distribution).
#' @name rlmb
#' @aliases
#' rlmb
#' print.rlmb
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x for \code{rlmb} a design matrix of dimension \code{m * p} and for 
#' \code{print.rlmb} the object to be printed. 
#' @param pfamily a description of the prior distribution and associated constants to be used in the model. This
#' should be a pfamily function (see \code{\link{pfamily}} for details of pfamily functions.)
#' @param digits the number of significant digits to use when printing.
#' @param progbar Logical. Whether to display a progress base during simulation.
#' @inheritParams lmb
#' @inheritParams glmb
#' @return \code{rlmb} returns a object of class \code{"rlmb"}.  The function \code{summary} 
#' (i.e., \code{\link{summary.rglmb}}) can be used to obtain or print a summary of the results.
#' The generic accessor functions \code{\link{coefficients}}, \code{\link{fitted.values}},
#' \code{\link{residuals}}, and \code{\link{extractAIC}} can be used to extract
#' various useful features of the value returned by \code{\link{rlmb}}.
#' An object of class \code{"rlmb"} is a list containing at least the following components:
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
#' Objects of class \code{"rlmb"} are normally of class \code{c("rlmb","rglmb","glmb","glm","lm")},
#' meaning they inherit from \code{rglmb}, \code{glmb}, \code{glm}, and \code{lm}. Well-designed
#' methods for these classes will be applied when appropriate, allowing \code{"rlmb"} objects to
#' benefit from existing infrastructure while supporting specialized behavior for restricted linear
#' model priors.
#' 
#' @details
#' The function \code{rlmb} is a minimalistic Bayesian simulation engine for Gaussian linear models. 
#' It bypasses classical model fitting and formula parsing, operating directly on numeric inputs such as 
#' the design matrix, response vector, and prior specification via the \code{\link{pfamily}} argument. 
#' Internally, \code{rlmb} generates independent draws from the posterior distribution using multivariate 
#' normal simulation when conjugate priors are specified.
#'
#' The modeling framework follows \insertCite{WilkinsonRogers1973}{glmbayes}, and the prior structure builds on the S system 
#' \insertCite{Chambers1992}{glmbayes}, Zellner's g-prior \insertCite{zellner1986gprior}{glmbayes}, and 
#' the conjugate prior formulation of Raiffa and Schlaifer \insertCite{Raiffa1961}{glmbayes}.
#'
#' Prior specification is handled via the \code{\link{pfamily}} argument, which defines the prior mean, 
#' covariance, and dispersion. The design of the \code{pfamily} family of functions was created by Kjell Nygren 
#' and is modeled on how \code{\link{glm}} uses \code{family} to specify the likelihood. A helper function, 
#' \code{\link{Prior_Setup}}, assists users in choosing prior parameters. It ships with sensible defaults but 
#' also allows full customization. Available priors include the  \code{dNormal}, \code{dNormalGamma} and 
#' \code{dIndependent_Normal_Gamma} priors. The last of these allows for more flexible prior structures 
#' including independent priors on variance components.
#'
#' Posterior draws are generated using standard simulation procedures for conjugate priors \insertCite{Raiffa1961}{glmbayes}. 
#' For non-conjugate setups, the function uses envelope-based accept-reject sampling via the 
#' likelihood-subgradient method \insertCite{Nygren2006}{glmbayes}. The \code{Gridtype} parameter controls 
#' how many tangent points are used to construct the envelope-trading off tightness against computational cost-
#' and the \code{iters} component reports the number of candidate samples generated before acceptance.
#'
#' The output includes posterior samples, prior specifications, dispersion estimates, and envelope diagnostics. 
#' While \code{rlmb} does not return a full model object or support generic methods like \code{predict} or 
#' \code{summary}, it is designed for efficient posterior simulation in Gaussian models where full model 
#' reconstruction is unnecessary.
#'
#' The \code{\link{rlmb}} function called from within \code{\link{lmb}}. 
#' It is intended for simulation-heavy workflows such as Gibbs sampling or posterior 
#' predictive checks where minimal overhead is preferred.
#'  
#' 
#' @family modelfuns
#' @seealso The classical modeling functions \code{\link[stats]{lm}} and \code{\link[stats]{glm}}.
#'
#' \code{\link{lmb}}, \code{\link{glmb}}, \code{\link{rglmb}}
#'   for related interfaces;
#' \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeOrchestrator}} for envelope stages
#' used in non-conjugate Gaussian sampling.
#' 
#' \code{\link{pfamily}} for documentation of pfamily functions used to specify priors.
#' 
#' \code{\link{Prior_Setup}}, \code{\link{Prior_Check}} for functions used to initialize and to check priors,  
#'
#' Further reading: \insertCite{Nygren2006}{glmbayes};
#' \insertCite{glmbayesSimmethods,glmbayesChapterA08,glmbayesIndNormGammaVignette}{glmbayes}.
#'
#' \code{\link{summary.glmb}}, \code{\link{predict.glmb}}, \code{\link{simulate.glmb}}, 
#' \code{\link{extractAIC.glmb}}, \code{\link{dummy.coef.glmb}} and methods(class="glmb") for methods 
#' inherited from class \code{glmb} and the methods and generic functions for classes \code{glm} and 
#' \code{lm} from which class \code{lmb} also inherits.
#'
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_rlmb.R
#' @export
#' @rdname rlmb
#' @order 1

rlmb <- function(
    n = 1,
    y,
    x,
    pfamily,
    offset = rep(0, nobs),
    weights = NULL,
    Gridtype = 2,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE,
    progbar = FALSE
){
  ## Pull in information from the pfamily  
  pf         <- pfamily$pfamily
  okfamilies <- c("gaussian")    # Only gaussian is okfamily for rlmb (different from rglmb)
  plinks     <- pfamily$plinks
  prior_list <- pfamily$prior_list 
  simfun     <- pfamily$simfun
  
  family <- gaussian()
  
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  
  if (length(n) > 1) n <- length(n)
  nobs  <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  
  if (is.null(offset)) 
    offset <- rep(0, nobs)
  
  if (is.null(weights)) weights <- rep(1, nobs)
  if (length(weights) == 1) weights <- rep(weights, nobs)
  nobs2 <- length(weights)
  nobs3 <- NROW(x)
  nobs4 <- NROW(offset)
  
  if (nobs2 != nobs) stop("weighting vector must have same number of elements as y")
  if (nobs3 != nobs) stop("matrix X must have same number of rows as y")
  if (nobs4 != nobs) stop("offset vector must have same number of rows as y")
  
  if (is.null(offset)) 
    offset <- rep.int(0, nobs)
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  if (family$family %in% okfamilies) {
    oklinks <- c("identity")
    if (!family$link %in% oklinks) {      
      stop(gettextf("link \"%s\" not available for selected family; available links are %s", 
                    family$link , paste(sQuote(oklinks), collapse = ", ")), 
           domain = NA)
    }
  } else {
    stop(gettextf("family \"%s\" not available for current pfamily; available families are %s", 
                  family$family , paste(sQuote(okfamilies), collapse = ", ")), 
         domain = NA)
  }
  
  ## normalize n_envopt (like rglmb)
  if (is.null(n_envopt)) n_envopt <- n
  n_envopt <- as.integer(n_envopt)
  
  ## Preserve simfun_args exactly as before, just extended
  simfun_args <- list(
    n          = n,
    y          = y,
    x          = x,
    prior_list = prior_list,
    offset     = offset,
    weights    = weights,
    family     = family,
    Gridtype   = Gridtype,
    n_envopt   = n_envopt,
    use_parallel = use_parallel,
    use_opencl  = use_opencl,
    verbose     = verbose,
    progbar     = progbar
  )
  
  ## Direct call, same style as original
  outlist <- simfun(
    n          = n,
    y          = y,
    x          = x,
    prior_list = prior_list,
    offset     = offset,
    weights    = weights,
    family     = family,
    Gridtype   = Gridtype,
    n_envopt   = n_envopt,
    use_parallel = use_parallel,
    use_opencl  = use_opencl,
    verbose     = verbose,
    progbar     = progbar
  )
  
  if (pfamily$pfamily == "dIndependent_Normal_Gamma") {
    if (!is.null(outlist$sim_bounds)) {
      pfamily$prior_list$disp_lower <- outlist$sim_bounds$low
      pfamily$prior_list$disp_upper <- outlist$sim_bounds$upp
    } else {
      cat("No simbounds returned in outlist.\n")
    }
  }

  if (pfamily$pfamily == "dGamma") {
    class(outlist) <- c("rGamma_reg", "rlmb", "rglmb", "glmb", "glm", "lm")
  } else {
    class(outlist) <- c("rlmb", "rglmb", "glmb", "glm", "lm")
  }
  
  simfun_call <- outlist$call 
  
  outlist$call        <- match.call()  # overwrite with the rlmb call
  outlist$pfamily     <- pfamily
  outlist$simfun_call <- simfun_call   # simulation call
  outlist$simfun_args <- simfun_args   # simulation arguments
  
  return(outlist)
}

## print.* S3 method provided by glmbayesCore (Imports); local registration commented out.
# #' @rdname rlmb
# #' @order 2
# #' @export
# #' @method print rlmb
# 
# 
# print.rlmb<-function (x, digits = max(3, getOption("digits") - 3), ...) 
# {
#  		
#     cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
#         "\n\n", sep = "")
#  if (length(coef(x))) {
#         cat("Simulated Coefficients")
#         cat(":\n")
#         print.default(format(x$coefficients, digits = digits), 
#             print.gap = 2, quote = FALSE)
#     }
#     else cat("No coefficients\n\n")
#  }


