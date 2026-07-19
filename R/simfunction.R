#' Simulation Functions for Bayesian Generalized Linear Models
#'
#' @description
#' Simulation functions provide a unified interface for generating posterior samples from Bayesian GLMs.
#' These functions are typically used within model fitting routines such as \code{\link{rglmb}} and \code{\link{rlmb}}, and
#' are also suitable for use in Block Gibbs sampling and other simulation-based inference techniques.
#'
#' @name simfuncs
#' @param object A fitted model object containing a \code{pfamily} component. The generic function \code{simfunction()} accesses the simulation metadata stored within such objects.
#' @param x An object of class \code{"simfunction"} or \code{"rGamma_reg"} to be printed.
#' @param n Number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y A vector of observations of length \code{m}.
#' @param x for the simulation functions a design matrix of dimension \code{m * p} and for 
#' the print functions the object to be printed. 
#' @param prior_list A list with prior parameters (e.g., shape, rate, beta) used in the simulation.
#' @param offset Optional numeric vector of length \code{m} specifying known components of the linear predictor.
#' @param weights Optional numeric vector of prior weights.
#' @param family A description of the error distribution and link function (see \code{\link{family}}).
#' @param Gridtype Optional integer specifying the method used to construct the envelope function.
#' @param n_envopt Effective sample size passed to EnvelopeOpt for grid
#'   construction. Defaults to match `n`. Larger values encourage tighter
#'   envelopes.
#' @param use_parallel Logical. Whether to use parallel processing.
#' @param use_opencl Logical. Whether to use OpenCL acceleration.
#' @param verbose Logical. Whether to print progress messages.
#' @param digits Number of significant digits to use for printed output.
#' @param progbar Logical. Whether to display a progress base during simulation.
#' @param \ldots Additional arguments passed to or from other methods.
#'
#' @return
#' \describe{
#'
#'   \item{\code{simfunction()}}{An object of class \code{"simfunction"} containing:
#'     \describe{
#'       \item{\code{name}}{Character string with the name of the simulation function.}
#'       \item{\code{call}}{The matched call used to generate the simulation.}
#'       \item{\code{args}}{A named list of arguments passed to the simulation function.}
#'     }
#'   }
#'
#'   \item{\code{rNormal_reg()}}{A list object with classes \code{"rglmb"}, \code{"glmb"}, \code{"glm"}, and \code{"lm"}.
#'     Elements include:
#'     \describe{
#'       \item{\code{coefficients}}{Matrix (\code{n * p}) of simulated regression coefficients, with column names from \code{x}.}
#'       \item{\code{coef.mode}}{Posterior mode of the coefficients. Gaussian: from \code{lm.fit}; non-Gaussian: BFGS mode shifted by prior mean.}
#'       \item{\code{dispersion}}{Scalar dispersion used. Poisson/Binomial: \code{1}; otherwise the supplied value. Quasi families: mean residual-based dispersion computed in the wrapper.}
#'       \item{\code{Prior}}{List with \code{mean} (prior mean vector) and \code{Precision} (prior precision matrix \code{P}).}
#'       \item{\code{prior.weights}}{Vector of prior weights used in the simulation (unscaled).}
#'       \item{\code{offset}}{Offset vector passed to the C++ sampler.}
#'       \item{\code{offset2}}{Offset used internally by the wrapper (copy of input or a zero vector).}
#'       \item{\code{y}}{Response vector.}
#'       \item{\code{x}}{Design matrix.}
#'       \item{\code{fit}}{Fitted/diagnostic object. Gaussian: result of \code{lm.fit} (class \code{"lm"}). Non-Gaussian: result of \code{glmb.wfit(...)}.}
#'       \item{\code{iters}}{Vector of iteration counts per sample. Gaussian: vector of ones; non-Gaussian: counts from the sampler.}
#'       \item{\code{Envelope}}{Envelope list used for accept-reject sampling (non-Gaussian); \code{NULL} for Gaussian.}
#'       \item{\code{family}}{Family object describing distribution and link.}
#'       \item{\code{famfunc}}{Processed family functions used internally (e.g., \code{f2}, \code{f3}).}
#'       \item{\code{call}}{Matched call to \code{rNormal_reg()}.}
#'       \item{\code{formula}}{Formula reconstructed from \code{y} and \code{x}.}
#'       \item{\code{model}}{Model frame corresponding to \code{formula}.}
#'       \item{\code{data}}{Data frame combining \code{y} and \code{x}.}
#'     }
#'   }
#'
#'   \item{\code{rNormalGamma_reg()}}{A list with class \code{"rglmb"} containing:
#'     \describe{
#'       \item{\code{coefficients}}{Matrix (\code{n * p}) of simulated regression coefficients; row \code{i} equals \code{Btilde + IR \%*\% rnorm(p) * sqrt(dispersion[i])}. Column names are set to \code{colnames(x)}.}
#'       \item{\code{coef.mode}}{Posterior mean/mode vector \code{Btilde} from \code{rNormal_reg.wfit()}.}
#'       \item{\code{dispersion}}{Numeric vector of length \code{n} with draws from the inverse-gamma posterior \code{1/rgamma(shape = shape + nobs/2, rate = rate + 0.5*S)}.}
#'       \item{\code{Prior}}{List with \code{mean} (as numeric vector \code{mu}) and \code{Precision} (matrix \code{P}).}
#'       \item{\code{offset}}{Offset vector as supplied.}
#'       \item{\code{prior.weights}}{Vector of prior weights \code{wt}.}
#'       \item{\code{y}}{Response vector.}
#'       \item{\code{x}}{Design matrix.}
#'       \item{\code{fit}}{Result from \code{rNormal_reg.wfit()}, including fields such as \code{Btilde}, \code{IR}, \code{S}, and \code{k}.}
#'       \item{\code{famfunc}}{Processed family functions for Gaussian models (from \code{glmbfamfunc(gaussian())}).}
#'       \item{\code{iters}}{Numeric vector (length \code{n}) of ones indicating per-draw iteration counts.}
#'       \item{\code{Envelope}}{\code{NULL}; no envelope is constructed in this conjugate setup.}
#'       \item{\code{call}}{Matched call to \code{rNormalGamma_reg()}.}
#'     }
#'   }
#'
#'   \item{\code{rindepNormalGamma_reg()}}{A list with class \code{"rglmb"} containing:
#'     \describe{
#'       \item{\code{coefficients}}{Matrix (\code{n * p}) of simulated regression coefficients, back-transformed to the original scale; column names set to \code{colnames(x)}.}
#'       \item{\code{coef.mode}}{Vector with the conditional posterior mode used for envelope anchoring (from the Gaussian fit).}
#'       \item{\code{dispersion}}{Numeric vector of length \code{n} with simulated dispersion draws.}
#'       \item{\code{Prior}}{List with prior components: \code{mean} (prior mean \code{mu}), \code{Sigma} (prior covariance), \code{shape} and \code{rate} (Gamma prior for dispersion), \code{Precision} (\code{solve(Sigma)}).}
#'       \item{\code{family}}{The \code{gaussian()} family object.}
#'       \item{\code{prior.weights}}{Vector of prior weights used in the simulation.}
#'       \item{\code{y}}{Response vector.}
#'       \item{\code{x}}{Design matrix.}
#'       \item{\code{call}}{Matched call to \code{rindepNormalGamma_reg()}.}
#'       \item{\code{famfunc}}{Processed family functions for Gaussian models (from \code{glmbfamfunc}).}
#'       \item{\code{iters}}{Vector with per-draw iteration counts returned by the joint sampler.}
#'       \item{\code{Envelope}}{\code{NULL}; envelope diagnostics are not returned by this function.}
#'       \item{\code{loglike}}{\code{NULL}; placeholder for log-likelihood values.}
#'       \item{\code{weight_out}}{Numeric vector of per-draw weights returned by the C++ routine.}
#'       \item{\code{sim_bounds}}{List with \code{low} and \code{upp}, the dispersion bounds used by the shared envelope.}
#'       \item{\code{offset2}}{Offset vector used internally (copy of input or a zero vector).}
#'     }
#'   }
#'
#'   \item{\code{rGamma_reg()}}{An object of class \code{"rGamma_reg"} containing:
#'     \describe{
#'       \item{\code{coefficients}}{A 1 * p matrix of assumed regression coefficients.}
#'       \item{\code{coef.mode}}{Currently \code{NULL}; reserved for future use.}
#'       \item{\code{dispersion}}{A vector of simulated dispersion values.}
#'       \item{\code{Prior}}{A list with prior parameters: \code{shape} and \code{rate}.}
#'       \item{\code{prior.weights}}{Vector of prior weights used in the simulation.}
#'       \item{\code{y}}{The response vector.}
#'            }
#'            }
#'            }
#'      
#'         
#' @details The low-level simulation functions **\code{rNormal_reg()}**, **\code{rNormalGamma_reg()}**, 
#' **\code{rindepNormalGamma_reg()}**, and **\code{rGamma_reg()}** generate iid samples from posterior 
#' distributions for specific model components. These model functions are used internally by the functions
#' **\code{rglmb()}** and **\code{rlmb()}** to generate samples.  
#'  
#' The \code{simfunction()} generic extracts metadata from simulation objects, including the function name, call, and arguments used. This is useful for introspection, reproducibility, and diagnostics.
#'
#' The lower-level simulation functions generate iid samples from posterior distributions for specific model components. 
#' These functions are used internally by \code{pfamily} constructors and model fitting routines.
#'
#' ## Simulation Functions
#'
#' - **\code{rNormal_reg()}**: Produces iid draws for regression coefficients in models with
#'   multivariate normal priors and log-concave likelihood functions. For Gaussian likelihoods,
#'   these are conjugate priors and standard simulation procedures for multivariate normal
#'   distributions are utilized \insertCite{LindleySmith1972,DiaconisYlvisaker1979}{glmbayes}.
#'   For all other families/link functions, the likelihood subgradient approach of
#'   \insertCite{Nygren2006}{glmbayes} is used to generate iid samples.
#'
#' - **\code{rNormalGamma_reg()}**: Produces iid draws for regression coefficients and the
#'   dispersion parameter in models with Normal-Gamma priors and Gaussian likelihoods, where
#'   this is a conjugate prior distribution. Standard simulation procedures for gamma and
#'   multivariate normal distributions are utilized
#'   \insertCite{Raiffa1961,LindleySmith1972}{glmbayes}.
#'
#' - **\code{rindepNormalGamma_reg()}**: Produces iid draws for regression coefficients
#'   and the dispersion parameter in models with independent Normal and truncated Gamma priors.
#'   This is a non-conjugate specification but can still be sampled using accept-reject procedures
#'   based on an enveloping approach (see vignette
#'   \insertCite{glmbayesIndNormGammaVignette}{glmbayes}).
#'
#' - **\code{rGamma_reg()}**: Simulates dispersion parameters for Gaussian and Gamma families
#'   using either standard gamma sampling or accept-reject methods based on likelihood
#'   subgradients \insertCite{Chen1979,glmbayesGammaVignette}{glmbayes}.
#'
#' @references
#' \insertAllCited{}
#' 
#' @author
#' The simulation framework was developed by Kjell Nygren as part of the \pkg{glmbayes} package. It builds on the likelihood subgradient approach described in \insertCite{Nygren2006}{glmbayes}, and extends classical Bayesian GLM sampling techniques.
#'
#' @seealso
#' \code{\link{pfamily}}, \code{\link{glmb}}, \code{\link{lmb}}, \code{\link{rglmb}}, \code{\link{rlmb}}
#' for modeling functions that consume simulation functions.
#'
#' \code{\link{rNormal_reg}}, \code{\link{rNormalGamma_reg}}, \code{\link{rGamma_reg}} for individual simulation functions.
#'
#' \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeEval}}, \code{\link{EnvelopeSize}} for envelope construction
#' and grid evaluation used in likelihood-subgradient sampling.
#'
#' Theory and implementation narrative: \insertCite{Nygren2006}{glmbayes};
#' \insertCite{glmbayesSimmethods,glmbayesChapterA08}{glmbayes}.
#'
#' @usage simfunction(object, ...)
#' @export
#' @rdname simfuncs
#' @order 1

simfunction <- function(object, ...) {
  UseMethod("simfunction")
}



#' @method simfunction default
#' @noRd
#' @export

simfunction.default <- function(object, ...) {
  if (is.null(object$pfamily)) stop("no pfamily object found")
  if (!inherits(object$pfamily, "pfamily")) stop("Object named pfamily is not of class pfamily")
  
  pf <- object$pfamily
  simfun <- pf$simfun
  
  simfun_name <- "anonymous or not found"
  fun_env <- environment(simfun)
  fun_names <- ls(fun_env)
  for (name in fun_names) {
    if (identical(simfun, get(name, envir = fun_env))) {
      simfun_name <- name
      break
    }
  }
  
  simfun_call <- if (!is.null(object$simfun_call)) object$simfun_call else NULL
  simfun_args <- if (!is.null(object$simfun_args)) object$simfun_args else list()
  
  structure(
    list(
      name = simfun_name,
      call = simfun_call,
      args = simfun_args
    ),
    class = "simfunction"
  )
}

#' @export
#' @method print simfunction
#' @rdname simfuncs
#' @order 9
 
print.simfunction <- function(x, ...) {
  cat("\nCall to Simulation Function:\n")
  if (!is.null(x$call)) {
    print(x$call)
  } else {
    cat("  [call not recorded]\n")
  }
  
  cat("\nSimulation Function Name:", x$name, "\n")
  
  if (!is.null(x$args) && length(x$args) > 0) {
    cat("\nArguments Passed:\n\n")
    for (argname in names(x$args)) {
      val <- x$args[[argname]]
      
      if (is.null(val)) {
        cat("  ", argname, ": [NULL]\n", sep = "")
      } else if (argname == "family") {
        cat("  ", argname, ":\n", sep = "")
        print(val)
      } else if (argname == "prior_list" && is.list(val)) {
        cat("  prior_list:\n")
        for (pname in names(val)) {
          pval <- val[[pname]]
          cat("    ", pname, ":\n", sep = "")
          if (is.null(pval)) {
            cat("      [NULL]\n")
          } else if (is.atomic(pval) || is.matrix(pval)) {
            print(pval)
          } else {
            cat("      [", class(pval), " with length ", length(pval), "]\n", sep = "")
          }
        }
      } else {
        cat("  ", argname, ":\n", sep = "")
        if (is.atomic(val) || is.matrix(val)) {
          print(val)
        } else {
          cat("    [", class(val), " with length ", length(val), "]\n", sep = "")
        }
      }
    }
  } else {
    cat("\nArguments Passed: [none recorded]\n")
  }
  
  invisible(x)
}







#' @family simfuncs
#' @example inst/examples/Ex_rGamma_reg.R
#' @usage rGamma_reg(n, y, x, prior_list, offset = NULL, weights = 1, family = gaussian(),
#'            Gridtype = 2,n_envopt = NULL,
#'             use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE,progbar=FALSE)
#' @export 
#' @rdname simfuncs
#' @order 5
#' @export



rGamma_reg <- function(
    n,
    y,
    x,
    prior_list,
    offset = NULL,
    weights = 1,
    family = gaussian(),
    Gridtype = 2,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE, progbar = FALSE
) {
  call <- match.call()
  
  ## Argument renaming and prior
  wt    <- weights
  alpha <- offset
  
  ## Basic checks on y and x
  n1 <- length(y)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (nrow(x) != n1)
    stop("Number of rows in x must match length of y.")
  
  ## 1) Validate and normalize offset (alpha)
  if (is.null(alpha)) {
    alpha <- rep(0, n1)
  } else {
    if (!is.numeric(alpha))
      stop("offset (alpha) must be numeric if supplied.")
    if (length(alpha) == 1L) {
      alpha <- rep(alpha, n1)
    } else if (length(alpha) != n1) {
      stop("offset (alpha) must be scalar or have length equal to length(y).")
    }
  }
  
  ## 2) Validate and normalize weights (wt)
  if (!is.numeric(wt))
    stop("weights must be numeric.")
  
  if (length(wt) == 1L) {
    wt <- rep(wt, n1)
  } else if (length(wt) != n1) {
    stop("weights must be either a scalar or have length equal to length(y).")
  }
  
  if (any(wt < 0))
    stop("weights must be nonnegative.")
  
  b     <- prior_list$beta
  shape <- prior_list$shape
  rate  <- prior_list$rate

  if (!is.null(prior_list$max_disp_perc)) {
    max_disp_perc <- prior_list$max_disp_perc
  } else {
    max_disp_perc <- 0.99
  }
  
  
  
  ## New: extract optional low/upp from prior_list
  if (!is.null(prior_list$disp_lower))  disp_lower <- prior_list$disp_lower  else disp_lower <- NULL
  if (!is.null(prior_list$disp_upper))  disp_upper <- prior_list$disp_upper  else disp_upper <- NULL
  
  ## Validation if both are provided
  if (!is.null(disp_lower) && !is.null(disp_upper)) {
    if (!is.numeric(disp_lower) || !is.numeric(disp_upper)) {
      stop("prior_list$disp_lower and prior_list$disp_upper must be numeric.")
    }
    if (disp_lower <= 0 || disp_upper <= 0) {
      stop("prior_list$disp_lower and prior_list$disp_upper must be positive.")
    }
    if (disp_upper <= disp_lower) {
      stop("prior_list$disp_upper must be strictly greater than prior_list$disp_lower.")
    }
  }
  
  ## Family handling
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  okfamilies <- c("gaussian", "Gamma")
  if (family$family %in% okfamilies) {
    if (family$family == "gaussian") oklinks <- c("identity")
    if (family$family == "Gamma")    oklinks <- c("log")
    if (!(family$link %in% oklinks)) {
      stop(gettextf(
        "link \"%s\" not available for selected family; available links are %s",
        family$link, paste(sQuote(oklinks), collapse = ", ")
      ), domain = NA)
    }
  } else {
    stop(gettextf(
      "family \"%s\" not available in glmbdisp; available families are %s",
      family$family, paste(sQuote(okfamilies), collapse = ", ")
    ), domain = NA)
  }
  
  ## ----------------------
  ## Gaussian case (updated)
  ## ----------------------
  if (family$family == "gaussian") {
    
    
    sim<-  .rGammaGaussian_cpp(
      n, y, x, b, wt, alpha, shape, rate,
      disp_lower, disp_upper, verbose = verbose
    )
    
    # Validate output
    if (!is.list(sim) || is.null(sim$dispersion) || is.null(sim$draws)) {
      stop("C++ .rGammaGaussian_cpp returned an invalid structure.")
    }
    
    out  <- sim$dispersion
    draws <- sim$draws
    
  }
  
  ## ----------------------
  ## Gamma case (MODIFIED)
  ## ----------------------
  if (family$family == "Gamma") {
    
    # Call C++ sampler — do NOT hide errors
    sim <- .rGammaGamma_cpp(
      n, y, x, b, wt, alpha, shape, rate,
      max_disp_perc, disp_lower, disp_upper,
      verbose = verbose
    )
    
    # Validate output
    if (!is.list(sim) || is.null(sim$dispersion) || is.null(sim$draws)) {
      stop("C++ .rGammaGamma_cpp returned an invalid structure.")
    }
    
    out  <- sim$dispersion
    draws <- sim$draws

  }
  
  outlist <- list(
    coefficients   = matrix(b, nrow = 1, ncol = length(b)),
    coef.mode      = NULL,
    dispersion     = out,
    Prior          = list(shape = shape, rate = rate),
    prior.weights  = wt,
    y              = y,
    x              = x,
    famfunc        = glmbfamfunc(family),
    iters          = draws,
    Envelope       = NULL
  )
  
  outlist$call <- call
  class(outlist) <- c(outlist$class, "rGamma_reg")
  
  return(outlist)
}


#' @export
#' @rdname simfuncs
#' @order 6
#' @method print rGamma_reg

print.rGamma_reg<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Simulated Dispersion")
    cat(":\n")
    print.default(format(x$dispersion, digits = digits), 
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
}



## Internal: intercept-only likelihood + ancillary restrictions for scalar closed-form conjugate updates
## (dGamma_Conjugate with Gamma / Poisson families). Not used for dNormal or other pfamilies.
.check_gamma_conjugate_scalar_design <- function(x, beta, alpha, wt, family_label = "", ctx = "") {
  xf <- if (is.matrix(x)) x else as.matrix(x)
  if (ncol(xf) != 1L) {
    stop(
      gettextf(
        "%sdGamma(Inv_Dispersion=FALSE) supports only a single linear predictor (intercept-only model; ncol(x) == 1). Got ncol(x) = %d for family %s. Use dNormal() for regression with multiple predictors.",
        if (nzchar(ctx)) paste0(ctx, ": ") else "",
        ncol(xf),
        family_label
      ),
      domain = NA
    )
  }
  bb <- tryCatch(as.matrix(beta, ncol = 1L), error = function(e) NULL)
  if (is.null(bb) || nrow(bb) != 1L || ncol(bb) != 1L) {
    stop(
      gettextf(
        "%sprior_list$beta must be a single fixed coefficient (1x1) matching the intercept-only design for dGamma(Inv_Dispersion=FALSE); use dNormal() for vector beta.",
        if (nzchar(ctx)) paste0(ctx, ": ") else ""
      ),
      domain = NA
    )
  }
  ## Offsets shift the linear predictor; scalar conjugate updates are not implemented with non-zero offset.
  if (length(alpha) != nrow(xf)) {
    stop(
      gettextf(
        "%sinternal length mismatch: offset vector and design matrix rows",
        if (nzchar(ctx)) paste0(ctx, ": ") else ""
      ),
      domain = NA
    )
  }
  if (any(!is.finite(alpha)) || max(abs(as.numeric(alpha))) > 1e-12) {
    stop(
      gettextf(
        "%snon-zero `offset` is not implemented for scalar dGamma_Conjugate models (family %s).",
        if (nzchar(ctx)) paste0(ctx, ": ") else "",
        family_label
      ),
      domain = NA
    )
  }
  ## Heterogeneous prior weights/exposures require a conjugate derivation not wired here yet.
  if (length(wt) == nrow(xf) && nrow(xf) > 1L && (max(as.numeric(wt)) - min(as.numeric(wt)) > 1e-12 * max(1, mean(as.numeric(wt))))) {
    stop(
      gettextf(
        "%snon-constant observation `weights` are not implemented for scalar dGamma conjugate (family %s).",
        if (nzchar(ctx)) paste0(ctx, ": ") else "",
        family_label
      ),
      domain = NA
    )
  }
  invisible(TRUE)
}




#' @family simfuncs
#' @usage rGamma_Conjugate_reg(n, y, x, prior_list, offset = NULL, weights = 1, family = gaussian(),
#'            Gridtype = 2,n_envopt = NULL,
#'             use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE,progbar=FALSE)
#' @export
#' @rdname simfuncs
#' @order 10
rGamma_Conjugate_reg <- function(
    n,
    y,
    x,
    prior_list,
    offset = NULL,
    weights = 1,
    family = gaussian(),
    Gridtype = 2,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE, progbar = FALSE
) {
  call <- match.call()

  ## rGamma_Conjugate_reg: scaffold for IID conjugate sampling (intercept-only safeguards below).

  ## Argument renaming and prior
  wt    <- weights
  alpha <- offset
  
  ## Basic checks on y and x
  n1 <- length(y)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (nrow(x) != n1)
    stop("Number of rows in x must match length of y.")
  
  ## 1) Validate and normalize offset (alpha)
  if (is.null(alpha)) {
    alpha <- rep(0, n1)
  } else {
    if (!is.numeric(alpha))
      stop("offset (alpha) must be numeric if supplied.")
    if (length(alpha) == 1L) {
      alpha <- rep(alpha, n1)
    } else if (length(alpha) != n1) {
      stop("offset (alpha) must be scalar or have length equal to length(y).")
    }
  }
  
  ## 2) Validate and normalize weights (wt)
  if (!is.numeric(wt))
    stop("weights must be numeric.")
  
  if (length(wt) == 1L) {
    wt <- rep(wt, n1)
  } else if (length(wt) != n1) {
    stop("weights must be either a scalar or have length equal to length(y).")
  }
  
  if (any(wt < 0))
    stop("weights must be nonnegative.")
  
  b         <- prior_list$beta
  shape     <- prior_list$shape
  rate      <- prior_list$rate
  lik_shape <- as.numeric(if (!is.null(prior_list$lik_shape)) prior_list$lik_shape else 1)[[1L]]

  if (!is.null(prior_list$max_disp_perc)) {
    max_disp_perc <- prior_list$max_disp_perc
  } else {
    max_disp_perc <- 0.99
  }
  
  ## Extract optional dispersion bounds from prior_list
  if (!is.null(prior_list$disp_lower))  disp_lower <- prior_list$disp_lower  else disp_lower <- NULL
  if (!is.null(prior_list$disp_upper))  disp_upper <- prior_list$disp_upper  else disp_upper <- NULL
  
  if (!is.null(disp_lower) && !is.null(disp_upper)) {
    if (!is.numeric(disp_lower) || !is.numeric(disp_upper)) {
      stop("prior_list$disp_lower and prior_list$disp_upper must be numeric.")
    }
    if (disp_lower <= 0 || disp_upper <= 0) {
      stop("prior_list$disp_lower and prior_list$disp_upper must be positive.")
    }
    if (disp_upper <= disp_lower) {
      stop("prior_list$disp_upper must be strictly greater than prior_list$disp_lower.")
    }
  }
  
  ## Family handling
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  ## Scalar conjugate guard: intercept-only, zero offset, constant weights required for
  ## closed-form Gamma–Poisson and Gamma–Gamma(identity) conjugate updates.
  conjugate_scalar_families <- c("poisson", "quasipoisson")
  if (family$family == "Gamma" && family$link == "identity")
    conjugate_scalar_families <- c(conjugate_scalar_families, "Gamma")
  if (family$family %in% conjugate_scalar_families) {
    .check_gamma_conjugate_scalar_design(
      x,
      prior_list$beta,
      alpha,
      wt,
      family_label = family$family,
      ctx = "`rGamma_Conjugate_reg()`"
    )
  }
  
  okfamilies <- c("gaussian", "Gamma", "poisson", "quasipoisson")
  if (family$family %in% okfamilies) {
    oklinks <- if (family$family == "gaussian") {
      c("identity")
    } else if (family$family == "Gamma") {
      ## identity: closed-form Gamma-Gamma conjugate (rate as coefficient).
      ## log: C++ dispersion sampler path.
      c("identity", "log")
    } else if (family$family %in% c("poisson", "quasipoisson")) {
      ## Closed-form conjugate: Gamma prior on λ with identity linear predictor η = λ.
      c("identity")
    }
    if (!(family$link %in% oklinks)) {
      stop(gettextf(
        "link \"%s\" not available for selected family; available links are %s",
        family$link, paste(sQuote(oklinks), collapse = ", ")
      ), domain = NA)
    }
  } else {
    stop(gettextf(
      "family \"%s\" not available in glmbdisp; available families are %s",
      family$family, paste(sQuote(okfamilies), collapse = ", ")
    ), domain = NA)
  }

  ## Number of requested posterior draws (matches `rglmb()` / roxygen convention for `n`).
  n_draw <- if (length(n) > 1L) length(n) else as.integer(n)
  if (!is.finite(n_draw) || n_draw < 1L) {
    stop("`n` must be a positive scalar or a vector whose length defines the number of draws.", call. = FALSE)
  }

  coef_out <- NULL
  disp_out <- NULL
  coef_mode_out <- NULL
  draws_out <- NULL

  ## ----------------------
  ## Gaussian case (updated)
  ## ----------------------
  if (family$family == "gaussian") {


    sim<-  .rGammaGaussian_cpp(
      n, y, x, b, wt, alpha, shape, rate,
      disp_lower, disp_upper, verbose = verbose
    )
    
    # Validate output
    if (!is.list(sim) || is.null(sim$dispersion) || is.null(sim$draws)) {
      stop("C++ .rGammaGaussian_cpp returned an invalid structure.")
    }

    disp_out  <- sim$dispersion
    draws_out <- sim$draws
    coef_out <- matrix(as.numeric(b), nrow = 1L, ncol = length(as.numeric(b)))

  }

  ## ----------------------
  ## Gamma(log): C++ dispersion sampler (existing path; coefficient fixed at beta)
  ## ----------------------
  if (family$family == "Gamma" && family$link == "log") {
    
    sim <- .rGammaGamma_cpp(
      n, y, x, b, wt, alpha, shape, rate,
      max_disp_perc, disp_lower, disp_upper,
      verbose = verbose
    )
    
    if (!is.list(sim) || is.null(sim$dispersion) || is.null(sim$draws)) {
      stop("C++ .rGammaGamma_cpp returned an invalid structure.")
    }
    
    disp_out  <- sim$dispersion
    draws_out <- sim$draws
    coef_out <- matrix(as.numeric(b), nrow = 1L, ncol = length(as.numeric(b)))

  }

  ## ----------------------
  ## Gamma(identity): closed-form Gamma–Gamma conjugate update (rate as coefficient)
  ## ----------------------
  ##
  ## Sampling model (intercept-only, identity link: η = β = rate of the Gamma responses):
  ##   Y_i | β ~ Gamma(shape = lik_shape, rate = β),  i = 1,...,n_obs.
  ## lik_shape k is the *known* Gamma likelihood shape (default 1 = exponential).
  ##
  ## Prior on β (R's shape-rate Gamma):
  ##   β ~ Gamma(shape = alpha0, rate = beta0)
  ##
  ## Likelihood kernel in β:
  ##   prod_i β^k * exp(-β * y_i)  =  β^(n*k) * exp(-β * sum_i y_i)
  ##
  ## Conjugate posterior:
  ##   β | y ~ Gamma(shape = alpha0 + n*k,  rate = beta0 + sum_i y_i).
  ##
  ## Note: the posterior *mean of β* = (alpha0 + n*k) / (beta0 + sum_y);
  ##       the posterior *mean of the Gamma response mean* = k / E[β] (not directly available
  ##       from β alone without further integration, but the mode of 1/β ≈ k / mode(β)).
  if (family$family == "Gamma" && family$link == "identity") {
    if (!isTRUE(all(y > 0))) {
      stop("`rGamma_Conjugate_reg()`: Gamma(identity) responses must be strictly positive.", call. = FALSE)
    }
    shape0     <- as.numeric(shape)[[1L]]
    rate0      <- as.numeric(rate)[[1L]]
    sum_y      <- sum(as.numeric(y))
    sum_w      <- sum(as.numeric(wt))
    shape_post <- shape0 + sum_w * lik_shape
    rate_post  <- rate0  + sum_y

    coef_out <- matrix(
      stats::rgamma(n_draw, shape = shape_post, rate = rate_post),
      nrow = n_draw,
      ncol = 1L
    )

    coef_mode_out <- matrix(
      if (shape_post > 1) ((shape_post - 1) / rate_post) else NA_real_,
      nrow = 1L,
      ncol = 1L
    )

    ## Dispersion for Gamma GLM = 1/lik_shape; stored as a constant across draws.
    disp_out  <- rep(1 / lik_shape, n_draw)
    draws_out <- matrix(1L, nrow = n_draw, ncol = 1L)

    coef_nm <- if (!is.null(colnames(x)) && nzchar(colnames(x)[[1L]])) colnames(x)[[1L]] else "(Intercept)"
    colnames(coef_out)      <- coef_nm
    colnames(coef_mode_out) <- coef_nm
    rownames(coef_mode_out) <- "(Mode)"
  }

  ## ----------------------
  ## Poisson / quasi-Poisson (identity link, Gamma prior on rate lambda): closed-form update
  ## ----------------------
  ##
  ## Scalar intercept-only linear predictor with identity link: eta = lambda (one column in `x`,
  ## offset zero; enforced by `.check_gamma_conjugate_scalar_design()`).
  ##
  ## Sampling model (constant rate across observations; prior weights act as identical exposures):
  ##   Y_i | lambda ~ Poisson(lambda * omega_i),  i = 1,...,n_obs.
  ## Constant omega (for n_obs > 1) gives sum_i omega_i = n_obs * omega.
  ##
  ## Prior on lambda (R's shape-rate Gamma: see `?rgamma`; pdf proportional to
  ## lambda^(alpha0-1) * exp(-beta0 * lambda)):
  ##   lambda ~ Gamma(shape = alpha0, rate = beta0)
  ## using prior_list$shape and prior_list$rate.
  ##
  ## Likelihood kernel in lambda (dropping factors not depending on lambda):
  ##   prod_i lambda^y_i exp(-omega_i lambda)
  ##     = lambda^(sum_i y_i) * exp(-lambda * sum_i omega_i).
  ##
  ## Conjugate posterior:
  ##   lambda | y ~ Gamma(shape = alpha0 + sum_i y_i,  rate = beta0 + sum_i omega_i).
  ##
  ## Each requested draw is IID: stats::rgamma(1, shape = shape_post, rate = rate_post), repeated
  ## n_draw times.  The vector `dispersion` mirrors `family$dispersion` (1 for Poisson; fixed
  ## quasi factor if present) and is not assigned a sampling role here.

  if (family$family %in% c("poisson", "quasipoisson")) {
    if (!isTRUE(all(y >= 0))) {
      stop("`rGamma_Conjugate_reg()`: Poisson / quasi-Poisson responses must be nonnegative.", call. = FALSE)
    }
    shape0 <- as.numeric(shape)[[1L]]
    rate0 <- as.numeric(rate)[[1L]]
    sum_y <- sum(as.numeric(y))
    sum_w <- sum(as.numeric(wt))
    shape_post <- shape0 + sum_y
    rate_post <- rate0 + sum_w

    coef_out <- matrix(
      stats::rgamma(n_draw, shape = shape_post, rate = rate_post),
      nrow = n_draw,
      ncol = 1L
    )

    coef_mode_out <- matrix(
      if (shape_post > 1) ((shape_post - 1) / rate_post) else NA_real_,
      nrow = 1L,
      ncol = 1L
    )

    phi <- suppressWarnings(as.numeric(family$dispersion))
    if (length(phi) < 1L || !is.finite(phi[[1L]]) || phi[[1L]] <= 0) {
      phi <- 1
    } else {
      phi <- phi[[1L]]
    }
    disp_out <- rep(phi, n_draw)

    draws_out <- matrix(1L, nrow = n_draw, ncol = 1L)

    coef_nm <- if (!is.null(colnames(x)) && nzchar(colnames(x)[[1L]])) colnames(x)[[1L]] else "(Intercept)"
    colnames(coef_out) <- coef_nm
    colnames(coef_mode_out) <- coef_nm
    rownames(coef_mode_out) <- "(Mode)"
  }

  if (is.null(coef_out) || is.null(disp_out) || is.null(draws_out)) {
    stop("`rGamma_Conjugate_reg()`: incomplete sampler output.", call. = FALSE)
  }

  closed_form_path <- family$family %in% c("poisson", "quasipoisson") ||
                      (family$family == "Gamma" && family$link == "identity")
  if (closed_form_path) {
    if (nrow(coef_out) != n_draw || length(disp_out) != n_draw || nrow(draws_out) != n_draw) {
      stop(
        "`rGamma_Conjugate_reg()`: draw dimension mismatch in closed-form conjugate path.",
        call. = FALSE
      )
    }
  }

  ## Assemble returned object (class `"rglmb"`):
  ## Gaussian / Gamma conjugate dispersion paths reuse fixed `prior_list$beta` in each row dimension
  ## (historic layout: one coefficient row paired with dispersion vector of length `n_draw`).
  ## Poisson / quasi-Poisson: each row is an IID λ draw (`identity` link ⇒ intercept equals λ);
  ## dispersion holds nominal φ (not a conjugate dispersion draw here).

  sh <- as.numeric(shape)[[1L]]
  rt <- as.numeric(rate)[[1L]]
  if (!is.finite(sh) || !is.finite(rt) || sh <= 0 || rt <= 0) {
    stop("`rGamma_Conjugate_reg()`: prior `shape` and `rate` must be finite and positive.", call. = FALSE)
  }
  ## Moments for Gamma(shape, rate) (stats parameterization): mean = sh/rt, variance = sh/rt^2,
  ## precision = rt^2/sh.  `summary.rglmb()` expects `Prior$mean` and `Prior$Precision` compatible
  ## with `ncol(x)` coefficient columns (intercept-only scalar conjugate repeats the same marginal).
  cn_x <- colnames(x)
  if (is.null(cn_x)) {
    cn_x <- paste0("V", seq_len(ncol(x)))
  }
  prior_mean <- stats::setNames(rep(sh / rt, length.out = ncol(x)), cn_x)
  prior_prec_diag <- rep((rt * rt) / sh, length.out = ncol(x))
  prior_precision <- diag(prior_prec_diag, nrow = ncol(x), ncol = ncol(x))
  dimnames(prior_precision) <- list(cn_x, cn_x)

  outlist <- list(
    coefficients   = coef_out,
    coef.mode      = coef_mode_out,
    dispersion     = disp_out,
    Prior          = list(
      shape      = shape,
      rate       = rate,
      mean       = prior_mean,
      Precision  = prior_precision
    ),
    prior.weights  = wt,
    y              = y,
    x              = x,
    famfunc        = glmbfamfunc(family, lik_shape = lik_shape),
    iters          = draws_out,
    Envelope       = NULL
  )

  colnames(outlist$coefficients) <- colnames(x)

  ## Structural parity with `rNormal_reg()`/`rNormalGLM` `rglmb` returns (`family`, `offset2`,
  ## reconstructed model objects). Stochastic draws above are unchanged.
  offset2 <- alpha
  rglmb_df <- as.data.frame(cbind(y, x))
  rglmb_f <- DF2formula(rglmb_df)
  rglmb_mf <- stats::model.frame(rglmb_f, rglmb_df)
  outlist$family <- family
  outlist$offset2 <- offset2
  outlist$formula <- rglmb_f
  outlist$model <- rglmb_mf
  outlist$data <- rglmb_df

  outlist$call <- call
  class(outlist) <- c(outlist$class, "rglmb", "glmb", "glm", "lm")

  return(outlist)
}


## ---------------------------------------------------------------------------
## rBeta_reg: IID conjugate sampler for Beta–Binomial(identity) models.
##
## Model:  Y_i | θ ~ Binomial(n_i, θ),  θ = Xβ  (identity link, intercept only)
##         θ ~ Beta(shape1, shape2)
## Posterior: θ | y ~ Beta(shape1 + Σ n_i y_i,  shape2 + Σ n_i (1 − y_i))
## ---------------------------------------------------------------------------

#' @family simfuncs
#' @usage rBeta_reg(n, y, x, prior_list, offset = NULL, weights = 1,
#'                 family = gaussian(), Gridtype = 2, n_envopt = NULL,
#'                 use_parallel = TRUE, use_opencl = FALSE,
#'                 verbose = FALSE, progbar = FALSE)
#' @export
#' @rdname simfuncs

rBeta_reg <- function(
    n,
    y,
    x,
    prior_list,
    offset       = NULL,
    weights      = 1,
    family       = gaussian(),
    Gridtype     = 2,
    n_envopt     = NULL,
    use_parallel = TRUE,
    use_opencl   = FALSE,
    verbose      = FALSE,
    progbar      = FALSE
) {
  call <- match.call()

  wt    <- weights
  alpha <- offset

  ## Basic checks on y and x
  n1 <- length(y)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (nrow(x) != n1)
    stop("Number of rows in x must match length of y.")

  ## Normalize offset
  if (is.null(alpha)) {
    alpha <- rep(0, n1)
  } else {
    if (!is.numeric(alpha))
      stop("offset (alpha) must be numeric if supplied.")
    if (length(alpha) == 1L) alpha <- rep(alpha, n1)
    else if (length(alpha) != n1)
      stop("offset (alpha) must be scalar or have length equal to length(y).")
  }

  ## Normalize weights (= number of trials for binomial)
  if (!is.numeric(wt)) stop("weights must be numeric.")
  if (length(wt) == 1L) {
    wt <- rep(wt, n1)
  } else if (length(wt) != n1) {
    stop("weights must be either a scalar or have length equal to length(y).")
  }
  if (any(wt < 0)) stop("weights must be nonneg.")

  ## Extract prior hyperparameters
  shape1 <- prior_list$shape1
  shape2 <- prior_list$shape2
  b      <- prior_list$beta

  if (is.null(shape1) || is.null(shape2))
    stop("`rBeta_reg()`: prior_list must contain `shape1` and `shape2`.")
  s1 <- as.numeric(shape1)[[1L]]
  s2 <- as.numeric(shape2)[[1L]]

  ## Family handling
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (is.null(family$family)) stop("'family' not recognized")

  ## Family / link guard
  if (!family$family %in% c("binomial", "quasibinomial"))
    stop("`rBeta_reg()`: only binomial / quasibinomial families are supported.")
  if (family$link != "identity")
    stop("`rBeta_reg()`: only the identity link is supported (theta = Xbeta).")

  ## Scalar design guard (intercept-only, zero offset, constant weights)
  if (ncol(x) != 1L)
    stop("`rBeta_reg()`: dBeta supports only intercept-only models (ncol(x) == 1). Use dNormal() for multiple predictors.")
  bb <- tryCatch(as.matrix(b, ncol = 1L), error = function(e) NULL)
  if (is.null(bb) || nrow(bb) != 1L || ncol(bb) != 1L)
    stop("`rBeta_reg()`: prior_list$beta must be a 1x1 matrix for the intercept-only model.")
  if (any(!is.finite(alpha)) || max(abs(as.numeric(alpha))) > 1e-12)
    stop("`rBeta_reg()`: non-zero offset is not implemented for dBeta models.")
  if (n1 > 1L &&
      (max(as.numeric(wt)) - min(as.numeric(wt))) > 1e-12 * max(1, mean(as.numeric(wt))))
    stop("`rBeta_reg()`: non-constant weights are not implemented for dBeta. Use dNormal() if trial counts vary.")

  ## Response guard: y must be proportions in [0, 1]
  y_num <- as.numeric(y)
  if (any(!is.finite(y_num)))
    stop("`rBeta_reg()`: response y contains non-finite values.")
  if (any(y_num < 0) || any(y_num > 1))
    stop("`rBeta_reg()`: response y must be in [0, 1] (proportions) for binomial(identity).")

  ## Conjugate posterior update
  ## With constant weights n (number of trials):
  ##   shape1_post = shape1 + sum(n * y) = shape1 + total successes
  ##   shape2_post = shape2 + sum(n * (1-y)) = shape2 + total failures
  sum_successes <- sum(y_num * wt)
  sum_failures  <- sum((1 - y_num) * wt)

  shape1_post <- s1 + sum_successes
  shape2_post <- s2 + sum_failures

  n_draw    <- as.integer(n)
  coef_out  <- matrix(
    stats::rbeta(n_draw, shape1 = shape1_post, shape2 = shape2_post),
    nrow = n_draw, ncol = 1L
  )
  disp_out  <- rep(1, n_draw)          ## binomial dispersion = 1
  draws_out <- matrix(1L, nrow = n_draw, ncol = 1L)

  ## Posterior mode: (shape1_post - 1) / (shape1_post + shape2_post - 2) when both > 1
  s12_post <- shape1_post + shape2_post
  if (shape1_post > 1 && shape2_post > 1) {
    mode_val <- (shape1_post - 1) / (s12_post - 2)
  } else {
    mode_val <- shape1_post / s12_post   ## fall back to posterior mean
  }
  coef_mode_out <- matrix(mode_val, nrow = 1L, ncol = 1L)
  rownames(coef_mode_out) <- "(Mode)"

  ## Prior moments for output Prior slot
  cn_x <- colnames(x)
  if (is.null(cn_x)) cn_x <- paste0("V", seq_len(ncol(x)))
  prior_mean_val <- s1 / (s1 + s2)
  prior_var_val  <- s1 * s2 / ((s1 + s2)^2 * (s1 + s2 + 1))
  prior_prec_val <- 1 / prior_var_val
  prior_mean_vec <- stats::setNames(rep(prior_mean_val, ncol(x)), cn_x)
  prior_precision <- diag(rep(prior_prec_val, ncol(x)), nrow = ncol(x), ncol = ncol(x))
  dimnames(prior_precision) <- list(cn_x, cn_x)

  outlist <- list(
    coefficients  = coef_out,
    coef.mode     = coef_mode_out,
    dispersion    = disp_out,
    Prior         = list(
      shape1     = shape1,
      shape2     = shape2,
      mean       = prior_mean_vec,
      Precision  = prior_precision
    ),
    prior.weights = wt,
    y             = y,
    x             = x,
    famfunc       = glmbfamfunc(family),
    iters         = draws_out,
    Envelope      = NULL
  )

  colnames(outlist$coefficients) <- colnames(x)

  offset2  <- alpha
  rglmb_df <- as.data.frame(cbind(y, x))
  rglmb_f  <- DF2formula(rglmb_df)
  rglmb_mf <- stats::model.frame(rglmb_f, rglmb_df)
  outlist$family  <- family
  outlist$offset2 <- offset2
  outlist$formula <- rglmb_f
  outlist$model   <- rglmb_mf
  outlist$data    <- rglmb_df

  outlist$call <- call
  class(outlist) <- c(outlist$class, "rglmb", "glmb", "glm", "lm")

  return(outlist)
}



#' @family simfuncs
#' @example inst/examples/Ex_rindepNormalGamma_reg.R
#' @usage rindepNormalGamma_reg(n, y, x, prior_list, offset = NULL, weights = 1,
#'                              family = gaussian(), Gridtype = 2,n_envopt = NULL,
#'                               use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, 
#'                              progbar = TRUE)
#' @export 
#' @rdname simfuncs
#' @order 4



rindepNormalGamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                                      Gridtype=2,n_envopt = NULL,
                                      use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE,
                                      progbar=TRUE){



  call<-match.call()
  
  offset2=offset
  wt=weights
  
  if(length(wt)==1) wt=rep(wt,length(y))
  
  ### Initial implementation of Likelihood subgradient Sampling 
  ### Currently uses as single point for conditional tangencis
  ### (at conditional posterior modes)
  ### Verify this yields correct results and then try to implement grid approach
  
  ## Use the prior list to set the prior elements if it is not missing
  ## Error checking to verify that the correct elements are present
  ## Shold be implemented
  
  
  ## Step 1: Validate Prior Specification
  
  if(missing(prior_list)) stop("Prior Specification Missing")
  if(!missing(prior_list)){
    if(!is.null(prior_list$mu)) mu=prior_list$mu
    if(!is.null(prior_list$Sigma)) Sigma=prior_list$Sigma
    if(!is.null(prior_list$dispersion)) dispersion=prior_list$dispersion
    else dispersion=NULL
    if(!is.null(prior_list$shape)) shape=prior_list$shape
    else shape=NULL
    if(!is.null(prior_list$rate)) rate=prior_list$rate
    else rate=NULL
    if (!is.null(prior_list$max_disp_perc)) {
      max_disp_perc <- prior_list$max_disp_perc
    } else {
      max_disp_perc <- 0.99
    }
    
    ## New: extract optional low/upp from prior_list
    if (!is.null(prior_list$disp_lower))  disp_lower <- prior_list$disp_lower  else disp_lower <- NULL
    if (!is.null(prior_list$disp_upper))  disp_upper <- prior_list$disp_upper  else disp_upper <- NULL
    
    ## Validation if both are provided
    if (!is.null(disp_lower) && !is.null(disp_upper)) {
      if (!is.numeric(disp_lower) || !is.numeric(disp_upper)) {
        stop("prior_list$disp_lower and prior_list$disp_upper must be numeric.")
      }
      if (disp_lower <= 0 || disp_upper <= 0) {
        stop("prior_list$disp_lower and prior_list$disp_upper must be positive.")
      }
      if (disp_upper <= disp_lower) {
        stop("prior_list$disp_upper must be strictly greater than prior_list$disp_lower.")
      }
    }
    
  }
  

  
  # Reconstruct P from Sigma
  R <- chol(Sigma)
  P <- chol2inv(R)
  P <- 0.5 * (P + t(P))
  
  
  ##########################  BEGIN *.CPP  MIGRATION   #########################################################

  ## --- NEW: Normalize and validate inputs before calling Rcore ---
  
  # Coerce basic types
  y  <- as.numeric(y)
  x  <- as.matrix(x)
  mu <- as.numeric(mu)
  wt <- as.numeric(wt)
  
  # Normalize offset
  if (is.null(offset2)) offset2 <- rep(0, length(y))
  offset2 <- as.numeric(offset2)
  
  # Normalize weights
  if (length(wt) == 1L) wt <- rep(wt, length(y))
  stopifnot(length(wt) == length(y))
  
  # Dimension checks
  stopifnot(nrow(x) == length(y))
  stopifnot(length(mu) == ncol(x))
  
  ## Prior-vs-data balance guard. The dispersion envelope caps the log-tilt
  ## at n_w/2 (the data contribution to shape2; see Remark 4.1.3 of the ING
  ## vignette), which presumes a likelihood-dominated regime. Under the
  ## Prior_Setup ING calibration shape = (n_prior + 1 + p)/2, so the implied
  ## effective prior sample size must not exceed the weighted observation
  ## count n_w = sum(wt) (equivalently pwt <= 0.5).
  if (!is.null(shape)) {
    n_w <- sum(wt)
    n_prior_implied <- 2 * shape - 1 - ncol(x)
    if (n_prior_implied > n_w) {
      stop(
        "dIndependent_Normal_Gamma prior implies n_prior = ",
        signif(n_prior_implied, 4), " effective prior observations, but the ",
        "data supply only n_w = sum(weights) = ", signif(n_w, 4), ". The ",
        "dispersion envelope requires n_prior <= n_w (prior weight pwt <= 0.5); ",
        "weaken the prior (smaller shape) or supply more data.",
        call. = FALSE
      )
    }
  }
  
  # Reconstruct P from Sigma and enforce SPD
  R    <- chol(Sigma)
  Pinv <- chol2inv(R)
  P    <- 0.5 * (Pinv + t(Pinv))
  
  stopifnot(isSymmetric(P))
  
  tol <- 1e-6
  ev  <- eigen(P, symmetric = TRUE)$values
  stopifnot(all(ev >= -tol * abs(ev[1L])))
  
  # dispersion must be numeric scalar or NULL
  if (!is.null(dispersion)) {
    dispersion <- as.numeric(dispersion)
    stopifnot(length(dispersion) == 1L, is.finite(dispersion))
  }
  
  if (is.null(n_envopt)) n_envopt <- n
  n_envopt <- as.integer(n_envopt)

  
  core_out <- .rIndepNormalGammaReg_cpp(
    n,
    y,
    x,
    mu,
    P,
    offset2,
    wt,
    shape,
    rate,
    max_disp_perc,
    disp_lower,
    disp_upper,
    Gridtype,
    n_envopt,
    use_parallel,
    use_opencl,
    verbose,
    progbar
  )
  


  
  out        <- core_out$out
  betastar   <- core_out$betastar
  disp_out   <- core_out$disp_out
  iters_out  <- core_out$iters_out
  weight_out <- core_out$weight_out
  low        <- core_out$low
  upp        <- core_out$upp
  
  famfunc=glmbfamfunc(gaussian())  
  f1=famfunc$f1
  
  R <- chol(Sigma)
  Prec <- chol2inv(R)
  Prec <- 0.5 * (Prec + t(Prec))   # enforce symmetry
  

  outlist=list(
    coefficients=t(out), 
    coef.mode=betastar,  ## For now, use the conditional mode (not universal)
    dispersion=disp_out,
    ## For now, name items in list like this-eventually make format/names
    ## consistent with true prior (current names needed by summary function)
    Prior=list(mean=mu,Sigma=Sigma,shape=shape,rate=rate,Precision=Prec), 
    family=gaussian(),
    prior.weights=wt,
    y=y,
    x=x,
    call=call,
    famfunc=famfunc,
    iters=iters_out,
    Envelope=NULL,
    loglike=NULL,
    weight_out=weight_out,
    sim_bounds=list(low=low,upp=upp)
    #,test_out=test_out
  )
  
  ## Build a minimal pfamily object so summary.rglmb can detect prior type
  pfamily_obj <- list(
    pfamily = "dIndependent_Normal_Gamma",
    prior_list = list(
      mu = mu,
      Sigma = Sigma,
      dispersion = dispersion,
      shape = shape,
      rate = rate,
      max_disp_perc = max_disp_perc,
      disp_lower = low,
      disp_upper = upp
    ),
    pfun = rIndependent_Normal_Gamma_prior
  )
  attr(pfamily_obj, "Prior Type") <- "dIndependent_Normal_Gamma"
  class(pfamily_obj) <- "pfamily"
  outlist$pfamily <- pfamily_obj
  
  colnames(outlist$coefficients)<-colnames(x)
  outlist$offset2<-offset2
  class(outlist)<-c(outlist$class,"rglmb")
  
  return(outlist)  
  
  
}


#' @family simfuncs 
#' @example inst/examples/Ex_rNormalGamma_reg.R
#' @usage rNormalGamma_reg(n, y, x, prior_list, offset = NULL, weights = 1, family = gaussian(),
#'                   Gridtype = 2,n_envopt = NULL, 
#'                   use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE,progbar=FALSE)
#' @export 
#' @rdname simfuncs
#' @order 3


rNormalGamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                            Gridtype=2,n_envopt = NULL,
                            use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, progbar = FALSE
){


  ## Added for consistency with earlier verion of function
  
  offset2=offset
  wt=weights
  
  ## Below code used precision matrix (not Sigma)
  ## Code checks for the presence of P in the prior
  ## if not present, it imputes Precision by inverting the Sigma matrix
  
  if(missing(prior_list)) stop("Prior Specification Missing")
  if(!missing(prior_list)){
    if(!is.null(prior_list$mu)) mu=prior_list$mu
    if (!is.null(prior_list$Sigma)) {
      Sigma <- prior_list$Sigma
      if (!isSymmetric(Sigma)) 
        stop("matrix Sigma must be symmetric")
    }
    
    if (!is.null(prior_list$P)) {
      P <- prior_list$P
      if (!isSymmetric(P)) 
        stop("matrix P must be symmetric")
    }
    
    if (is.null(prior_list$P)) {
      R <- chol(prior_list$Sigma)
      P <- chol2inv(R)
      P <- 0.5 * (P + t(P))   # enforce symmetry
    }
    if(!is.null(prior_list$dispersion)) dispersion=prior_list$dispersion
    else dispersion=NULL
    if(!is.null(prior_list$shape)) shape=prior_list$shape
    else shape=NULL
    if(!is.null(prior_list$rate)) rate=prior_list$rate
    else rate=NULL
  }
  
  if(is.numeric(n)==FALSE||is.numeric(y)==FALSE||is.numeric(x)==FALSE||
     is.numeric(mu)==FALSE||is.numeric(P)==FALSE) stop("non-numeric argument to numeric function")
  
  x <- as.matrix(x)
  mu<-as.matrix(as.vector(mu))
  P<-as.matrix(P)    
  
  ## Allow function to be called without offset2
  
  if(length(n)>1) n<-length(n)	   
  nobs <- NROW(y)
  nvars <- ncol(x)
  
  if(is.null(offset2)) offset2=rep(0,nobs)
  nvars2<-length(mu)	
  if(!nvars==nvars2) stop("incompatible dimensions")
  if (!all(dim(P) == c(nvars2, nvars2))) 
    stop("incompatible dimensions")
  if(!isSymmetric(P))stop("matrix P must be symmetric")
  if(length(wt)==1) wt=rep(wt,nobs)
  if(any(wt < 0)) stop("weights must be non-negative")
  nobs2=NROW(wt)
  nobs3=NROW(x)
  nobs4=NROW(offset2)
  if(nobs2!=nobs) stop("weighting vector must have same number of elements as y")
  if(nobs3!=nobs) stop("matrix X must have same number of rows as y")
  if(nobs4!=nobs) stop("offset vector must have same number of rows as y")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(P, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'P' is not positive definite")
  

  ## Add Call to new function here or after famfunc / f1 is set (not sure if f1 is actually used)
  
  sim <- .rNormalGammaReg_cpp(
    n          = n,
    y          = y,
    x          = x,
    mu         = mu,
    P          = P,
    offset     = offset2,
    wt         = wt,
    shape      = shape,
    rate       = rate,
    max_disp_perc = NULL,
    disp_lower = NULL,
    disp_upper = NULL,
    verbose    = verbose
  )
  
  # Extract components returned by C++
  out1        <- sim$coefficients
  Btilde      <- sim$coef.mode
  dispersion  <- sim$dispersion
  fit         <- sim$fit
  famfunc     <- sim$famfunc
  draws       <- sim$iters
  
  # Build final outlist (mirrors original structure)
  outlist <- list(
    coefficients   = out1,
    coef.mode      = Btilde,
    dispersion     = dispersion,
    family         = family,
    offset         = offset,
    Prior          = list(mean = as.numeric(mu), Precision = P),
    prior.weights  = wt,
    y              = y,
    x              = x,
    fit            = fit,
    famfunc        = famfunc,
    iters          = draws,
    Envelope       = NULL
  )
  
  ## Build a minimal pfamily object so summary.rglmb can detect prior type
  R_pl <- chol(P)
  Sigma_pl <- chol2inv(R_pl)
  Sigma_pl <- 0.5 * (Sigma_pl + t(Sigma_pl))
  pfamily_obj <- list(
    pfamily = "dNormal_Gamma",
    prior_list = list(
      mu = as.numeric(mu),
      Sigma = Sigma_pl,
      P = P,
      shape = shape,
      rate = rate
    ),
    pfun = rNormal_Gamma_prior
  )
  attr(pfamily_obj, "Prior Type") <- "dNormal_Gamma"
  class(pfamily_obj) <- "pfamily"
  outlist$pfamily <- pfamily_obj
  
  colnames(outlist$coefficients) <- colnames(x)
  
  outlist$call <- match.call()
  class(outlist) <- c(outlist$class, "rglmb")
  
  return(outlist)
  

}



#' @family simfuncs 
#' @example inst/examples/Ex_rNormal_reg.R
#' @usage rNormal_reg(n, y, x, prior_list, offset = NULL, weights = 1,
#'             family = gaussian(), Gridtype = 2, n_envopt = NULL,
#'             use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE,progbar=FALSE)
#' @export 
#' @rdname simfuncs
#' @order 2



rNormal_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                      Gridtype=2,n_envopt = NULL,
                      use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE,progbar=FALSE){
  
  ## Added for consistency with earlier verion of function
  ## Useful to copy offset and weight and then to modify to non-null as needed  
  
  offset2=offset
  wt=weights
  
  ## Missing control variables (add option to pass these)
  ## Setting for default Gridtype might be important
  
  #Gridtype=2
  
  
  ## Below code used precision matrix (not Sigma)
  ## Code checks for the presence of P in the prior
  ## if not present, it imputes Precision by inverting the Sigma matrix
  
  if(missing(prior_list)) stop("Prior Specification Missing")
  if(!missing(prior_list)){
    if(!is.null(prior_list$mu)) mu=prior_list$mu
    if(!is.null(prior_list$Sigma)) Sigma=prior_list$Sigma
    if(!is.null(prior_list$P)) P=prior_list$P

    if (is.null(prior_list$P)) {
      R <- chol(prior_list$Sigma)
      Pinv <- chol2inv(R)
      P <- 0.5 * (Pinv + t(Pinv))   # enforce symmetry
    }

    if(!is.null(prior_list$dispersion)) dispersion=prior_list$dispersion
    else dispersion=NULL
    if(!is.null(prior_list$shape)) shape=prior_list$shape
    else shape=NULL
    if(!is.null(prior_list$rate)) rate=prior_list$rate
    else rate=NULL
  }
  
  
  if(is.numeric(n)==FALSE||is.numeric(y)==FALSE||is.numeric(x)==FALSE||
     is.numeric(mu)==FALSE||is.numeric(P)==FALSE) stop("non-numeric argument to numeric function")
  
  # normalize n_envopt
  if (is.null(n_envopt)) n_envopt <- n
  n_envopt <- as.integer(n_envopt)
  
  
  x <- as.matrix(x)
  mu<-as.matrix(as.vector(mu))
  P<-as.matrix(P)    
  
  ## Start value should be contingent on the family and link
  
  start <- mu
  
  ## Allow function to be called without offset2
  
  if(length(n)>1) n<-length(n)	   
  nobs <- NROW(y)
  nvars <- ncol(x)
  
  if(is.null(offset2)) offset2=rep(0,nobs)
  nvars2<-length(mu)	
  if(!nvars==nvars2) stop("incompatible dimensions")
  if (!all(dim(P) == c(nvars2, nvars2))) 
    stop("incompatible dimensions")
  if(!isSymmetric(P))stop("matrix P must be symmetric")
  if(length(wt)==1) wt=rep(wt,nobs)
  nobs2=NROW(wt)
  nobs3=NROW(x)
  nobs4=NROW(offset2)
  if(nobs2!=nobs) stop("weighting vector must have same number of elements as y")
  if(nobs3!=nobs) stop("matrix X must have same number of rows as y")
  if(nobs4!=nobs) stop("offset vector must have same number of rows as y")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(P, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'P' is not positive definite")
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  okfamilies <- c("gaussian","poisson","binomial","quasipoisson","quasibinomial","Gamma")
  if(family$family %in% okfamilies){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-c("log")		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-c("logit","probit","cloglog")		
    if(family$family=="Gamma") oklinks<-c("log")		
    if(family$link %in% oklinks){
      
      famfunc<-glmbfamfunc(family)
      f1<-famfunc$f1
      f2<-famfunc$f2
      f3<-famfunc$f3
      #      f5<-famfunc$f5
      #      f6<-famfunc$f6
    }
    else{
      stop(gettextf("link \"%s\" not available for selected family; available links are %s", 
                    family$link , paste(sQuote(oklinks), collapse = ", ")), 
           domain = NA)
      
    }	
    
  }		
  else {
    stop(gettextf("family \"%s\" not available in glmb; available families are %s", 
                  family$family , paste(sQuote(okfamilies), collapse = ", ")), 
         domain = NA)
  }
  
  ## ddef: from dNormal() if present; else infer from whether dispersion was supplied on prior_list
  if ("ddef" %in% names(prior_list)) {
    ddef <- prior_list$ddef
  } else {
    ddef <- is.null(prior_list$dispersion)
  }
  if (family$family %in% c("gaussian", "Gamma") && isTRUE(ddef)) {
    stop(paste0(
      "For gaussian() and Gamma() models, dNormal() requires an explicit dispersion ",
      "(e.g. dispersion = ps$dispersion from Prior_Setup()). ",
      "Omitted or NULL dispersion is not allowed"
    ))
  }
  
  if(family$family=="gaussian"){ 
    outlist<-.rNormalReg_cpp(n=n,y=y,x=x,mu=mu,P=P,offset=offset2,wt=wt,dispersion=dispersion,
                            ##                      famfunc=famfunc,f1=f1,
                            f2=f2,f3=f3,start=mu)
    class(outlist$fit)="lm"
    
  }
  else{
    if(is.null(dispersion)){dispersion2=1}
    else{dispersion2=dispersion}
    
    #  stop("Inputs to function above")
    outlist<-.rNormalGLM_cpp(n=n,y=y,x=x,mu=mu,P=P,offset=offset2,wt=wt,
                             dispersion=dispersion2,
                             ##famfunc=famfunc,f1=f1,
                             f2=f2,f3=f3,
                             start=start,family=family$family,link=family$link,Gridtype=Gridtype,
                             n_envopt = n_envopt,       # pass through
                             use_parallel = use_parallel,
                             use_opencl = use_opencl,
                             verbose = verbose)
    

    betastar=outlist$coef.mode  # Posterior mode from optim
    x=outlist$x
    y=outlist$y
    weights=outlist$prior.weights
    
    
    if(family$family=="quasipoisson"||family$family=="quasibinomial"){
      
      
      linkinv<-family$linkinv
      
      ## Compute dispersion and then rerun
      disp_temp=rep(0,n)
      m=length(y)  
      k=ncol(x)    
      res_temp=matrix(0,nrow=n,ncol=m)
      fit_temp=x%*%t(outlist$coefficients)
      for(l in 1:n){
        fit_temp[1:m,l]=linkinv(offset2+fit_temp[1:m,l])
        res_temp[l,1:m]=(y-fit_temp[1:m,l])
        disp_temp[l]=(1/(m-k))*sum(res_temp[l,1:m]^2*wt/fit_temp[1:m,l])
        
      }
      
      outlist<-.rNormalGLM_cpp(n=n,y=y,x=x,mu=mu,P=P,offset=offset2,
                               wt=wt,
                               dispersion=mean(disp_temp),
                               f2=f2,f3=f3,
                               start=start,family=family$family,link=family$link,Gridtype=Gridtype,
                               n_envopt = n_envopt,       # pass through
                               use_parallel = use_parallel,
                               use_opencl = use_opencl,
                               verbose = verbose)

            
      outlist$call <- match.call()  # overwrite with the rNormal_reg call
      outlist$dispersion=mean(disp_temp)
      
    }
    

    
    ## get influence info for original model
    outlist$fit=glmb.wfit(x,y,weights,offset=offset2,family=family,Bbar=mu,P,betastar)
    
    
  }
  
  
  colnames(outlist$coefficients)<-colnames(x)

  ## Build pfamily object so summary.rglmb etc. can detect prior type
  pl_disp <- if (!is.null(dispersion)) dispersion else outlist$dispersion
  pfamily_obj <- list(
    pfamily = "dNormal",
    prior_list = list(mu = as.numeric(mu), Sigma = Sigma, dispersion = pl_disp),
    pfun = rNormal_prior
  )
  attr(pfamily_obj, "Prior Type") <- "dNormal"
  class(pfamily_obj) <- "pfamily"
  outlist$pfamily <- pfamily_obj
  
  # include family in final list
  
  rglmb_df=as.data.frame(cbind(y,x))
  rglmb_f=DF2formula(rglmb_df)
  rglmb_mf=model.frame(rglmb_f,rglmb_df)
  
  outlist$family=family
  outlist$famfunc=famfunc
  outlist$call<-match.call()
  outlist$offset2<-offset2
  outlist$formula<-rglmb_f
  outlist$model<-rglmb_mf
  outlist$data<-rglmb_df
  
  class(outlist)<-c(outlist$class,c("rglmb","glmb","glm","lm"))
  outlist
  
}



# Helpers --------------------------------------------------------------------

################################## Utility functions used by the above  #################



# Helper: log(exp(b) - exp(a)) for a < b, in a numerically stable way
logdiffexp <- function(a, b) {
  # assumes a < b (strict), returns log(exp(b) - exp(a))
  # exp(b) - exp(a) = exp(b) * (1 - exp(a - b))
  # so log(...) = b + log(1 - exp(a - b))
  b + log(-expm1(a - b))
}



# .rNormalGLM_std_cpp --> rNormalGLM_std
# .rnorm_reg_cpp --> rNormalReg_cpp
# .rindep_norm_gamma_reg_cpp --> rIndepNormalGammaReg_cpp
# .rindep_norm_gamma_reg_std_cpp -->rIndepNormalGammaReg_std_cpp
# .rindep_norm_gamma_reg_std_parallel_cpp --> rIndepNormalGammaReg_std_parallel_cpp
# .rnnorm_reg_cpp --> .rNormalGLM_cpp
