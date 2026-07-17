#' Fitting Bayesian Linear Models
#'
#' \code{lmb} is used to fit Bayesian linear models, specified by giving a symbolic descriptions of the linear 
#' predictor and the prior distribution.
#' @name lmb
#' @aliases
#' lmb
#' print.lmb
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param pfamily a description of the prior distribution and associated constants to be used in the model.
#' For a single-response formula this should be a single \code{\link{pfamily}} object.
#' For a multi-response formula (e.g. \code{cbind(y1, y2) ~ x}) this must be a \code{list} of
#' \code{\link{pfamily}} objects with exactly one entry per response column; passing a single
#' \code{pfamily} object is an error.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain \code{NA}s.  The default is set by 
#' the \code{na.action} setting of \code{\link{options}}, and is \code{\link[stats]{na.fail}} 
#' if that is unset.  The \sQuote{factory-fresh} default is \code{stats{na.omit}}.  
#' Another possible value is \code{NULL}, no action.  Value \code{stats{na.exclude}} 
#' can be useful.
#' @param contrasts an optional list. See the \code{contrasts.arg} of 
#' \code{model.matrix.default}.
#' @param offset this can be used to specify an a priori known component to be 
#' included in the linear predictor during fitting. This should be \code{NULL} or a numeric 
#' vector or matrix of extents
#' matching those of the response.  One or more \code{offset} terms can be
#' included in the formula instead or as well, and if more than one are specified their 
#' sum is used.  See \code{model.offset}.
#' @param method the method to be used in fitting the classical model during a call to \code{\link{glm}}. The default method \code{glm.fit} 
#' uses iteratively reweighted least squares (IWLS): the alternative \code{"model.frame"} returns the model frame and does no fitting.
#' User-supplied fitting functions can be supplied either as a function or a character string naming a 
#' function, with a function which takes the same arguments as \code{glm.fit}. If specified as a character string it is looked up from within the \pkg{stats} namespace.
#' @param digits the number of significant digits to use when printing.
#' @inheritParams stats::lm
#' @inheritParams glmb
#' @return \code{lmb} returns an object of class \code{"lmb"}. The function \code{summary} (i.e., 
#' \code{\link{summary.glmb}}) can be used to obtain or print a summary of the results.  The generic accessor functions 
#' \code{\link{coefficients}}, \code{\link{fitted.values}}, \code{\link{residuals}}, and \code{\link{extractAIC}} can be used 
#' to extract various useful features of the value returned by \code{\link{lmb}}.
#' 
#' An object of class \code{"lmb"} is a list containing at least the following components:
#' \item{lm}{an object of class \code{"lm"} containing the output from a call to the function \code{\link{lm}}}
#' \item{coefficients}{a matrix of dimension \code{n} by \code{length(mu)} with one sample in each row}
#' \item{coef.means}{a vector of \code{length(mu)} with the estimated posterior mean coefficients}
#' \item{coef.mode}{a vector of \code{length(mu)} with the estimated posterior mode coefficients}
#' \item{dispersion}{Either a constant provided as part of the call, or a vector of length \code{n} with one sample in each row.}
#' \item{Prior}{A list with the priors specified for the model in question. Items in
#' list may vary based on the type of prior}
#' \item{residuals}{a matrix of dimension \code{n} by \code{length(y)} with one sample for the deviance residuals in each row}
#' \item{fitted.values}{a matrix of dimension \code{n} by \code{length(y)} with one sample for the fitted values in each row}
#' \item{linear.predictors}{an \code{n} by \code{length(y)} matrix with one sample for the linear fit on the link scale in each row}
#' \item{deviance}{an \code{n} by \code{1} matrix with one sample for the deviance in each row}
#' \item{pD}{An Estimate for the effective number of parameters}
#' \item{Dbar}{Expected value for minus twice the log-likelihood function}
#' \item{Dthetabar}{Value of minus twice the log-likelihood function evaluated at the mean value for the coefficients}   
#' \item{DIC}{Estimated Deviance Information criterion} 
#' \item{weights}{a vector of weights specified or implied by the model} 
#' \item{prior.weights}{a vector of weights specified or implied by the model} 
#' \item{y}{a vector of observations of length \code{m}.} 
#' \item{x}{a design matrix of dimension \code{m * p}} 
#' \item{model}{if requested (the default),the model frame} 
#' \item{call}{the matched call} 
#' \item{formula}{the formula supplied} 
#' \item{terms}{the \code{\link{terms}} object used} 
#' \item{data}{the \code{data argument}} 
#' \item{famfunc}{family functions used during estimation and post processing}
#' \item{iters}{an \code{n} by \code{1} matrix giving the number of candidates generated before acceptance for each sample.}
#' \item{contrasts}{(where relevant) the contrasts used.}
#' \item{xlevels}{(where relevant) a record of the levels of the factors used in fitting}
#' \item{pfamily}{the prior family specified}
#' \item{digits}{the number of significant digits to use when printing.}
#' In addition, non-empty fits will have (yet to be implemented) components \code{qr}, \code{R}
#' and \code{effects} relating to the final weighted linear fit for the posterior mode.  
#' Objects of class \code{"lmb"} are normally of class \code{c("lmb","glmb","glm","lm")},
#' that is inherit from classes \code{glmb}. \code{glm} and \code{lm} and well-designed
#' methods from those classed will be applied when appropriate.
#' 
#' @details
#' The function \code{lmb} is a Bayesian extension of the classical \code{\link[stats]{lm}} function. 
#' It retains the familiar formula interface and model setup used in \code{lm}, while introducing 
#' posterior simulation and prior specification via the \code{\link{pfamily}} argument. Internally, 
#' \code{lmb} calls \code{lm} to obtain the classical least squares fit, then generates independent 
#' draws from the posterior distribution using either multivariate normal simulation (for Gaussian priors) 
#' or accept-reject sampling via likelihood-subgradient envelopes \insertCite{Nygren2006}{glmbayes}.
#'
#' The symbolic formula interface follows \insertCite{WilkinsonRogers1973}{glmbayes}, and the overall design of \code{lm} was inspired by the S system 
#' \insertCite{Chambers1992}{glmbayes}. \code{lmb} comes with many of the same types of generic methods 
#' that are available to \code{lm} and \code{glm}, including \code{predict}, \code{residuals}, \code{extractAIC}, 
#' and \code{summary}. Many of these are inherited from \code{glmb}.
#'
#' Prior specification is handled via the \code{\link{pfamily}} argument, which defines the prior mean, 
#' covariance, and dispersion. The design of the \code{pfamily} family of functions was created by Kjell Nygren 
#' and is modeled on how \code{\link{glm}} uses \code{family} to specify the likelihood. A helper function, 
#' \code{\link{Prior_Setup}}, assists users in choosing prior parameters. It ships with sensible defaults but 
#' also allows full customization. All models support the \code{dNormal} prior; the Gaussian family also supports 
#' \code{dNormalGamma} and \code{dIndependent_Normal_Gamma}, which allow for more flexible prior structures 
#' including independent priors on variance components.
#'
#' Posterior draws are generated using the prior specification provided via \code{pfamily}. For Gaussian models 
#' with conjugate priors, draws are obtained directly from the posterior distribution using standard simulation 
#' procedures for multivariate normal densities \insertCite{Raiffa1961}{glmbayes}. For non-conjugate setups, the 
#' function uses envelope-based accept-reject sampling, where the \code{Gridtype} parameter controls the granularity 
#' of the envelope construction. The number of candidates generated before acceptance is returned in the 
#' \code{iters} component.
#'
#' The output includes both the classical \code{lm} fit and Bayesian diagnostics such as the Deviance 
#' Information Criterion (DIC), effective number of parameters (pD), and posterior summaries. 
#' The DIC, introduced by \insertCite{Spiegelhalter2002}{glmbayes}, provides a Bayesian analog to AIC 
#' by balancing model fit and complexity using posterior expectations. This dual structure allows users 
#' to compare classical and Bayesian fits side-by-side, and to leverage familiar modeling workflows 
#' while gaining access to richer inferential tools.
#'
#' The \code{\link{lmb}} function is a specialized version of \code{\link{glmb}} for Gaussian models, 
#' and does not require a \code{family} argument. For conjugate models, it uses standard simulation methods for 
#' posterior draws, avoiding the need for envelope construction or subgradient sampling. 
#' Like \code{glmb}, it returns objects compatible with many standard methods from \code{lm} and \code{glm}, 
#' including \code{\link{extractAIC}}, \code{\link{fitted.values}}, and \code{\link{residuals}}.
#'
#' For more minimalistic workflows, \code{\link{rlmb}} and \code{\link{rglmb}} offer stripped-down interfaces 
#' for posterior sampling without the overhead of full model objects. \code{rlmb} is called from within 
#' \code{lmb}. The functions \code{rlmb} might be useful in Gibbs sampling or simulation-heavy contexts.
#'  
#' 
#' @author The \R implementation of \code{lmb} has been written by Kjell Nygren and
#' was built to be a Bayesian version of the \code{lm} function and hence tries
#' to mirror the features of the \code{lm} function to the greatest extent possible while also taking advantage
#' of some of the method functions developed for the \code{glmb} function. For details
#' on the author(s) for the \code{lm} function see the documentation for \code{\link[stats]{lm}}.
#'     
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @family modelfuns
#' @seealso The classical modeling functions \code{\link[stats]{lm}} and \code{\link[stats]{glm}}.
#'
#' \code{\link{glmb}}, \code{\link{rglmb}}, \code{\link{rlmb}}
#'   for related Bayesian GLM/linear interfaces;
#' \code{\link{EnvelopeBuild}} for envelope construction when accept--reject sampling is used.
#' 
#' \code{\link{pfamily}} for documentation of pfamily functions used to specify priors.
#' 
#' \code{\link{Prior_Setup}}, \code{\link{Prior_Check}} for functions used to initialize and to check priors,  
#'
#' Further reading: \insertCite{Nygren2006}{glmbayes};
#' \insertCite{glmbayesSimmethods,glmbayesChapterA08}{glmbayes};
#' independent Normal--Gamma sampler: \insertCite{glmbayesIndNormGammaVignette}{glmbayes}.
#'
#' \code{\link{summary.glmb}}, \code{\link{predict.glmb}}, \code{\link{simulate.glmb}}, 
#' \code{\link{extractAIC.glmb}}, \code{\link{dummy.coef.glmb}} and methods(class="glmb") for methods 
#' inherited from class \code{glmb} and the methods and generic functions for classes \code{glm} and 
#' \code{lm} from which class \code{lmb} also inherits.
#'
#' @example inst/examples/Ex_lmb.R
#' 
#' @export
lmb <- function(
    formula,
    pfamily,
    n = 1000,
    data,
    subset,
    weights,
    na.action,
    method = "qr",
    model = TRUE,
    x = TRUE,
    y = TRUE,
    qr = TRUE,
    singular.ok = TRUE,
    contrasts = NULL,
    offset,
    Gridtype = 2,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE,
    ...
) {
  mc <- match.call(expand.dots = FALSE)
  call <- match.call()
  if (missing(data)) {
    data <- environment(formula)
  }

  meta_args <- list(
    formula = formula,
    data = data,
    contrasts = contrasts,
    has_subset = !missing(subset),
    has_weights = !missing(weights),
    has_na.action = !missing(na.action),
    has_offset = !missing(offset)
  )
  if (!missing(subset)) meta_args$subset <- subset
  if (!missing(weights)) meta_args$weights <- weights
  if (!missing(na.action)) meta_args$na.action <- na.action
  if (!missing(offset)) meta_args$offset <- offset
  meta <- do.call(.mlmb_formula_meta, meta_args)

  lmb_args <- list(
    n = n,
    data = data,
    method = method,
    model = model,
    x = x,
    y = y,
    qr = qr,
    singular.ok = singular.ok,
    contrasts = contrasts,
    Gridtype = Gridtype,
    n_envopt = n_envopt,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = verbose
  )
  if (!missing(subset)) lmb_args$subset <- subset
  if (!missing(weights)) lmb_args$weights <- weights
  if (!missing(na.action)) lmb_args$na.action <- na.action
  if (!missing(offset)) lmb_args$offset <- offset
  if (length(list(...))) {
    lmb_args <- c(lmb_args, list(...))
  }

  if (meta$l1 == 1L) {
    fit <- do.call(
      .uni_lmb,
      c(list(formula = formula, pfamily = pfamily), lmb_args)
    )
    fit$call <- .mlmb_lmb_display_call(mc, formula)
    return(fit)
  }

  if (!is.list(pfamily) || inherits(pfamily, "pfamily") ||
      (!is.null(pfamily$pfamily) && !is.null(pfamily$prior_list))) {
    stop(
      "lmb(): 'pfamily' must be a list of ", meta$l1, " pfamily objects ",
      "(one per response column) when y has multiple columns.\n",
      "  Provide: pfamily = list(dNormal(...), dNormal(...), ...)",
      call. = FALSE
    )
  }
  if (length(pfamily) != meta$l1) {
    stop(
      "lmb(): length(pfamily) = ", length(pfamily),
      " but y has ", meta$l1, " columns. ",
      "Supply one pfamily per response column.",
      call. = FALSE
    )
  }
  pfamily_list <- pfamily
  pfamily_lists <- .mrglmb_normalize_pfamily_lists(
    pfamily_list,
    meta$l1,
    meta$p,
    .validate_pfamily_for_rlmb
  )

  block_results <- vector("list", meta$l1)
  for (j in seq_len(meta$l1)) {
    f_j <- stats::reformulate(meta$termlabels, response = meta$coef_names[j])
    fit_j <- do.call(
      .uni_lmb,
      c(list(formula = f_j, pfamily = pfamily_lists[[j]]), lmb_args)
    )
    fit_j$call <- .mlmb_lmb_display_call(mc, f_j)
    block_results[[j]] <- fit_j
  }

  .mlmb_assemble(
    block_results = block_results,
    coef_names = meta$coef_names,
    call = call,
    formula = formula,
    l1 = meta$l1,
    p = meta$p,
    pred_names = meta$pred_names,
    pfamily_lists = pfamily_lists
  )
}


## Portions of the model-frame/model-matrix setup, fitted-object structure,
## and method conventions follow or adapt stats::lm() and related stats
## methods. See inst/COPYRIGHTS for derived-code notices.
if (FALSE) {
lmb <- function (
    formula,
    pfamily,
    n = 1000,
    data,
    subset,
    weights,
    na.action,
    method = "qr",
    model = TRUE,
    x = TRUE,
    y = TRUE,
    qr = TRUE,
    singular.ok = TRUE,
    contrasts = NULL,
    offset,
    Gridtype = 2,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE,
    ...
){
  ret.x <- x
  ret.y <- y
  cl <- match.call(expand.dots = FALSE)
  if (length(cl) >= 1L && is.function(cl[[1L]])) {
    cl[[1L]] <- as.name("lmb")
  }
  mf <- cl
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  ## need stats:: for non-standard evaluation
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (method == "model.frame")
    return(mf)
  else if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'", method),
            domain = NA)
  mt <- attr(mf, "terms") # allow model.frame to update it
  y <- model.response(mf, "numeric")
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }
  
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (is.matrix(y))
      matrix(, 0, 3) else numeric(), residuals = y,
      fitted.values = 0 * y, weights = w, rank = 0L,
      df.residual = if (!is.null(w)) sum(w != 0) else
        if (is.matrix(y)) nrow(y) else length(y))
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    z <- if (is.null(w)) lm.fit(x, y, offset = offset,
                                singular.ok = singular.ok, ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, ...)
  }
  class(z) <- c(if (is.matrix(y)) "mlm", "lm")
  
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- x
  if (ret.y)
    z$y <- y
  if (!qr) z$qr <- NULL
  
  if (!is.null(x)) {
    z$assign <- attr(x, "assign")
  }
  
  ######   End of lm function
  # Verify inputs and Initialize
  
  ## Pull in information from the pfamily
  prior_list <- pfamily$prior_list
  y <- z$y
  x <- z$x
  b <- z$coefficients

  if (pfamily$pfamily != "dGamma") {
    mu <- as.matrix(as.vector(prior_list$mu))
    Sigma <- as.matrix(prior_list$Sigma)
    dispersion <- prior_list$dispersion

    R <- chol(Sigma)
    P <- chol2inv(R)
    P <- 0.5 * (P + t(P))
  }
  
  if (is.null(z$weights)) wtin <- rep(1, length(y))
  else wtin <- z$weights
  
  ## normalize n_envopt (mirror glmb)
  if (is.null(n_envopt)) n_envopt <- n
  n_envopt <- as.integer(n_envopt)
  
  sim <- rlmb(
    n       = n,
    y       = y,
    x       = x,
    pfamily = pfamily,
    offset  = offset,
    weights = wtin,
    Gridtype   = Gridtype,
    n_envopt   = n_envopt,
    use_parallel = use_parallel,
    use_opencl  = use_opencl,
    verbose     = verbose
  )
  
  if (pfamily$pfamily == "dIndependent_Normal_Gamma") {
    if (!is.null(sim$sim_bounds)) {
      pfamily$prior_list$disp_lower <- sim$sim_bounds$low
      pfamily$prior_list$disp_upper <- sim$sim_bounds$upp
    } else {
      cat("No simbounds returned in sim.\n")
    }
  }
  
  dispersion2 <- sim$dispersion
  famfunc <- sim$famfunc

  if (pfamily$pfamily == "dGamma") {
    Prior <- list(shape = prior_list$shape, rate = prior_list$rate)
  } else {
    Prior <- list(mean = as.numeric(mu), Variance = Sigma)
    names(Prior$mean) <- colnames(z$x)
    colnames(Prior$Variance) <- colnames(z$x)
    rownames(Prior$Variance) <- colnames(z$x)
  }
  
  if (!is.null(offset)) {
    
    if (length(dispersion2) == 1) {
      
      # Scale weights by dispersion for consistent deviance computation
      wt_scaled <- wtin / dispersion2
      
      DICinfo <- DIC_Info(sim$coefficients, y = y, x = x, alpha = offset,
                          f1 = famfunc$f1, f4 = famfunc$f4,
                          wt = wt_scaled, dispersion = 1)
    }
    
    if (length(dispersion2) > 1) {
      
      DICinfo <- DIC_Info(sim$coefficients, y = y, x = x, alpha = 0,
                          f1 = famfunc$f1, f4 = famfunc$f4,
                          wt = wtin, dispersion = dispersion2)
    }
    
    linear.predictors <- t(offset + x %*% t(sim$coefficients))
    fitted.values <- linear.predictors
    
  }
  
  if (is.null(offset)) {
    
    if (length(dispersion2) == 1) {
      
      # Scale weights by dispersion for consistent deviance computation
      wt_scaled <- wtin / dispersion2
      
      DICinfo <- DIC_Info(sim$coefficients, y = y, x = x, alpha = 0,
                          f1 = famfunc$f1, f4 = famfunc$f4,
                          wt = wt_scaled, dispersion = 1)
    }
    
    if (length(dispersion2) > 1) {
      
      DICinfo <- DIC_Info(sim$coefficients, y = y, x = x, alpha = 0,
                          f1 = famfunc$f1, f4 = famfunc$f4,
                          wt = wtin, dispersion = dispersion2)
    }
    
    linear.predictors <- t(x %*% t(sim$coefficients))
    fitted.values <- linear.predictors
    
  }
  
  # For dGamma, coefficients have 1 row; replicate to n rows for consistent structure
  if (nrow(sim$coefficients) == 1L) {
    fitted.values <- matrix(rep(fitted.values, n), nrow = n, byrow = TRUE)
    linear.predictors <- matrix(rep(linear.predictors, n), nrow = n, byrow = TRUE)
  }
  
  residuals <- fitted.values
  
  for (i in 1:n) {
    residuals[i, 1:length(y)] <- y - residuals[i, 1:length(y)]
  }
  
  outlist <- list(
    lm = z,
    coefficients = sim$coefficients,
    coef.means = colMeans(sim$coefficients),
    coef.mode = sim$coef.mode,
    dispersion = dispersion2,
    residuals = residuals,
    Prior = Prior,
    fitted.values = fitted.values,
    linear.predictors = linear.predictors,
    deviance = DICinfo$Deviance,
    pD = DICinfo$pD,
    Dbar = DICinfo$Dbar,
    Dthetabar = DICinfo$Dthetabar,
    DIC = DICinfo$DIC,
    prior.weights = wtin,
    weights = wtin,
    offset = offset,
    y = z$y,
    x = z$x,
    model = z$model,
    call = z$call,
    formula = z$formula,
    terms = z$terms,
    data = mf,
    fit = sim$fit,
    famfunc = famfunc,
    iters = sim$iters,
    contrasts = z$contrasts,
    xlevels = z$xlevels,
    pfamily = pfamily,
    simfun_call = sim$simfun_call,
    simfun_args = sim$simfun_args
  )
  
  outlist$call <- cl
  
  if (pfamily$pfamily == "dGamma") {
    class(outlist) <- c("rGamma_reg", outlist$class, "lmb", "glmb", "glm", "lm")
  } else {
    class(outlist) <- c(outlist$class, "lmb", "glmb", "glm", "lm")
  }
  outlist
}
} ## end if(FALSE) -- old lmb

#' @rdname lmb
#' @method print lmb
#' @export


print.lmb<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  
  cat("\nCall:  \n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Posterior Mean Coefficients")
    cat(":\n")
    print.default(format(x$coef.means, digits = digits), 
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
}

#' @keywords internal
.uni_lmb <- function (
    formula,
    pfamily,
    n = 1000,
    data,
    subset,
    weights,
    na.action,
    method = "qr",
    model = TRUE,
    x = TRUE,
    y = TRUE,
    qr = TRUE,
    singular.ok = TRUE,
    contrasts = NULL,
    offset,
    Gridtype = 2,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE,
    ...
){
  ret.x <- x
  ret.y <- y
  cl <- match.call(expand.dots = FALSE)
  if (length(cl) >= 1L && is.function(cl[[1L]])) {
    cl[[1L]] <- as.name("lmb")
  }
  mf <- cl
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  ## need stats:: for non-standard evaluation
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (method == "model.frame")
    return(mf)
  else if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'", method),
            domain = NA)
  mt <- attr(mf, "terms") # allow model.frame to update it
  y <- model.response(mf, "numeric")
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }
  
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (is.matrix(y))
      matrix(, 0, 3) else numeric(), residuals = y,
      fitted.values = 0 * y, weights = w, rank = 0L,
      df.residual = if (!is.null(w)) sum(w != 0) else
        if (is.matrix(y)) nrow(y) else length(y))
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    z <- if (is.null(w)) lm.fit(x, y, offset = offset,
                                singular.ok = singular.ok, ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, ...)
  }
  class(z) <- c(if (is.matrix(y)) "mlm", "lm")
  
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- x
  if (ret.y)
    z$y <- y
  if (!qr) z$qr <- NULL
  
  if (!is.null(x)) {
    z$assign <- attr(x, "assign")
  }
  
  ######   End of lm function
  # Verify inputs and Initialize
  
  ## Pull in information from the pfamily
  prior_list <- pfamily$prior_list
  y <- z$y
  x <- z$x
  b <- z$coefficients

  if (pfamily$pfamily != "dGamma") {
    mu <- as.matrix(as.vector(prior_list$mu))
    Sigma <- as.matrix(prior_list$Sigma)
    dispersion <- prior_list$dispersion

    R <- chol(Sigma)
    P <- chol2inv(R)
    P <- 0.5 * (P + t(P))
  }
  
  if (is.null(z$weights)) wtin <- rep(1, length(y))
  else wtin <- z$weights
  
  ## normalize n_envopt (mirror glmb)
  if (is.null(n_envopt)) n_envopt <- n
  n_envopt <- as.integer(n_envopt)
  
  sim <- rlmb(
    n       = n,
    y       = y,
    x       = x,
    pfamily = pfamily,
    offset  = offset,
    weights = wtin,
    Gridtype   = Gridtype,
    n_envopt   = n_envopt,
    use_parallel = use_parallel,
    use_opencl  = use_opencl,
    verbose     = verbose
  )
  
  if (pfamily$pfamily == "dIndependent_Normal_Gamma") {
    if (!is.null(sim$sim_bounds)) {
      pfamily$prior_list$disp_lower <- sim$sim_bounds$low
      pfamily$prior_list$disp_upper <- sim$sim_bounds$upp
    } else {
      cat("No simbounds returned in sim.\n")
    }
  }
  
  dispersion2 <- sim$dispersion
  famfunc <- sim$famfunc

  if (pfamily$pfamily == "dGamma") {
    Prior <- list(shape = prior_list$shape, rate = prior_list$rate)
  } else {
    Prior <- list(mean = as.numeric(mu), Variance = Sigma)
    names(Prior$mean) <- colnames(z$x)
    colnames(Prior$Variance) <- colnames(z$x)
    rownames(Prior$Variance) <- colnames(z$x)
  }
  
  if (!is.null(offset)) {
    
    if (length(dispersion2) == 1) {
      
      # Scale weights by dispersion for consistent deviance computation
      wt_scaled <- wtin / dispersion2
      
      DICinfo <- DIC_Info(sim$coefficients, y = y, x = x, alpha = offset,
                          f1 = famfunc$f1, f4 = famfunc$f4,
                          wt = wt_scaled, dispersion = 1)
    }
    
    if (length(dispersion2) > 1) {
      
      DICinfo <- DIC_Info(sim$coefficients, y = y, x = x, alpha = 0,
                          f1 = famfunc$f1, f4 = famfunc$f4,
                          wt = wtin, dispersion = dispersion2)
    }
    
    linear.predictors <- t(offset + x %*% t(sim$coefficients))
    fitted.values <- linear.predictors
    
  }
  
  if (is.null(offset)) {
    
    if (length(dispersion2) == 1) {
      
      # Scale weights by dispersion for consistent deviance computation
      wt_scaled <- wtin / dispersion2
      
      DICinfo <- DIC_Info(sim$coefficients, y = y, x = x, alpha = 0,
                          f1 = famfunc$f1, f4 = famfunc$f4,
                          wt = wt_scaled, dispersion = 1)
    }
    
    if (length(dispersion2) > 1) {
      
      DICinfo <- DIC_Info(sim$coefficients, y = y, x = x, alpha = 0,
                          f1 = famfunc$f1, f4 = famfunc$f4,
                          wt = wtin, dispersion = dispersion2)
    }
    
    linear.predictors <- t(x %*% t(sim$coefficients))
    fitted.values <- linear.predictors
    
  }
  
  # For dGamma, coefficients have 1 row; replicate to n rows for consistent structure
  if (nrow(sim$coefficients) == 1L) {
    fitted.values <- matrix(rep(fitted.values, n), nrow = n, byrow = TRUE)
    linear.predictors <- matrix(rep(linear.predictors, n), nrow = n, byrow = TRUE)
  }
  
  residuals <- fitted.values
  
  for (i in 1:n) {
    residuals[i, 1:length(y)] <- y - residuals[i, 1:length(y)]
  }
  
  outlist <- list(
    lm = z,
    coefficients = sim$coefficients,
    coef.means = colMeans(sim$coefficients),
    coef.mode = sim$coef.mode,
    dispersion = dispersion2,
    residuals = residuals,
    Prior = Prior,
    fitted.values = fitted.values,
    linear.predictors = linear.predictors,
    deviance = DICinfo$Deviance,
    pD = DICinfo$pD,
    Dbar = DICinfo$Dbar,
    Dthetabar = DICinfo$Dthetabar,
    DIC = DICinfo$DIC,
    prior.weights = wtin,
    weights = wtin,
    offset = offset,
    y = z$y,
    x = z$x,
    model = z$model,
    call = z$call,
    formula = z$formula,
    terms = z$terms,
    data = mf,
    fit = sim$fit,
    famfunc = famfunc,
    iters = sim$iters,
    contrasts = z$contrasts,
    xlevels = z$xlevels,
    pfamily = pfamily,
    simfun_call = sim$simfun_call,
    simfun_args = sim$simfun_args
  )
  
  outlist$call <- cl
  
  if (pfamily$pfamily == "dGamma") {
    class(outlist) <- c("rGamma_reg", outlist$class, "lmb", "glmb", "glm", "lm")
  } else {
    class(outlist) <- c(outlist$class, "lmb", "glmb", "glm", "lm")
  }
  outlist
}
## Multi-response Bayesian linear models (\code{lmb})
##
## @description
## Fits one \code{\link{lmb}} model per column of a multivariate response
## (formula left-hand side with \code{cbind(...)}), sharing the same predictors
## on the right-hand side. Returns a named list of \code{"lmb"} objects with
## class \code{"mlmb"}.
##
## @details
## This is the formula / \code{data} interface counterpart to
## \code{\link{multi_rlmb}} (matrix \code{y}, \code{x}). Each response column
## uses its own \code{pfamily_list[[j]]}. Use \code{\link{multi_prior_setup}}
## to build aligned priors, then \code{\link{summary.mlmb}} or
## \code{\link{print.mlmb}} for output styled like \code{\link[stats]{summary.mlm}}
## (cf.\ \code{\link[stats]{lm}} with a matrix response via \code{cbind(...)}).
##
## @param formula A \code{\link{formula}} with a matrix response on the left-hand
##   side (typically \code{cbind(...)}).
## @param pfamily_list Named or unnamed list of length equal to the number of
##   response columns; each element is a \code{\link{pfamily}} object for
##   \code{\link{lmb}}.
## @inheritParams lmb
## @inheritParams multi_rlmb
## @return A named list of class \code{"mlmb"}. Element \code{j} is an
##   \code{"lmb"} fit for response \code{j}. Attributes include \code{call},
##   \code{formula}, \code{coef_names}, \code{pred_names}, and
##   \code{pfamily_lists}.
## @seealso \code{\link{lmb}}, \code{\link{multi_rlmb}}, \code{\link{multi_prior_setup}},
##   \code{\link{summary.mlmb}}, \code{\link{print.mlmb}},
##   \code{\link[stats]{lm}} with \code{cbind} responses.
## @family modelfuns
## @example inst/examples/Ex_multi_lmb.R
## @export
#' @keywords internal
#' Build a short \code{lmb()} call for printing (omits \code{pfamily}).
.mlmb_lmb_display_call <- function(mc_multi, formula_j) {
  cl <- call("lmb", formula = formula_j)
  pass <- c(
    "n", "data", "subset", "weights", "na.action", "offset",
    "use_parallel", "use_opencl", "verbose", "Gridtype", "n_envopt"
  )
  for (nm in pass) {
    if (nm %in% names(mc_multi)) {
      cl[[nm]] <- mc_multi[[nm]]
    }
  }
  cl
}

#' @keywords internal
.mlmb_formula_meta <- function(
    formula,
    data,
    subset = NULL,
    weights = NULL,
    na.action = NULL,
    offset = NULL,
    contrasts = NULL,
    has_subset = FALSE,
    has_weights = FALSE,
    has_na.action = FALSE,
    has_offset = FALSE
) {
  mf_args <- list(
    formula = formula,
    data = data,
    drop.unused.levels = TRUE
  )
  if (has_subset) mf_args$subset <- subset
  if (has_weights) mf_args$weights <- weights
  if (has_na.action) mf_args$na.action <- na.action
  if (has_offset) mf_args$offset <- offset
  if (!is.null(contrasts)) mf_args$contrasts <- contrasts
  mf <- do.call(stats::model.frame, mf_args)
  mt <- attr(mf, "terms")
  Y <- as.matrix(stats::model.response(mf, "any"))
  l1 <- ncol(Y)
  if (l1 < 1L) {
    stop(
      "formula must specify at least one response column (e.g. cbind(...)).",
      call. = FALSE
    )
  }
  coef_names <- colnames(Y)
  if (is.null(coef_names) || length(coef_names) != l1) {
    coef_names <- paste0("Y", seq_len(l1))
  }

  x_mat <- stats::model.matrix(mt, mf, contrasts)
  p <- ncol(x_mat)
  pred_names <- colnames(x_mat)
  if (is.null(pred_names) || length(pred_names) != p) {
    pred_names <- paste0("X", seq_len(p))
  }

  list(
    coef_names = coef_names,
    termlabels = attr(mt, "term.labels"),
    l1 = l1,
    p = p,
    pred_names = pred_names
  )
}

#' @keywords internal
.mlmb_assemble <- function(
    block_results,
    coef_names,
    call,
    formula,
    l1,
    p,
    pred_names,
    pfamily_lists = NULL
) {
  outlist <- setNames(block_results, coef_names)
  attr(outlist, "call") <- call
  attr(outlist, "formula") <- formula
  attr(outlist, "l1") <- l1
  attr(outlist, "p") <- p
  attr(outlist, "coef_names") <- coef_names
  attr(outlist, "pred_names") <- pred_names
  if (!is.null(pfamily_lists)) {
    attr(outlist, "pfamily_lists") <- pfamily_lists
  }
  class(outlist) <- "mlmb"
  outlist
}

#' @keywords internal
.mlmb_coef_means_matrix <- function(object) {
  nm <- names(object)
  if (length(nm) < 1L) {
    return(NULL)
  }
  cm <- do.call(cbind, lapply(object, function(fit) fit$coef.means))
  rn <- names(object[[1L]]$coef.means)
  if (is.null(rn)) {
    rn <- colnames(object[[1L]]$x)
  }
  if (!is.null(rn) && nrow(cm) == length(rn)) {
    rownames(cm) <- rn
  }
  colnames(cm) <- nm
  cm
}

#' @keywords internal
.mlmb_dic_table <- function(object) {
  nm <- names(object)
  if (length(nm) < 1L) {
    return(NULL)
  }
  pD <- vapply(object, function(fit) fit$pD, numeric(1))
  dic <- vapply(object, function(fit) fit$DIC, numeric(1))
  cbind(pD = pD, DIC = dic)
}

#' Prior setup for multiple Gaussian responses
#'
#' @inheritParams Prior_Setup
#' @return A named list of class \code{"multi_PriorSetup"}. Each element is a
#'   \code{\link{Prior_Setup}} result for one column of the response (names from
#'   \code{colnames(y)} or \code{Y1}, \code{Y2}, \ldots).
#' @family prior
#' @export
multi_prior_setup <- function(
    formula,
    family = gaussian(),
    data = NULL,
    weights = NULL,
    subset = NULL,
    na.action = na.fail,
    offset = NULL,
    contrasts = NULL,
    pwt = NULL,
    pwt_default_low = 0.01,
    pwt_default_high = 0.05,
    n_prior = NULL,
    sd = NULL,
    dispersion = NULL,
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model"),
    mu = NULL,
    k = 1,
    ...
) {
  call <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family) || family$family != "gaussian") {
    stop(
      "multi_prior_setup() currently supports family = gaussian() only.",
      call. = FALSE
    )
  }

  if (missing(data)) {
    data <- environment(formula)
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(
    c("formula", "data", "subset", "weights", "na.action", "offset"),
    names(mf),
    0L
  )
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  Y <- as.matrix(model.response(mf, "any"))
  l1 <- ncol(Y)
  if (l1 < 1L) {
    stop("formula must specify at least one response column.", call. = FALSE)
  }
  coef_names <- colnames(Y)
  if (is.null(coef_names) || length(coef_names) != l1) {
    coef_names <- paste0("Y", seq_len(l1))
  }

  termlabels <- attr(mt, "term.labels")
  ps_args <- list(
    family = gaussian(),
    data = data,
    weights = weights,
    subset = subset,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    pwt = pwt,
    pwt_default_low = pwt_default_low,
    pwt_default_high = pwt_default_high,
    n_prior = n_prior,
    sd = sd,
    dispersion = dispersion,
    intercept_source = intercept_source,
    effects_source = effects_source,
    mu = mu,
    k = k
  )

  setups <- setNames(vector("list", l1), coef_names)
  for (j in seq_len(l1)) {
    f_j <- stats::reformulate(termlabels, response = coef_names[j])
    setups[[j]] <- do.call(
      Prior_Setup,
      c(list(formula = f_j), ps_args, list(...))
    )
  }

  attr(setups, "call") <- call
  attr(setups, "formula") <- formula
  class(setups) <- c("multi_PriorSetup", "list")
  setups
}

#' @keywords internal
.mrglmb_normalize_pfamily_lists <- function(pfamily_list, l1, p, validate_fn) {
  if (!is.list(pfamily_list)) {
    stop(
      "pfamily_list must be a list of length ncol(y) of per-column pfamily objects.",
      call. = FALSE
    )
  }
  if (!is.null(pfamily_list$pfamily) && !is.null(pfamily_list$prior_list)) {
    stop(
      "pfamily_list must be a list of pfamily objects (one per column of y), ",
      "not a single pfamily with components pfamily and prior_list.",
      call. = FALSE
    )
  }
  if (length(pfamily_list) != l1) {
    stop("length(pfamily_list) must equal ncol(y) = ", l1, ".", call. = FALSE)
  }
  lapply(seq_len(l1), function(j) {
    validate_fn(pfamily_list[[j]], j = j, p = p)
  })
}

#' @keywords internal
.validate_pfamily_for_rlmb <- function(pl, j, p) {
  if (!inherits(pl, "pfamily")) {
    stop("pfamily_list[[", j, "]] must inherit from class \"pfamily\".", call. = FALSE)
  }
  if (is.null(pl$pfamily) || is.null(pl$prior_list) || is.null(pl$simfun)) {
    stop(
      "pfamily_list[[", j, "]] must contain 'pfamily', 'prior_list', and 'simfun'.",
      call. = FALSE
    )
  }
  mu <- pl$prior_list$mu
  if (is.null(mu)) {
    stop("pfamily_list[[", j, "]]$prior_list must contain 'mu'.", call. = FALSE)
  }
  mu <- as.numeric(mu)
  if (length(mu) != p) {
    stop(
      "pfamily_list[[", j, "]]$prior_list$mu must have length ncol(x) = ", p, ".",
      call. = FALSE
    )
  }
  pl
}