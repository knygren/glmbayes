## Phase 3: formula layer in glmbayes; sampling via glmbayesCore::rglmb().
## Portions follow stats::glm(); see inst/COPYRIGHTS.

#' Fitting Bayesian Generalized Linear Models
#'
#' \code{glmb} is used to fit Bayesian generalized linear models, specified by giving a symbolic descriptions of 
#' the linear predictor,  the error distribution, and the prior distribution.
#' @name glmb
#' @aliases
#' glmb
#' print.glmb
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param pfamily a description of the prior distribution and associated constants to be used in the model. This
#' should be a pfamily function (see \code{\link{pfamily}} for details of pfamily functions).
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain \code{NA}s.  The default is set by 
#' the \code{na.action} setting of \code{\link{options}}, and is \code{\link[stats]{na.fail}} 
#' if that is unset.  The \sQuote{factory-fresh} default is \code{stats{na.omit}}.  
#' Another possible value is \code{NULL}, no action.  Value \code{stats{na.exclude}} 
#' can be useful.
#' @param offset this can be used to specify an \emph{a priori} known component to be included in the linear 
#' predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of 
#' cases.  One or more \code{\link[stats]{offset}} terms can be included in the formula instead or as well, and 
#' if more than one is specified their sum is used.  See documentation for \code{model.offset} at 
#' \code{\link[stats]{model.extract}}.
#' @param Gridtype an optional argument specifying the method used to determine the number of tangent points used to construct the enveloping function.
#' @param n_envopt Effective sample size passed to EnvelopeOpt for grid
#'   construction. Defaults to match `n`. Larger values encourage tighter
#'   envelopes.
#' @param digits the number of significant digits to use when printing.
#' @inheritParams   stats::glm 
#' @param use_parallel Logical. Whether to use parallel processing during simulation.
#' @param use_opencl Logical. Whether to use OpenCL acceleration during Envelope construction.
#' @param verbose Logical. Whether to print progress messages.
#' @return \code{glmb} returns an object of class \code{"glmb"}. The function \code{summary} (i.e., 
#' \code{\link{summary.glmb}}) can be used to obtain or print a summary of the results.  The generic accessor functions 
#' \code{\link{coefficients}}, \code{\link{fitted.values}}, \code{\link{residuals}}, and \code{\link{extractAIC}} can be used 
#' to extract various useful features of the value returned by \code{\link{glmb}}.
#' 
#' An object of class \code{"glmb"} is a list containing at least the following components:
#' \item{glm}{an object of class \code{"glm"} containing the output from a call to the function \code{\link{glm}}}
#' \item{coefficients}{a matrix of dimension \code{n} by \code{length(mu)} with one sample in each row}
#' \item{coef.means}{a vector of \code{length(mu)} with the estimated posterior mean coefficients}
#' \item{coef.mode}{a vector of \code{length(mu)} with the estimated posterior mode coefficients}
#' \item{dispersion}{Either a constant provided as part of the call, or a vector of length \code{n} with one sample in each row.}
#' \item{Prior}{A list with the priors specified for the model in question. Items in
#' list may vary based on the type of prior}
#' \item{fitted.values}{a matrix of dimension \code{n} by \code{length(y)} with one sample for the mean fitted values in each row}
#' \item{family}{the \code{\link{family}} object used.}
#' \item{linear.predictors}{an \code{n} by \code{length(y)} matrix with one sample for the linear fit on the link scale in each row}
#' \item{deviance}{an \code{n} by \code{1} matrix with one sample for the deviance in each row}
#' \item{pD}{An Estimate for the effective number of parameters}
#' \item{Dbar}{Expected value for minus twice the log-likelihood function}
#' \item{Dthetabar}{Value of minus twice the log-likelihood function evaluated at the mean value for the coefficients}   
#' \item{DIC}{Estimated Deviance Information criterion} 
#' \item{prior.weights}{a vector of weights specified or implied by the model} 
#' \item{y}{a vector with the dependent variable} 
#' \item{x}{a matrix with the implied design matrix for the model} 
#' \item{model}{if requested (the default),the model frame} 
#' \item{call}{the matched call} 
#' \item{formula}{the formula supplie} 
#' \item{terms}{the \code{\link{terms}} object used} 
#' \item{data}{the \code{data argument}} 
#' \item{famfunc}{Family functions used during estimation process}
#' \item{iters}{an \code{n} by \code{1} matrix giving the number of candidates generated before acceptance for each sample.}
#' \item{contrasts}{(where relevant) the contrasts used.}
#' \item{xlevels}{(where relevant) a record of the levels of the factors used in fitting}
#' \item{digits}{the number of significant digits to use when printing.}
#' In addition, non-empty fits will have (yet to be implemented) components \code{qr}, \code{R}
#' and \code{effects} relating to the final weighted linear fit for the posterior mode.  
#' Objects of class \code{"glmb"} are normall of class \code{c("glmb","glm","lm")},
#' that is inherit from classes \code{glm} and \code{lm} and well-designed
#' methods from those classed will be applied when appropriate.
#' 
#' If a \link{binomial} \code{glmb} model was specified by giving a two-column 
#' response, the weights returned by \code{prior.weights} are the total number of
#' cases (factored by the supplied case weights) and the component of \code{y}
#' of the result is the proportion of successes.
#' 
#' @details
#' The function \code{glmb} is a Bayesian version of the classical
#' \code{\link{glm}} function. The original R implementation of
#' \code{glm} was written by Simon Davies (under Ross Ihaka at the
#' University of Auckland) and has since been extensively rewritten by
#' members of the R Core Team; its design was inspired by the S
#' function described in \insertCite{Hastie1992}{glmbayes}, which in turn relies on the 
#' formula framework described in \insertCite{WilkinsonRogers1973}{glmbayes}.
#'
#' Setup (including the use of formulas and families) mirrors that of \code{\link{glm}} but adds a 
#' required \code{pfamily} argument to specify the prior distribution. The design of the \code{pfamily} family of 
#' functions was created by Kjell Nygren and is modeled on how \code{glm} 
#' uses \code{family} to specify the likelihood.
#'
#' For any implemented combination of family, link, and \code{pfamily}, 
#' \code{glmb} generates independent draws from the posterior density-
#' no MCMC chains are required. Results can be printed or summarized
#' with methods that mirror those for \code{\link{glm}} (e.g.\ \code{\link{print.glmb}},
#' \code{\link{summary.glmb}}), as well as all the usual \code{glm}/\code{lm}
#' generics (\code{\link{predict}}, \code{\link{residuals}}, etc.).
#'
#' A helper, \code{\link{Prior_Setup}}, assists users in choosing prior
#' parameters. It ships with sensible defaults but also allows full
#' customization. In particular, the default for \code{dNormal} is a
#' reparameterization of Zellner's g-prior \insertCite{zellner1986gprior}{glmbayes}.
#'
#' Currently supported response families are
#' \code{gaussian} (identity link), \code{poisson} and \code{quasipoisson}
#' (log link), \code{gamma} (log link), and \code{binomial} and
#' \code{quasibinomial} (logit, probit, cloglog). All families support a
#' \code{dNormal} prior; the Gaussian family also offers
#' \code{dNormalGamma} and \code{dIndependent_Normal_Gamma}.
#' Two conjugate pfamilies add closed-form IID sampling for intercept-only
#' models with an identity link: \code{\link{dBeta}} for
#' \code{binomial(link = "identity")} (Beta--Binomial conjugacy) and
#' \code{\link{dGamma}(Inv_Dispersion = FALSE)} for
#' \code{poisson(link = "identity")} and \code{Gamma(link = "identity")}
#' (Gamma--Poisson and Gamma--Gamma rate conjugacy).
#' For dispersion estimation with fixed coefficients, \code{\link{dGamma}}
#' (default \code{Inv_Dispersion = TRUE}) places a Gamma prior on the
#' inverse dispersion for Gaussian and Gamma(log) models.
#'
#' For the Gaussian family, draws under \code{dNormal} and
#' \code{dNormalGamma} come from posterior distributions resulting from conjugate 
#' prior distributions \insertCite{Raiffa1961}{glmbayes}. For all other priors or response families, 
#' we use an accept-reject sampler built on the likelihood-subgradient envelope
#' method \insertCite{Nygren2006}{glmbayes}. The
#' \code{Gridtype} argument controls how many tangent points are used
#' in the envelope-trading off envelope tightness against construction
#' cost-and \code{iters} reports candidate counts before acceptance.
#'
#' By default, \code{glmb} draws \code{n = 1000} samples, uses parallel
#' CPU simulation, and-if \code{use_opencl = TRUE}-GPU-accelerated
#' envelope building. \code{"glmb"} comes with many of the same kinds of method functions
#' that come with \code{"glm"} and \code{"lm"}, so you can still call \code{\link{extractAIC}},
#' \code{\link{fitted.values}}, or any other standard method.
#'
#' The \code{\link{lmb}} function is a Bayesian version of the \code{\link[stats]{lm}} function that can 
#' be used to estimate models from the Gaussian family without the need for a \code{family} argument.
#'
#' \code{\link{rglmb}} and \code{\link{rlmb}} are functions with more minimalistic interfaces for estimating the same 
#' models without most of the internal overhead (these functions are called internally by \code{glmb} and \code{lmb}). 
#' The reduced overhead may be beneficial for Gibbs sampling implementations.
#' 
#'  
#' @author The \R implementation of \code{glmb} has been written by Kjell Nygren and
#' was built to be a Bayesian version of the \code{glm} function and hence tries
#' to mirror the features of the \code{glm} function to the greatest extent possible. For details
#' on the author(s) for the \code{glm} function see the documentation for \code{\link[stats]{glm}}.    
#' @family modelfuns
#' 
#' @seealso
#'   \code{\link[stats]{lm}}, \code{\link[stats]{glm}}, \code{\link[stats]{family}}, \code{\link[stats]{formula}} 
#'     for classical modeling functions, family objects, and formula syntax
#'
#'  \code{\link{pfamily}} for documentation of pfamily functions used to specify priors.
#' 
#'  \code{\link{Prior_Setup}}, \code{\link{Prior_Check}} for functions used to initialize and to check priors,  
#'
#'  \code{\link[glmbayesCore]{EnvelopeBuild}} for envelope construction methods.
#'
#'  Further reading: \insertCite{Nygren2006}{glmbayes};
#'  \insertCite{glmbayesChapter00,glmbayesChapterA02,glmbayesSimmethods,glmbayesChapterA08}{glmbayes};
#'  OpenCL/GPU: \insertCite{glmbayesChapter12,glmbayesChapterA10}{glmbayes}.
#'     
#'   \code{\link{summary.glmb}}, \code{\link{predict.glmb}}, \code{\link{residuals.glmb}}, \code{\link{simulate.glmb}},  
#'   \code{\link{extractAIC.glmb}}, \code{\link{dummy.coef.glmb}} and methods(class="glmb") for \code{glmb}  
#'   and the methods and generic functions for classes \code{glm} and \code{lm} from which class \code{glmb} inherits.
#'
#' @references
#' \insertAllCited{}
#' \insertRef{Dobson1990}{glmbayes}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_glmb.R
#' 
#' @export

## Portions of the model-frame/model-matrix setup, fitted-object structure,
## and method conventions follow or adapt stats::glm() and related stats
## methods. See inst/COPYRIGHTS for derived-code notices.
glmb<-function (formula, family = binomial,pfamily=dNormal(mu,Sigma,dispersion=1),n=1000,data, weights,
                use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE,
                subset,
                offset,na.action, Gridtype=2,
                n_envopt = NULL,       # envelope sizing proxy (defaults to n)
                start = NULL, etastart, 
                mustart,  control = list(...), model = TRUE, 
                method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
                ...) 
{
  call <- match.call()

  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (identical(method, "model.frame")) 
    return(mf)
  if (!is.character(method) && !is.function(method)) 
    stop("invalid 'method' argument")
  if (identical(method, "glm.fit")) 
    control <- do.call("glm.control", control)
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) 
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) 
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                    length(offset), NROW(Y)), domain = NA)
  }
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  fit <- eval(call(if (is.function(method)) "method" else method, 
                   x = X, y = Y, weights = weights, start = start, etastart = etastart, 
                   mustart = mustart, offset = offset, family = family, 
                   control = control, intercept = attr(mt, "intercept") > 
                     0L))
  if (length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <- eval(call(if (is.function(method)) "method" else method, 
                      x = X[, "(Intercept)", drop = FALSE], y = Y, weights = weights, 
                      offset = offset, family = family, control = control, 
                      intercept = TRUE))
    if (!fit2$converged) 
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    fit$null.deviance <- fit2$deviance
  }
  if (model) 
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  fit$x <- X
  fit <- c(fit, list(call = call, formula = formula, terms = mt, 
                     data = data, offset = offset, control = control, method = method, 
                     contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
                                                                             mf)))
  class(fit) <- c(fit$class, c("glm", "lm"))
  
  # Verify inputs and Initialize
  
  ## Use the prior list to set the prior elements if it is not missing
  ## Error checking to verify that the correct elements are present
  
  y<-fit$y	
  x<-fit$x
  b<-fit$coefficients
  
  prior_list=pfamily$prior_list 
  
  if (pfamily$pfamily != "dGamma") {
    mu<-as.matrix(as.vector(prior_list$mu))
    Sigma<-as.matrix(prior_list$Sigma)    

    R <- chol(Sigma)
    P <- chol2inv(R)
    P <- 0.5 * (P + t(P))
  }
  
    wtin<-fit$prior.weights	

  # NEW: normalize n_envopt (default to n, ensure integer >= 1)
  if (is.null(n_envopt)) n_envopt <- n
  if (length(n_envopt) != 1L || is.na(n_envopt) || n_envopt < 1) {
    stop("`n_envopt` must be a positive integer scalar.")
  }
  n_envopt <- as.integer(n_envopt)
  
  
  #sim<-rglmb(n=n,y=y,x=x,mu=mu,P=P,wt=wtin,dispersion=dispersion,shape=shape,rate=rate,offset2=offset,family=family,
  #           start=b,Gridtype=Gridtype)
  
  
    sim <- glmbayesCore::rglmb(n=n,y=y,x=x,family=family,pfamily=pfamily,offset=offset,
             weights=wtin,
             Gridtype=Gridtype,      n_envopt = n_envopt,     # NEW: pass through
             use_parallel = use_parallel, use_opencl = use_opencl, verbose = verbose)
  

  #dispersion2<-dispersion
  dispersion2<-sim$dispersion

  famfunc<-sim$famfunc

  
    
  # Determine family type for dispersion handling
  is_gaussian <- fit$family$family == "gaussian"
  is_quasi <- fit$family$family %in% c("quasipoisson", "quasibinomial")
  
  # Scale weights and fix dispersion if needed
  if (is_gaussian &&  length(dispersion2)==1 ) {
    
    
    wt_scaled <- wtin / dispersion2
    dispersion_fixed <- 1
    
    
  } else {
    wt_scaled <- wtin
    dispersion_fixed <- dispersion2
  }
  
  
    
  if (!is.null(prior_list$Inv_Dispersion) && isTRUE(prior_list$Inv_Dispersion)) {
    ## dGamma(Inv_Dispersion = TRUE): prior on precision/shape — no mu/Sigma
    Prior <- list(shape = prior_list$shape, rate = prior_list$rate)
  } else if (!is.null(prior_list$Inv_Dispersion) && identical(prior_list$Inv_Dispersion, FALSE)) {
    ## dGamma(Inv_Dispersion = FALSE): conjugate rate prior — has mu/Sigma surrogate
    cn <- colnames(fit$x)
    if (is.null(cn)) cn <- colnames(prior_list$Sigma)
    if (is.null(cn)) cn <- rownames(prior_list$Sigma)
    if (is.null(cn)) cn <- paste0("V", seq_along(as.vector(prior_list$mu)))
    Prior <- list(
      shape    = prior_list$shape,
      rate     = prior_list$rate,
      mean     = stats::setNames(as.vector(prior_list$mu), cn),
      Variance = prior_list$Sigma
    )
    colnames(Prior$Variance) <- cn
    rownames(Prior$Variance) <- cn
  } else {
    ## dNormal, dNormal_Gamma, dBeta, etc.: mu/Sigma only
    Prior <- list(mean = prior_list$mu, Variance = prior_list$Sigma)
    names(Prior$mean)        <- colnames(fit$x)
    colnames(Prior$Variance) <- colnames(fit$x)
    rownames(Prior$Variance) <- colnames(fit$x)
  }
  
  
  ### Set dispersion to null for quasi-families to prevent DIC from calculating


  if (!is.null(offset)) {
    if(length(dispersion2)==1){
      #    DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin/dispersion2,dispersion=dispersion2)
      
      if(fit$family$family=="quasipoisson"||fit$family$family=="quasibinomial"){
        DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin/sim$dispersion,dispersion=1)
        res=residuals(summary(sim))
        DICinfo$Deviance=rowSums(res*res)    

        DICinfo$DIC=DICinfo$DIC*sim$dispersion    
        
      }
      #else  DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
      else{DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin/sim$dispersion,dispersion=1)
      res=residuals(summary(sim))
      DICinfo$Deviance=rowSums(res*res)   

      DICinfo$DIC=DICinfo$DIC*sim$dispersion    }
      
    }
    
    if(length(dispersion2)>1){
      #  DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
      DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
    }
    
    linear.predictors<-t(offset+x%*%t(sim$coefficients))}
  

  
  if (is.null(offset)) {
    

    if(length(dispersion2)==1){
      
      if(fit$family$family=="quasipoisson"||fit$family$family=="quasibinomial"){
        DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin/sim$dispersion,dispersion=1)

        res=residuals(summary(sim))
        DICinfo$Deviance=rowSums(res*res)    

        DICinfo$DIC=DICinfo$DIC*sim$dispersion    
        #  DICinfo$DIC=NULL      
      }
      #else  DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
      else{
        

        DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin/sim$dispersion,dispersion=1)


      res=residuals(summary(sim))
      

      DICinfo$Deviance=rowSums(res*res)    

  ##    DICinfo$DIC=DICinfo$DIC*sim$dispersion
      } 
      
      
    }
    

    
    if(length(dispersion2)>1){
      #  DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
      DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
    }
    
    linear.predictors<-t(x%*%t(sim$coefficients))
    
  }
  

  
  # Only update this here so that DIC calculation above works
  dispersion2<-sim$dispersion
  
  linkinv<-fit$family$linkinv
  fitted.values<-linkinv(linear.predictors)
  
  
  outlist<-list(
    glm=fit,
    coefficients=sim$coefficients,
    coef.means=colMeans(sim$coefficients),
    coef.mode=sim$coef.mode,
    dispersion=dispersion2,
    Prior=Prior,
    fitted.values=fitted.values,
    family=fit$family,
    linear.predictors=linear.predictors,
    deviance=DICinfo$Deviance,
    pD=DICinfo$pD,
    Dbar=DICinfo$Dbar,
    Dthetabar=DICinfo$Dthetabar,
    DIC=DICinfo$DIC,
    prior.weights=fit$prior.weights,
    y=fit$y,
    x=fit$x,
    model=fit$model,
    call=fit$call,
    formula=fit$formula,
    terms=fit$terms,
    data=fit$data,
    fit=sim$fit,
    famfunc=famfunc,
    iters=sim$iters,
    contrasts=fit$contrasts,	  
    xlevels=fit$xlevels,
    pfamily=pfamily,
    simfun_call=sim$simfun_call,
    simfun_args=sim$simfun_args
    
  )
  
  outlist$call<-match.call()
  
  if (pfamily$pfamily == "dGamma" && isTRUE(pfamily$prior_list$Inv_Dispersion)) {
    class(outlist) <- c("rGamma_reg", outlist$class, "glmb", "glm", "lm")
  } else {
    class(outlist) <- c(outlist$class, "glmb", "glm", "lm")
  }
  outlist
}
#' @rdname glmb
#' @method print glmb
#' @export

print.glmb<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Posterior Mean Coefficients")
    cat(":\n")
    print.default(format(x$coef.means, digits = digits), 
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  cat("\nEffective Number of Parameters:",x$pD,"\n")
  cat("Expected Residual Deviance:",mean(x$deviance),"\n")
  cat("DIC:",x$DIC,"\n\n")
}
