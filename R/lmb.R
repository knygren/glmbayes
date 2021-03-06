#' Fitting Bayesian Linear Models
#'
#' \code{lmb} is used to fit Bayesian linear models, specified by giving a symbolic descriptions of the linear 
#' predictor and the prior distribution.
#' @aliases
#' lmb
#' print.lmb
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param pfamily a description of the prior distribution and associated constants to be used in the model. This
#' should be a pfamily function (see \code{\link{pfamily}} for details of pfamily functions).
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
#' @return \code{glmb} returns an object of class \code{"glmb"}. The function \code{summary} (i.e., 
#' \code{\link{summary.glmb}}) can be used to obtain or print a summary of the results.  The generic accessor functions 
#' \code{\link{coefficients}}, \code{\link{fitted.values}}, \code{\link{residuals}}, and \code{\link{extractAIC}} can be used 
#' to extract various useful features of the value returned by code{\link{glmb}}.
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
#' Objects of class \code{"glmb"} are normall of class \code{c("glmb","glm","lm")},
#' that is inherit from classes \code{glm} and \code{lm} and well-designed
#' methods from those classed will be applied when appropriate.
#' 
#' If a \link{binomial} \code{glmb} model was specified by giving a two-column 
#' response, the weights returned by \code{prior.weights} are the total number of
#' cases (factored by the supplied case weights) and the component of \code{y}
#' of the result is the proportion of successes.
#' 
#' @details The function \code{lmb} is a Bayesian version of the classical \code{\link{lm}} function.  Setup of
#' the models mirrors those for the \code{\link{lm}} function but add the required argument \code{\link{pfamily}}
#' with a prior formulation. The function generates random iid samples from the posterior density associated 
#' with the model. 
#'  
#' The function returns the output from a call to the function \code{\link{lm}} as well as the simulated 
#' Bayesian coefficients and associated outputs. In addition, the function returns model diagnostic 
#' information related to the \code{DIC}, a Bayesian Information Criterion similar to the \code{AIC} for 
#' classical models. 
#' 
#' For the gaussian family, iid samples from the posterior density are genererated using 
#' standard simulation procedures for multivariate normal densities. For all other 
#' families, the samples are generated using accept-reject procedures by leveraging the 
#' likelihood subgradient approach of Nygren and Nygren (2006). This approach relies on
#' tight enveloping functions that bound the posterior density from above. The Gridtype 
#' parameter is used to select the method used for determining the size of a Grid used 
#' to build the enveloping function. See \code{\link{EnvelopeBuild}} for details. 
#' Depending on the selection, the time to build the envelope and the acceptance rate 
#' during the simulation process may vary. The returned value \code{iters} contains the 
#' number of candidates generated before acceptance for each draw.
#'
#' @author The \R implementation of \code{lmb} has been written by Kjell Nygren and
#' was built to be a Bayesian version of the \code{lm} function and hence tries
#' to mirror the features of the \code{lm} function to the greatest extent possible while also taking advantage
#' of some of the method functions developed for the \code{glmb} function. For details
#' on the author(s) for the \code{lm} function see the documentation for \code{\link[stats]{lm}}.
#'     
#' @references 
#' Chambers, J.M.(1992) \emph{Linear models.} Chapter 4 of \emph{Statistical Models in S}
#' eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#' 
#' Wilkinson, G.N. and Rogers, C.E. (1973). Symbolic descriptions of factorial models for 
#' analysis of variance. \emph{Applied Statistics}, \bold{22}, 392-399.
#' doi: \href{https://doi.org/10.2307/2346786}{10.2307/2346786}.
#' 
#' Nygren, K.N. and Nygren, L.M (2006)
#' Likelihood Subgradient Densities. \emph{Journal of the American Statistical Association}.
#' vol.101, no.475, pp 1144-1156.
#' doi: \href{https://doi.org/10.1198/016214506000000357}{10.1198/016214506000000357}.
#' 
#' Raiffa, Howard and Schlaifer, R (1961)
#' \emph{Applied Statistical Decision Theory.}
#' Boston: Clinton Press, Inc.
#'
#' @family modelfuns
#' @seealso The classical modeling functions \code{\link[stats]{lm}} and \code{\link[stats]{glm}}.
#' 
#' \code{\link{pfamily}} for documentation of pfamily functions used to specify priors.
#' 
#' \code{\link{Prior_Setup}}, \code{\link{Prior_Check}} for functions used to initialize and to check priors,  
#'
#' \code{\link{summary.glmb}}, \code{\link{predict.glmb}}, \code{\link{simulate.glmb}}, 
#' \code{\link{extractAIC.glmb}}, \code{\link{dummy.coef.glmb}} and methods(class="glmb") for methods 
#' inherited from class \code{glmb} and the methods and generic functions for classes \code{glm} and 
#' \code{lm} from which class \code{lmb} also inherits.
#'
#' @example inst/examples/Ex_lmb.R
#' 
#' @export lmb
#' @exportClass lmb


lmb <- function ( formula, pfamily, n=1000,data, subset, weights, na.action,method = "qr", model = TRUE, x = TRUE, 
                  y = TRUE,qr = TRUE, singular.ok = TRUE, contrasts = NULL,offset, ...)
{
  
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
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
  if(!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }
  
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (is.matrix(y))
      matrix(,0,3) else numeric(), residuals = y,
      fitted.values = 0 * y, weights = w, rank = 0L,
      df.residual = if(!is.null(w)) sum(w != 0) else
        if (is.matrix(y)) nrow(y) else length(y))
    if(!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    z <- if(is.null(w)) lm.fit(x, y, offset = offset,
                               singular.ok=singular.ok, ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok=singular.ok, ...)
  }
  class(z) <- c(if(is.matrix(y)) "mlm", "lm")
  
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
  
  ######   End of lm function
  # Verify inputs and Initialize
  
  
  #    if(is.numeric(n)==FALSE||is.numeric(mu)==FALSE||is.numeric(Sigma)==FALSE) stop("non-numeric argument to numeric function")
  ## Pull in information from the pfamily  
  
  prior_list=pfamily$prior_list 
  y<-z$y	
  x<-z$x
  b<-z$coefficients
  mu<-as.matrix(as.vector(prior_list$mu))
  Sigma<-as.matrix(prior_list$Sigma)    
  dispersion=prior_list$dispersion
  P<-solve(prior_list$Sigma) 
  
  if(is.null(z$weights))     wtin<-rep(1,length(y))
  else wtin=z$weights	
  
  sim<-rlmb(n=n,y=y,x=x,pfamily=pfamily,weights=wtin,
            offset=offset)
  
  
  
  dispersion2<-sim$dispersion
  
  famfunc<-sim$famfunc
  
  
  Prior<-list(mean=as.numeric(mu),Variance=Sigma)
  names(Prior$mean)<-colnames(z$x)
  colnames(Prior$Variance)<-colnames(z$x)
  rownames(Prior$Variance)<-colnames(z$x)
  
  
  if (!is.null(offset)) {
    
    if(length(dispersion2)==1){
      #DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin/dispersion2,dispersion=dispersion2)
      DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
      
          }
    
    if(length(dispersion2)>1){
      #DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
      DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
      
          }
    
    linear.predictors<-t(offset+x%*%t(sim$coefficients))
    fitted.values=linear.predictors
    
  }
  
  
  if (is.null(offset)) {
    
    if(length(dispersion2)==1){
      #DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin/dispersion2,dispersion=dispersion2)
      DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
      
          }
    
    if(length(dispersion2)>1){
      DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)

          }
    
    
    linear.predictors<-t(x%*%t(sim$coefficients))
    fitted.values=linear.predictors
    
  }
  
  residuals=fitted.values
  
  for(i in 1:n){
    
    residuals[i,1:length(y)]=y-residuals[i,1:length(y)]
  }	
  
  #linkinv<-z$family$linkinv
  outlist<-list(
    lm=z,
    coefficients=sim$coefficients,
    coef.means=colMeans(sim$coefficients),
    coef.mode=sim$coef.mode,
    dispersion=dispersion2,
    residuals=residuals,
    Prior=Prior,
    fitted.values=fitted.values,
    #family=fit$family,
    linear.predictors=linear.predictors,
    deviance=DICinfo$Deviance,
    pD=DICinfo$pD,
    Dbar=DICinfo$Dbar,
    Dthetabar=DICinfo$Dthetabar,
    DIC=DICinfo$DIC,
    prior.weights=wtin,
    weights=wtin,
    offset=offset,
    y=z$y,
    x=z$x,
    model=z$model,
    call=z$call,
    formula=z$formula,
    terms=z$terms,
    data=z$data,
    fit=sim$fit,
    famfunc=famfunc,
    iters=sim$iters,
    contrasts=z$contrasts,	  
    xlevels=z$xlevels,
    pfamily=pfamily
    
  )
  
  outlist$call<-match.call()
  
  class(outlist)<-c(outlist$class,"lmb","glmb","glm","lm")
  outlist
}

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

