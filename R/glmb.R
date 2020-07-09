#' Fitting Bayesian Generalized Linear Models
#'
#' \code{glmb} is used to fit Bayesian generalized linear models, specified by giving a symbolic descriptions of 
#' the linear predictor,  the error distribution, and the prior distribution.
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
#' @param digits the number of significant digits to use when printing.
#' @inheritParams stats::glm
#' @return \code{glmb} returns an object of class \code{"glmb"}. The function \code{summary} (i.e., 
#' \code{\link{summary.glmb}}) can be used to obtain or print a summary of the results.  The generic accessor functions 
#' \code{\link{coefficients}}, \code{\link{fitted.values}}, \code{\link{residuals}}, and \code{\link{extractAIC}} can be used 
#' to extract various useful features of the value returned by code{\link{glmb}}.
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
#' @details The function \code{glmb} is a Bayesian version of the classical \code{\link{glm}} function.  Setup of
#' the models mirrors those for the \code{\link{glm}} function but add additional required arguments
#' \code{mu} and \code{Sigma}, representing a multivariate normal prior. In addition, the dispersion 
#' parameter must currently be provided for the gaussian, Gamma, quasipoisson, and quasibinomial families 
#' (future implementations may incoporate priors for these into the \code{glmb} function).  The function 
#' generates random iid samples from the posterior density associated with the model (assuming a fixed 
#' dispersion parameter). 
#' 
#' Current implemented models include the gaussian family (identity link function), the poisson and
#' quasipoisson families (log link function), the gamma family (log link function), as well as the 
#' binomial and quasibinomial families (logit, probit, and cloglog link functions).
#' 
#' The function returns the output from a call to the function \code{\link{glm}} as well as the simulated 
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
#' @author The \R implementation of \code{glmb} has been written by Kjell Nygren and
#' was built to be a Bayesian version of the \code{glm} function and hence tries
#' to mirror the features of the \code{glm} function to the greatest extent possible. For details
#' on the author(s) for the \code{glm} function see the documentation for \code{\link[stats]{glm}}.    
#' @family modelfuns
#' @seealso \code{\link[stats]{lm}} and \code{\link[stats]{glm}} for classical modeling functions.
#' 
#' \code{\link{pfamily}} for documentation of pfamily functions used to specify priors.
#' 
#' \code{\link{Prior_Setup}}, \code{\link{Prior_Check}} for functions used to initialize and to check priors,  
#'
#' \code{\link{summary.glmb}}, \code{\link{predict.glmb}}, \code{\link{residuals.glmb}}, \code{\link{simulate.glmb}}, 
#' \code{\link{extractAIC.glmb}}, \code{\link{dummy.coef.glmb}} and methods(class="glmb") for \code{glmb} 
#' and the methods and generic functions for classes \code{glm} and \code{lm} from which class \code{glmb} inherits.
#' 
#'     
#' @references 
#' Dobson, A. J. (1990)
#' \emph{An Introduction to Generalized Linear Models.}
#' London: Chapman and Hall.
#' 
#' Hastie, T. J. and Pregibon, D. (1992)
#' \emph{Generalized linear models.}
#' Chapter 6 of \emph{Statistical Models in S}
#' eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#' McCullagh P. and Nelder, J. A. (1989)
#' \emph{Generalized Linear Models.}
#' London: Chapman and Hall.
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
#' Venables, W. N. and Ripley, B. D. (2002)
#' \emph{Modern Applied Statistics with S.}
#' New York: Springer.
#' 
#' @example inst/examples/Ex_glmb.R
#' 
#' @export glmb
#' @exportClass glmb

glmb<-function (formula, family = binomial,pfamily=dNormal(mu,Sigma,dispersion=1),n=1000,data, weights, subset,
                offset,na.action, Gridtype=1,start = NULL, etastart, 
                mustart,  control = list(...), model = TRUE, 
                method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
                ...) 
{
  # Note, renamed subset argument to subset2 as it caused conflict with subset function
  
  call <- match.call()
  
  # In earlier iteration, subset did not work so used subset2 as argument.
  # Trying to rename argument and then setting subset2=subset as part of process of changing
  
  #subset2=subset # line 153, 
  
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
  
  mu<-as.matrix(as.vector(prior_list$mu))
  Sigma<-as.matrix(prior_list$Sigma)    
  P<-solve(prior_list$Sigma) 
  wtin<-fit$prior.weights	
  
  
  #sim<-rglmb(n=n,y=y,x=x,mu=mu,P=P,wt=wtin,dispersion=dispersion,shape=shape,rate=rate,offset2=offset,family=family,
  #           start=b,Gridtype=Gridtype)
  
  
  
  sim<-rglmb(n=n,y=y,x=x,family=family,pfamily=pfamily,offset=offset,weights=wtin)
  
  
  
  dispersion2<-sim$dispersion
  famfunc<-sim$famfunc
  
  Prior<-list(mean=prior_list$mu,Variance=prior_list$Sigma)
  names(Prior$mean)<-colnames(fit$x)
  colnames(Prior$Variance)<-colnames(fit$x)
  rownames(Prior$Variance)<-colnames(fit$x)
  
  
  if (!is.null(offset)) {
    if(length(dispersion2)==1){
        #    DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin/dispersion2,dispersion=dispersion2)
            DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
            
            }
    
    if(length(dispersion2)>1){
    #  DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
      DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
    }
    
    linear.predictors<-t(offset+x%*%t(sim$coefficients))}
  if (is.null(offset)) {
    if(length(dispersion2)==1){
      #DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin/dispersion2,dispersion=dispersion2)
      DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
      
          }
    
    if(length(dispersion2)>1){
    #  DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
      DICinfo<-DIC_Info(sim$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
    }
    
    linear.predictors<-t(x%*%t(sim$coefficients))
    
  }
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
    famfunc=famfunc,
    iters=sim$iters,
    contrasts=fit$contrasts,	  
    xlevels=fit$xlevels,
    pfamily=pfamily
  )
  
  outlist$call<-match.call()
  
  class(outlist)<-c(outlist$class,"glmb","glm","lm")
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



#' Calculates DIC and Deviance Information
#'
#' Caculates DIC and Deviance Information for a fitted model.
#' @param coefficients A matrix with coefficients from the rglmb function
#' @param y a vector of observations of length m
#' @param x a design matrix of dimension m*p
#' @param alpha an offset parameter or vector
#' @param f1 A function used to calculate the negative log-likelihood
#' @param f4 A function used to calculate the deviance
#' @param wt a vector of weights
#' @param dispersion dispersion parameter
#' 
#' @details Calculates DIC and Deviance Information
#' @return A list with the following components
#' \item{Deviance}{A \code{n * 1} matrix with the deviance for each draw}
#' \item{Dbar}{Mean for negative 2 times negative log-likelihood}
#' \item{Dthetabar}{Negative 2 times log-likelihood evaluated at mean parameters}
#' \item{pD}{Effective number of parameters}
#' \item{DIC}{DIC statistic}
#' @example inst/examples/Ex_glmbdic.R
#' @noRd 

DIC_Info<-function(coefficients,y,x,alpha=0,f1,f4,wt=1,dispersion=1){
  
  l1<-length(coefficients[1,])
  l2<-length(coefficients[,1])
  
  D<-matrix(0,nrow=l2,ncol=1)
  D2<-matrix(0,nrow=l2,ncol=1)
  
  #disp_temp=dispersion
  #if(length(dispersion)==1) disp_temp=rep(dispersion,length(y))
  
  
  
  if(length(dispersion)==1){
    for(i in 1:l2){
      b<-as.vector(coefficients[i,])
      # No offset! [Edit to multiply by dispersion to get residual deviance]
      
      # This calculation likely has cancellation effects
      #  Passing dispersion as both     
      #D[i,1]<-dispersion*f4(b=b,y=y,x=x,alpha=alpha,wt=wt/dispersion,dispersion=dispersion)
      D[i,1]<-f4(b=b,y=y,x=x,alpha=alpha,wt=wt,dispersion=dispersion)
      
      D2[i,1]<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt)
    }
    
    Dbar<-mean(D2)
    
    b<-colMeans(coefficients)
    
    Dthetabar<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt)
    
  }
  
  if(length(dispersion)>1){
    for(i in 1:l2){
      b<-as.vector(coefficients[i,])
      #D[i,1]<-f4(b=b,y=y,x=x,alpha=alpha,wt=wt/dispersion[i],dispersion=dispersion[i])
      D[i,1]<-f4(b=b,y=y,x=x,alpha=alpha,wt=wt,dispersion=dispersion[i])
      
      D2[i,1]<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt/dispersion[i])
    }
    
    Dbar<-mean(D2)
    
    b<-colMeans(coefficients)
    # Use the mean for the dispersion in Dthetabar calculation
    # using mean for 1/dispersion (precision) would change the DIC
    dispbar<-mean(dispersion)
    Dthetabar<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt/dispbar)
    
    
  }
  
  
  pD<-Dbar-Dthetabar
  DIC<-pD+Dbar
  list(Deviance=D,Dbar=Dbar,Dthetabar=Dthetabar,pD=pD,DIC=DIC)
}





