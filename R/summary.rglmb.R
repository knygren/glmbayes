#' Summarizing Bayesian Generalized Linear Model Distribution Functions
#'
#' These functions are all \code{\link{methods}} for class \code{rglmb},
#' \code{rlmb}, or \code{summary.rglmb} objects.
#' 
#' @aliases 
#' summary.rglmb
#' summary.rlmb
#' print.summary.rglmb
#' @param object an object of class \code{"rglmb"} or \code{"rlmb"} for which a 
#' summary is desired.
#' @param x an object of class \code{"summary.rglmb"} for which a printed output is desired.
#' @param digits the number of significant digits to use when printing.
#' @param \ldots Additional optional arguments
#' @return \code{summary.rglmb} returns a object of class \code{"summary.rglmb"}, a 
#' list with components: 
#' \item{n}{number of draws generated}
#' \item{coefficients1}{Matrix with the prior mean and approximate weight for the prior relative to the data}
#' \item{coefficients}{Matrix with columns for the posterior mode, posterior mean, posterior standard
#' deviation, monte carlo error, and tail probabilities (posterior probability of observing a 
#' value for the coefficient as extreme as the prior mean)}
#' \item{Percentiles}{Matrix with estimated percentiles associated with the posterior density}
#' @details The \code{summary.rglmb} function summarizes the output from the 
#' \code{\link{rglmb}} or \code{\link{rlmb}} function. It takes an object of class \code{rglmb}
#' or \code{rlmb} as an input. For \code{dGamma} rate-prior fits, it delegates to
#' \code{summary} on the dispersion draws directly.
#' The output to a large extent mirrors the output from the \code{\link{summary.glm}} 
#' function. This is particularly true for the print output from the function 
#' (i.e. the output of the function \code{\link{print.summary.rglmb}}).
#' 
#' @seealso \code{\link{rglmb}}; \code{\link{glmb}}, \code{\link{glmbayes-package}},
#'   \code{\link{rlmb}}, \code{\link{lmb}}, \code{\link{summary}},
#'   \code{\link[stats]{summary.lm}}, \code{\link[stats]{summary.glm}}.
#' @example inst/examples/Ex_summary.rglmb.R
#' @export
#' @method summary rglmb


summary.rglmb<-function(object,...){


  ## Detect prior type
  prior_type <- attr(object$pfamily, "Prior Type")
  
  
  if (!is.null(prior_type) && prior_type=="dGamma") {
    return(summary(object$dispersion))
  }
  
  ### Pull in information needed to compute linear.predictors, 
  ## fitted.values, and DICInfo  

  offset=.rglmb_get_offset(object)
  dispersion2=object$dispersion
  famfunc=object$famfunc
  y=object$y
  x=object$x
  ## Need to check if prior weights adjust for dispersion  
  ## This should calculate but need to verify DIC is correct
  
  wtin=object$prior.weights  


  ##############################################


    
  if (!is.null(offset)) {
    if(length(dispersion2)==1){
      
      DICinfo<-DIC_Info(object$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
    }
    
    if(length(dispersion2)>1){
      DICinfo<-DIC_Info(object$coefficients,y=y,x=x,alpha=offset,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
    }
    
    linear.predictors<-t(offset+x%*%t(object$coefficients))}
  if (is.null(offset)) {
    if(length(dispersion2)==1){
      
      DICinfo<-DIC_Info(object$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
    }
    
    if(length(dispersion2)>1){
      DICinfo<-DIC_Info(object$coefficients,y=y,x=x,alpha=0,f1=famfunc$f1,f4=famfunc$f4,wt=wtin,dispersion=dispersion2)
    }
    
    linear.predictors<-t(x%*%t(object$coefficients))
    
  }
    
  linkinv<-object$family$linkinv
  fitted.values<-linkinv(linear.predictors)


  ##################################################
  
  n<-length(object$coefficients[,1])  
  l1<-length(object$coef.mode)
  percentiles<-matrix(0,nrow=l1,ncol=7)
  se<-sqrt(diag(var(object$coefficients)))
  mc<-se/sqrt(n)
  

  Prec <- .rglmb_prior_precision(object$Prior)
  R <- chol(Prec)
  Prec_inv <- chol2inv(R)
  Prec_inv <- 0.5 * (Prec_inv + t(Prec_inv))   # enforce symmetry

    priorrank<-matrix(0,nrow=l1,ncol=1)
  pval1<-matrix(0,nrow=l1,ncol=1)
  pval2<-matrix(0,nrow=l1,ncol=1)
  

  for(i in 1:l1){
    percentiles[i,]<-quantile(object$coefficients[,i],probs=c(0.01,0.025,0.05,0.5,0.95,0.975,0.99))
    test<-append(object$coefficients[,i],object$Prior$mean[i])
    test2<-rank(test)
    priorrank[i,1]<-test2[n+1]
    pval1[i,1]<-priorrank[i,1]/(n+1)
    pval2[i,1]<-min(pval1[i,1],1-pval1[i,1])
    
  }

  ## Add call to glm to recover the mle information
  glm_mle=glm(y~x-1,family=object$family,weights=wtin)
  ml<-coef(glm_mle) 
  se1<-sqrt(diag(vcov(glm_mle)))
    
  Tab1 <- cbind(
    "Prior Mean" = as.numeric(object$Prior$mean),
    "Prior.sd"   = as.numeric(sqrt(diag(Prec_inv))),
    "Max Like."  = as.numeric(ml),
    "Like.sd"    = as.numeric(se1)
  )
      
  TAB<-cbind("Post.Mode"=as.numeric(object$coef.mode),
             "Post.Mean"=as.numeric(colMeans(coef(object))),
             "Post.Sd"=se,"MC Error"=as.numeric(mc),
             "Pr(tail)"=as.numeric(pval2))
  TAB2<-cbind("1.0%"=percentiles[,1],
              "2.5%"=percentiles[,2],
              "5.0%"=percentiles[,3],
              Median=as.numeric(percentiles[,4]),
              "95.0%"=percentiles[,5],
              "97.5%"=as.numeric(percentiles[,6]),"99.0%"=
                as.numeric(percentiles[,7])
              )

  # Recover coefficient names from the most reliable source
  coef_names <- colnames(object$coefficients)
  
  if (is.null(coef_names)) {
    coef_names <- names(object$coef.mode)
  }
  
  if (is.null(coef_names)) {
    coef_names <- colnames(object$x)
  }
  
  if (is.null(coef_names)) {
    coef_names <- paste0("V", seq_len(ncol(object$coefficients)))
  }
  
  # Apply names to all summary tables
  rownames(Tab1) <- coef_names
  rownames(TAB)  <- coef_names
  rownames(TAB2) <- coef_names

  res<-list(
    coefficients=object$coefficients,
    coef.means=colMeans(object$coefficients),
    coef.mode=object$coef.mode,
    dispersion=mean(object$dispersion),
    Prior=object$Prior,
    fitted.values=fitted.values,
    family=family(glm_mle),
    linear.predictors=linear.predictors,
    deviance=DICinfo$Deviance,
    pD=DICinfo$pD,
    Dbar=DICinfo$Dbar,
    Dthetabar=DICinfo$Dthetabar,
    DIC=DICinfo$DIC,
    prior.weights=object$prior.weights,
    y=object$y,
    x=object$x,
    model=model.frame(glm_mle),
    call = object$call,
    formula=object$formula,
    data=object$data,
    famfunc=object$famfunc,
    iters=object$iters,
    Envelope=object$Envelope,
    loglike=object$loglike,
    n=n,
    coefficients.Tab0=Tab1,
    coefficients.Tab1=TAB,
    Percentiles=TAB2)    
  
  class(res)<-c(res$class,"summary.rglmb","glmb","glm","lm")
  
  res
  
}


#' @rdname summary.rglmb
#' @export
#' @method summary rlmb

summary.rlmb <- function(object, ...) {
  summary.rglmb(object, ...)
}


#' @rdname summary.rglmb
#' @export
#' @method print summary.rglmb

print.summary.rglmb<-function(x,digits = max(3, getOption("digits") - 3),...){
  cat("Call\n")
  print(x$call)
    cat("\nExpected Deviance Residuals:\n")
    mres<-colMeans(residuals(x))
  if (length(mres) > 5) {
    fn   <- fivenum(mres)
    names(fn) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(fn)
  } else {
    print(mres)
  }
  
  cat("\nPrior Estimates with Standard Deviations\n\n")
  printCoefmat(x$coefficients.Tab0,digits=digits)
  cat("\nBayesian Estimates Based on",x$n,"iid draws\n\n")
  printCoefmat(x$coefficients.Tab1,digits=digits,P.values=TRUE,
               has.Pvalue=TRUE)
  cat("\nDistribution Percentiles\n\n")
  printCoefmat(x$Percentiles,digits=digits)
  cat("\nEffective Number of Parameters:",x$pD,"\n")
  cat("Expected Residual Deviance:",mean(x$deviance),"\n")
  cat("DIC:",x$DIC,"\n\n")
  cat("Expected Mean dispersion:",x$dispersion,"\n")
  cat("Sq.root of Expected Mean dispersion:",sqrt(x$dispersion),"\n\n")
  cat("Mean Likelihood Subgradient Candidates Per iid sample:",mean(x$iters),"\n\n")

}


## Offset recovery for summary.rglmb(): most rglmb()/rlmb() outputs carry
## offset2 (a copy of the input offset, or a zero vector -- see
## R/simfunction.R); simfun_args$offset is a defensive fallback for any other
## simfun-produced object that stores it there instead. Matches
## glmbayesCore's `.rglmb_get_offset()`.
.rglmb_get_offset <- function(object) {
  if (!is.null(object$offset2)) {
    return(object$offset2)
  }
  if (!is.null(object$simfun_args$offset)) {
    return(object$simfun_args$offset)
  }
  rep(0, NROW(object$y))
}

## Prior-precision recovery for summary.rglmb(): most priors carry
## object$Prior$Precision directly; this falls back to inverting Sigma/
## Variance for priors that only store the covariance. Matches
## glmbayesCore's `.rglmb_prior_precision()`.
.rglmb_prior_precision <- function(prior) {
  if (!is.null(prior$Precision)) {
    P <- prior$Precision
    return(0.5 * (P + t(P)))
  }
  if (!is.null(prior$Sigma)) {
    V <- prior$Sigma
  } else if (!is.null(prior$Variance)) {
    V <- prior$Variance
  } else {
    stop("Could not recover prior precision from object$Prior.", call. = FALSE)
  }
  P <- chol2inv(chol(V))
  0.5 * (P + t(P))
}
