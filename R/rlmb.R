#' The Bayesian Linear Model Distribution
#'
#' \code{rlmb} is used to generate iid samples from Bayesian Linear Models with multivariate normal priors. 
#' The model is specified by providing a data vector, a design matrix, and at least 2 prior constants.
#' @name 
#' rlmb
#' @aliases
#' rlmb
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x a design matrix of dimension \code{m * p}.
#' @param family a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)
#' @param pfamily the prior family to use for the model (including the constants passed to prior). 
#' @param offset an offset parameter
#' @param weights a weighting variable
#' @param digits the number of significant digits to use when printing.
#' @param \ldots additional optional arguments.
#' @return \code{rlmb} returns a object of class \code{"rglmb"}.  The function \code{summary} 
#' (i.e., \code{\link{summary.rglmb}}) can be used to obtain or print a summary of the results.
#' The generic accessor functions \code{\link{coefficients}}, \code{\link{fitted.values}},
#' \code{\link{residuals}}, and \code{\link{extractAIC}} can be used to extract
#' various useful features of the value returned by \code{\link{rglmb}}.
#' An object of class \code{"rglmb"} is a list containing at least the following components:
#' \item{coefficients}{a \code{n} by \code{length(mu)} matrix with one sample in each row}
#' \item{PostMode}{a vector of \code{length(mu)} with the estimated posterior mode coefficients}
#' \item{Prior}{A list with two components. The first being the prior mean vector and the second the prior precision matrix}
#' \item{iters}{an \code{n} by \code{1} matrix giving the number of candidates generated before acceptance for each sample.}
#' \item{famfunc}{an object of class \code{"famfunc"}}
#' \item{Envelope}{an object of class \code{"envelope"}  }
#' \item{dispersion}{the dispersion parameter used in the model}
#' \item{loglike}{a \code{n} by \code{1} matrix containing the negative loglikelihood for each sample.}
#' @details The \code{rglmb} function produces iid samples for Bayesian generalized linear 
#' models. It has a more minimialistic interface than than the \code{\link{glmb}} 
#' function. Core required inputs for the function include the data vector, the design  
#' matrix and a prior specification. In addition, the dispersion parameter must 
#' currently be provided for the gaussian, Gamma, quasipoisson, and quasibinomial 
#' families (future implementations may incorporate a prior for these into the 
#' \code{rglmb} function).
#' 
#' Current implemented models include the gaussian family (identity link function), the
#' poisson and quasipoisson families (log link function), the gamma family (log link 
#' function), as well as the binomial and quasibinomial families (logit, probit, and 
#' cloglog link functions). The function returns the simulated Bayesian coefficients 
#' and some associated outputs.
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
#' 
#' Raiffa, Howard and Schlaifer, R (1961)
#' \emph{Applied Statistical Decision Theory.}
#' Boston: Clinton Press, Inc.
#' 
#' Venables, W. N. and Ripley, B. D. (2002)
#' \emph{Modern Applied Statistics with S.}
#' New York: Springer.
#' 
#' @example inst/examples/Ex_rglmb.R
#' @export
#' @rdname rlmb
#' @order 1


rlmb<-function(n=1,y,x,family=gaussian(),pfamily,
               offset=rep(0,nobs),weights=1)
  {

  ## Pull in information from the pfamily  
  pf=pfamily$pfamily
  #okfamilies=pfamily$okfamilies  
  okfamilies <- c("gaussian")    # Only gaussian okfamily for rlmb
  plinks=pfamily$plinks
  prior_list=pfamily$prior_list 
  simfun=pfamily$simfun
  
    
#  if(is.numeric(n)==FALSE||is.numeric(y)==FALSE||is.numeric(x)==FALSE||
#     is.numeric(mu)==FALSE||is.numeric(P)==FALSE) stop("non-numeric argument to numeric function")
  
  x <- as.matrix(x)
  mu<-as.matrix(as.vector(mu))
  P<-as.matrix(P)    
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) 
    rownames(y)
  else names(y)
  if(length(n)>1) n<-length(n)	   
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(offset2)) 
    offset2 <- rep(0, nobs)
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
  
  if (is.null(start)) 
    start <- mu
  if (is.null(offset2)) 
    offset2 <- rep.int(0, nobs)
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  #  if(family$family %in% okfamilies){
  #    if(family$family=="gaussian") oklinks<-c("identity")
  
  if(family$family %in% okfamilies){
    oklinks<-c("identity")
    if(!family$link %in% oklinks){      
      stop(gettextf("link \"%s\" not available for selected family; available links are %s", 
                    family$link , paste(sQuote(oklinks), collapse = ", ")), 
           domain = NA)
    }
  }
  else{
    stop(gettextf("family \"%s\" not available for current pfamily; available families are %s", 
                  family$family , paste(sQuote(okfamilies), collapse = ", ")), 
         domain = NA)
    
  }
  
  
  ## Call relevant simulation function (for now without control2 list)
  
  outlist=simfun(n=n,y=y,x=x,prior_list=prior_list,offset=offset,weights=weights,family=family)
  
  return(outlist)
  

#    dispersion2=dispersion
#    if(is.null(dispersion)){dispersion2=0}
#    if(dispersion2>0){outlist<-.rnorm_reg_cpp(n=n,y=y,x=x,mu=mu,P=P,offset2=offset2,wt=wt,dispersion=dispersion,famfunc=famfunc,f1=f1,f2=f2,f3=f3,start=mu)}
#    else{
#      ## For now, build list here - will be created by pfamily
#      prior_list2=list(mu=mu,Sigma=solve(P),P=P,shape=shape,rate=rate)        
#      outlist=rnorm_gamma_reg(n=n,y=y,x=x,prior_list=prior_list2,offset=offset2,weights=wt)  
#      }
    

#  colnames(outlist$coefficients)<-colnames(x)
  
#  outlist$call<-match.call()
  
#  class(outlist)<-c(outlist$class,"rlmb")
#  outlist
  
}

#' @rdname rlmb
#' @order 2
#' @export
#' @method print rlmb


print.rlmb<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
 		
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
 if (length(coef(x))) {
        cat("Simulated Coefficients")
        cat(":\n")
        print.default(format(x$coefficients, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
 }


