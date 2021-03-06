#' The Bayesian Linear Model Distribution
#'
#' \code{rlmb} is used to generate iid samples from Bayesian Linear Models with multivariate normal priors. 
#' The model is specified by providing a data vector, a design matrix, and a pfamily (determining the 
#' prior distribution).
#' @name 
#' rlmb
#' @aliases
#' rlmb
#' rlmb.print
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x for \code{rlmb} a design matrix of dimension \code{m * p} and for 
#' \code{print.rlmb} the object to be printed. 
#' @param pfamily a description of the prior distribution and associated constants to be used in the model. This
#' should be a pfamily function (see \code{\link{pfamily}} for details of pfamily functions.)
#' @param digits the number of significant digits to use when printing.
#' @inheritParams lmb
#' @return \code{rlmb} returns a object of class \code{"rlmb"}.  The function \code{summary} 
#' (i.e., \code{\link{summary.rglmb}}) can be used to obtain or print a summary of the results.
#' The generic accessor functions \code{\link{coefficients}}, \code{\link{fitted.values}},
#' \code{\link{residuals}}, and \code{\link{extractAIC}} can be used to extract
#' various useful features of the value returned by \code{\link{rglmb}}.
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
#' @details The \code{rlmb} function produces iid samples for Bayesian generalized linear 
#' models. It has a more minimialistic interface than than the \code{\link{lmb}} 
#' function. Core required inputs for the function include the data vector, the design  
#' matrix and a prior specification. In addition, the dispersion parameter must 
#' currently be provided for the gaussian, Gamma, quasipoisson, and quasibinomial 
#' families (future implementations may incorporate a prior for these into the 
#' \code{rlmb} function).
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
#' @example inst/examples/Ex_rlmb.R
#' @export
#' @rdname rlmb
#' @order 1


rlmb<-function(n=1,y,x,pfamily,offset=rep(0,nobs),weights=NULL)
  {


  ## Pull in information from the pfamily  
  pf=pfamily$pfamily
  #okfamilies=pfamily$okfamilies  
  okfamilies <- c("gaussian")    # Only gaussian is okfamily for rlmb (different from rglmb)
  plinks=pfamily$plinks
  prior_list=pfamily$prior_list 
  simfun=pfamily$simfun

  family=gaussian()
  
#  if(is.numeric(n)==FALSE||is.numeric(y)==FALSE||is.numeric(x)==FALSE||
#     is.numeric(mu)==FALSE||is.numeric(P)==FALSE) stop("non-numeric argument to numeric function")

  #P=solve(prior_list$Sigma)
    
  x <- as.matrix(x)
  #mu<-as.matrix(as.vector(prior_list$mu))
  #P<-as.matrix(P)    
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) 
    rownames(y)
  else names(y)
  if(length(n)>1) n<-length(n)	   
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(offset)) 
    offset <- rep(0, nobs)
  
  #nvars2<-length(mu)	
  #if(!nvars==nvars2) stop("incompatible dimensions")
  #if (!all(dim(P) == c(nvars2, nvars2))) 
  #  stop("incompatible dimensions")
  
  
  #if(!isSymmetric(P))stop("matrix P must be symmetric")

    
  if(is.null(weights)) weights=rep(1,nobs)
  if(length(weights)==1) weights=rep(weights,nobs)
  nobs2=length(weights)
  nobs3=NROW(x)
  nobs4=NROW(offset)
  
  if(nobs2!=nobs) stop("weighting vector must have same number of elements as y")
  if(nobs3!=nobs) stop("matrix X must have same number of rows as y")
  if(nobs4!=nobs) stop("offset vector must have same number of rows as y")
  
  #tol<- 1e-06 # Link this to Magnitude of P	
  #eS <- eigen(P, symmetric = TRUE,only.values = FALSE)
  #ev <- eS$values
  #if (!all(ev >= -tol * abs(ev[1L]))) 
  #  stop("'P' is not positive definite")

  #if (is.null(start)) 
  #  start <- mu
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
  


  outlist=simfun(n=n,y=y,x=x,prior_list=prior_list,offset=offset,weights=weights,family=family)


  outlist$pfamily=pfamily

  return(outlist)
  


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


