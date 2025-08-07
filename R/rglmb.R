#' The Bayesian Generalized Linear Model Distribution
#'
#' \code{rglmb} is used to generate iid samples for Bayesian Generalized Linear Models.
#' The model is specified by providing a data vector, a design matrix, 
#' the family (determining the likelihood function) and the pfamily (determining the 
#' prior distribution).
#' @param y a vector of observations of length \code{m}.
#' @param x for \code{rglmb} a design matrix of dimension \code{m * p} and for \code{print.rglmb} the object to be printed. 
#' @inheritParams glmb
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
#' @author The \R implementation of \code{rglmb} has been written by Kjell Nygren and
#' was built to be a Bayesian version of the \code{glm} function but with a more minimalistic interface 
#' than the \code{glmb} function. It also borrows some of its structure from other random generating function 
#' like \code{\link{rnorm}} and hence the \code{r} prefix. 
#' 
#' @family modelfuns
#' @seealso \code{\link[stats]{lm}} and \code{\link[stats]{glm}} for classical modeling functions.
#' 
#' \code{\link{family}} for documentation of family functions used to specify priors.

#' \code{\link{pfamily}} for documentation of pfamily functions used to specify priors.
#' 
#' \code{\link{Prior_Setup}}, \code{\link{Prior_Check}} for functions used to initialize and to check priors,  
#'
#' \code{\link{summary.glmb}}, \code{\link{predict.glmb}}, \code{\link{residuals.glmb}}, \code{\link{simulate.glmb}}, 
#' \code{\link{extractAIC.glmb}}, \code{\link{dummy.coef.glmb}} and methods(class="glmb") for \code{glmb} 
#' and the methods and generic functions for classes \code{glm} and \code{lm} from which class \code{glmb} inherits.
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
#' @example inst/examples/Ex_rglmb.R
#' @export 
#' @rdname rglmb
#' @order 1



rglmb<-function(n=1,y,x,family=gaussian(),pfamily,offset=NULL,weights=1,
                use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE){
  
  ## Pull in information from the pfamily  
  pf=pfamily$pfamily
  okfamilies=pfamily$okfamilies  
  plinks=pfamily$plinks
  prior_list=pfamily$prior_list 
  simfun=pfamily$simfun
  
  ## Pull in information on families  
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  ## Check that the family is implemented for the pfamily
  
  if(family$family %in% okfamilies){
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
  
  
  ## Call relevant simulation function (for now without control2 list)
  
#  outlist=simfun(n=n,y=y,x=x,prior_list=prior_list,offset=offset,weights=weights,family=family)

  outlist = simfun(n = n, y = y, x = x, prior_list = prior_list,offset = offset, weights = weights, family = family, use_parallel = use_parallel, use_opencl = use_opencl, verbose = verbose)
  
  outlist$pfamily=pfamily
  
  return(outlist)
  
}


#' @rdname rglmb
#' @order 2
#' @method print rglmb
#' @export

print.rglmb<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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





