#' The Bayesian Generalized Linear Model Distribution
#'
#' \code{rglmb} is used to generate iid samples for Bayesian Generalized Linear Models.
#' The model is specified by providing a data vector, a design matrix, 
#' the family (determining the likelihood function) and the pfamily (determining the 
#' prior distribution).
#' @param y a vector of observations of length \code{m}.
#' @param x for \code{rglmb} a design matrix of dimension \code{m * p} and for \code{print.rglmb} the object to be printed. 
#' @param pfamily a description of the prior distribution and associated constants to be used in the model. This
#' should be a pfamily function (see \code{\link{pfamily}} for details of pfamily functions.)
#' @param offset an offset parameter
#' @param weights a weighting variable
#' @inheritParams glmb
#' @return Currently mainly the draws for the dispersion and the regression coefficients
#' will be updated to return outputs consistent with other function
#' @family modelfuns
#' @seealso The classical modeling functions \code{\link[stats]{lm}} and \code{\link[stats]{glm}}.
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



rglmb<-function(n=1,y,x,family=gaussian(),pfamily,offset=NULL,weights=1){
  
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
  
  outlist=simfun(n=n,y=y,x=x,prior_list=prior_list,offset=offset,weights=weights,family=family)

  
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





