#' The Bayesian Normal Regression Distribution
#'
#' \code{rnorm_reg} is used to generate iid samples from Bayesian linear models with a normal prior. 
#' The model is specified by providing a data vector, a design matrix, and 2 prior constants.
#' @aliases
#' rnorm_reg
#' rnorm_reg_cpp
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x a design matrix of dimension \code{m * p}.
#' @param mu a vector of length \code{p} giving the prior means of the variables in the design matrix.
#' @param P a positive-definite symmetric matrix of dimension \code{p * p} specifying the prior precision matrix of the variable.
#' @param offset2 this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.
#' @param wt an optional vector of \sQuote{prior weights} to be used in the fitting process. Should be NULL or a numeric vector.
#' @param dispersion the dispersion parameter. Either a single numerical value or NULL (the default). Must be provided here. Use \code{\link{rnorm_gamma_reg}} or Block-Gibbs sampling using \code{\link{rglmb_dispersion}} to give the dispersion a prior.
#' @return \code{rnorm_reg} returns a object of class \code{"rglmb"}.  The function \code{summary} 
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
#' cloglog link functions).  The function returns the simulated Bayesian coefficients 
#' and some associated outputs.
#' 
#' For the gaussian family, iid samples from the posterior density is genererated using 
#' standard simulation procedures for multivariate normal densities. For all other 
#' families, the samples are generated using accept-reject procedures by leveraging the 
#' likelihood subgradient approach of Nygren and Nygren (2006). This approach relies on
#' tight enveloping functions that bound the posterior density from above. The Gridtype 
#' parameter is used to select the method used for determining the size of a Grid used 
#' to build the enveloping function. See \code{\link{EnvelopeBuild_c}} for details. 
#' Depending on the selection, the time to build the envelope and the acceptance rate 
#' during the simulation process may vary. The returned value \code{iters} contains the 
#' number of candidates generated before acceptance for each draw.
#' @examples
#' 1+1
#' 10+1


rnorm_reg<-function(n=1,y,x,mu,P,wt=1,dispersion=1,offset2=NULL)
  {
  
  if(is.numeric(n)==FALSE||is.numeric(y)==FALSE||is.numeric(x)==FALSE||
     is.numeric(mu)==FALSE||is.numeric(P)==FALSE) stop("non-numeric argument to numeric function")
  
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
  
  okfamilies <- c("gaussian")
  if(family$family %in% okfamilies){
    if(family$family=="gaussian") oklinks<-c("identity")
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
  
  

    dispersion2=dispersion
    if(is.null(dispersion)){dispersion2=0}
    if(dispersion2>0){
      
      outlist<-rnorm_reg_cpp(n=n,y=y,x=x,mu=mu,P=P,offset2=offset2,wt=wt,dispersion=dispersion,famfunc=famfunc,f1=f1,f2=f2,f3=f3,start=mu)
    }
    
    else{
      
      stop("dispersion must be>0")
      
    }
    


  colnames(outlist$coefficients)<-colnames(x)
  
  outlist$call<-match.call()
  
  class(outlist)<-c(outlist$class,"rglmb")
  outlist
  
}





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


summary.rglmb<-function(object,...){


n<-length(object$coefficients[,1])  
l1<-length(object$PostMode)
percentiles<-matrix(0,nrow=l1,ncol=7)
se<-sqrt(diag(var(object$coefficients)))
mc<-se/n
Priorwt<-(se/sqrt(diag(solve(object$Prior$Precision))))^2
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

Tab1<-cbind("Prior Mean"=object$Prior$mean,"Prior.sd"=as.numeric(sqrt(diag(solve(object$Prior$Precision)))),"Approx.Prior.wt"=Priorwt)
TAB<-cbind("Post.Mode"=as.numeric(object$PostMode),"Post.Mean"=colMeans(coef(object)),"Post.Sd"=se,"MC Error"=as.numeric(mc),"Pr(tail)"=as.numeric(pval2))
TAB2<-cbind("1.0%"=percentiles[,1],"2.5%"=percentiles[,2],"5.0%"=percentiles[,3],Median=as.numeric(percentiles[,4]),"95.0%"=percentiles[,5],"97.5%"=as.numeric(percentiles[,6]),"99.0%"=as.numeric(percentiles[,7]))

rownames(Tab1)<-rownames(TAB)
rownames(TAB2)<-rownames(TAB)

res<-list(call=object$call,n=n,coefficients1=Tab1,coefficients=TAB,Percentiles=TAB2)

class(res)<-"summary.rglmb"
res

}

print.summary.rglmb<-function(x,...){
cat("Call\n")
print(x$call)
cat("\nPrior Estimates with Standard Deviations\n\n")
printCoefmat(x$coefficients1,digits=4)
cat("\nBayesian Estimates Based on",x$n,"iid draws\n\n")
printCoefmat(x$coefficients,digits=4,P.values=TRUE,has.Pvalue=TRUE)
cat("\nDistribution Percentiles\n\n")
printCoefmat(x$Percentiles,digits=4)

}










