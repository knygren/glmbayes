#' The Bayesian Generalized Linear Model Dispersion Distribution
#'
#' \code{rGamma_reg} is used to generate iid samples for the dispersion parameter in Bayesian Generalized Linear models. 
#' The model is specified by providing the model structure, regression coefficient, and a Gamma prior.
#' @aliases
#' rGamma_reg
#' print.rGamma_reg
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x a design matrix of dimension \code{m * p}.
#' @param prior_list a list with the prior parameters (shape and rate) for the 
#' prior distribution and a vector with the assumed value (beta) for the regression coefficients.
#' @param offset this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.
#' @param weights an optional vector of \sQuote{prior weights} to be used in the fitting process. Should be NULL or a numeric vector.
#' @param family a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)
#' @param digits Number of significant digits to use for printed outpute
#' @param object an object of class \code{"glmb_dispersion"} that is to be summarized
#' @param \ldots further arguments passed to or from other methods
#' @return \code{rGamma_reg} returns a object of class \code{"rglmbdisp"}.  The function \code{summary} 
#' (i.e., \code{\link{summary.rglmb}}) can be used to obtain or print a summary of the results.
#' The generic accessor functions \code{\link{coefficients}}, \code{\link{fitted.values}},
#' \code{\link{residuals}}, and \code{\link{extractAIC}} can be used to extract
#' various useful features of the value returned by \code{\link{rGamma_reg}}.
#' An object of class \code{"rglmbdisp"} is a list containing at least the following components:
#' \item{dispersion}{an \code{n} by \code{1} matrix with simulated values for the dispersion}
#' \item{Prior}{A list with two components. The first being the prior mean vector and the second the prior precision matrix}
#' @details The \code{rGamma_reg} function produces iid samples for 
#' the dispersion parameter in Bayesian generalized linear models (gaussian and Gamma families
#' only). Core required inputs for the function include the data vector, the design  
#' matrix, an estimate for the regression coefficients, and a prior gamma distribution specification. 
#' The function returns the simulated Bayesian conditional dispersion and the prior specification.
#' The most practical use of this function is likely to be in conjunction with the 
#' \code{\link{rglmb}} function as part of block Gibbs sampling procedues.
#' 
#' For the gaussian family, iid samples from the posterior density are genererated using 
#' standard simulation procedures for gamma distributions. For the Gamma family,  
#' , the samples are generated using accept-reject procedures by leveraging the 
#' likelihood subgradient approach of Nygren and Nygren (2006). 
#' @family simfuncs 
#' @references A reference
#' @example inst/examples/Ex_rglmb_dispersion.R
#' @export 
#' @rdname rGamma_reg
#' @order 1
#' @export

## Offset in this model may need to be handled better


rGamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian()){
    
  call <- match.call()

  ## Renaming for consistency with earlier version
  
  wt=weights
  alpha=offset
  
  b=prior_list$beta
  shape=prior_list$shape
  rate=prior_list$rate
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  okfamilies <- c("gaussian","Gamma")
  if(family$family %in% okfamilies){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="Gamma") oklinks<-c("log")		
    if(family$link %in% oklinks)  {}
    else{stop(gettextf("link \"%s\" not available for selected family; available links are %s", 
                                                         family$link , paste(sQuote(oklinks), collapse = ", ")), 
                                                domain = NA)
    }
  }
  
  else{
    stop(gettextf("family \"%s\" not available in glmbdisp; available families are %s", 
                  family$family , paste(sQuote(okfamilies), collapse = ", ")), 
         domain = NA)
    
  }
  
  n1<-length(y)
  
  if(family$family=="gaussian"){
  y1<-as.matrix(y)-alpha
  xb<-x%*%b
  res<-y1-xb
  SS<-res*res
  
  a1<-shape+n1/2
  b1<-rate+sum(SS)/2
  
  out<-1/rgamma(n,shape=a1,rate=b1) 
  }

if(family$family=="Gamma")
  {
  
  mu1<-t(exp(alpha+x%*%b))
  
  testfunc<-function(v,wt){  
    -sum(lgamma(wt*v)+0.5*log(wt*v)+wt*v-wt*v*log(wt*v))
  }

  
  shape2=shape + 0.5 *n1
  rate1=rate +sum(wt*((y/mu1)-log(y/mu1)-1))
  
  vstar1<-shape2/rate1
  
  vout<-function(v){
    vstar1-(v/rate1)*sum((wt*digamma(wt*v) -wt*log(wt*v) + 0.5/v) )  
  }
  
  # Initialize vstar2
  vstar<-vstar1
  
  ## Optimize vstar?
  for(j in 1:20){
    vstar<-vout(vstar)
  }

  testbar<-testfunc(vstar,wt)
  cbar<--sum((wt*digamma(wt*vstar) -wt*log(wt*vstar) + 0.5/vstar))
  
  
  
  rate2=  rate +sum(wt*((y/mu1)-log(y/mu1)-1))-sum((wt*digamma(wt*vstar) -wt*log(wt*vstar) + 0.5/vstar) )
  
  out<-matrix(0,n)
  test<-matrix(0,n)
  a<-matrix(0,n)
  
  ## Implements rejection sampling for dispersion (likelihood subgradient approach)
  ## Likely should have a short paper with this derivation
  ## Not sure if approach extends to other densities besides gamma
  
  for(i in 1:n)
  {
    while(a[i]==0){
      out[i]<-rgamma(1,shape=shape2,rate=rate2)
      
      test[i]<-testfunc(out[i],wt)-(testbar+cbar*(out[i]-vstar))-log(runif(1,0,1))
      if(test[i]>0) a[i]<-1
    }
  }
  
  out<-1/out
  
}
  
  outlist=list(
    coefficients=matrix(b,nrow=1,ncol=length(b)),
    coef.mode=NULL,
    dispersion=out,
    Prior=list(shape=shape,rate=rate),
    prior.weights=weights,
    y=y,
    x=x,
    famfunc=famfunc(family),
    iters=rep(1,n),
    Envelope=NULL
    )
  

  outlist$call<-match.call()
  
  class(outlist)<-c(outlist$class,"rGamma_reg")

  return(outlist)

}


#' @export
#' @rdname rGamma_reg
#' @order 3
#' @method print rGamma_reg

print.rGamma_reg<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Simulated Dispersion")
    cat(":\n")
    print.default(format(x$dispersion, digits = digits), 
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
}

#' @export
#' @rdname rGamma_reg
#' @order 3
#' @method summary rGamma_reg


summary.rGamma_reg<-function(object,...){
  
  n<-length(object$dispersion)  
  percentiles<-matrix(0,nrow=1,ncol=7)
  me=mean(object$dispersion)
  se<-sqrt(var(object$dispersion))
  mc<-se/n
  Priorwt<-(se/(sqrt(object$Prior$shape)/object$Prior$rate))^2
    percentiles[1,]<-quantile(object$dispersion,probs=c(0.01,0.025,0.05,0.5,0.95,0.975,0.99))
    test<-append(object$dispersion,object$Prior$shape/object$Prior$rate)
    test2<-rank(test)
    priorrank<-test2[n+1]
   pval1<-priorrank/(n+1)
  pval2<-min(pval1,1-pval1)
    

  Tab1<-cbind("Prior.Mean"=object$Prior$shape/object$Prior$rate,"Prior.Sd"=sqrt(object$Prior$shape)/object$Prior$rate
              ,"Approx.Prior.wt"=Priorwt
              )
  TAB<-cbind(
    #"Post.Mode"=as.numeric(object$PostMode),
    "Post.Mean"=me,
    "Post.Sd"=se,
    "MC Error"=as.numeric(mc)
    ,"Pr(tail)"=as.numeric(pval2)
    )
  TAB2<-cbind("1.0%"=percentiles[,1],"2.5%"=percentiles[,2],"5.0%"=percentiles[,3],Median=as.numeric(percentiles[,4]),"95.0%"=percentiles[,5],"97.5%"=as.numeric(percentiles[,6]),"99.0%"=as.numeric(percentiles[,7]))
  
rownames(TAB)=c("dispersion")
rownames(Tab1)=c("dispersion")
rownames(TAB2)=c("dispersion")

  res<-list(call=object$call,
            n=n,
            coefficients1=Tab1,
            coefficients=TAB,
            Percentiles=TAB2
            )
  
  # Reuse summary.rglmb class
  
  class(res)<-"summary.rglmb"
  
  res
  
}


