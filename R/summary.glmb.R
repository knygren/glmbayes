#' Summarizing Bayesian Generalized Linear Model Fits
#'
#' These functions are all \code{\link{methods}} for class \code{glmb} or \code{summary.glmb} objects.
#' 
#' @aliases 
#' summary.glmb
#' print.summary.glmb
#' @param object an object of class \code{"glmb"} for which a summary is desired.
#' @param x an object of class \code{"summary.glmb"} for which a printed output is desired.
#' @param digits the number of significant digits to use when printing.
#' @param \ldots Additional optional arguments
#' @return \code{summary.glmb} returns a object of class \code{"summary.glmb"}, a 
#' list with components: 
#' \item{call}{the component from \code{object}}
#' \item{n}{number of draws generated}
#' \item{residuals}{vector of mean deviance residuals}
#' \item{coefficients1}{Matrix with the prior mean and maximum likelihood coefficients with associated standard deviations}
#' \item{coefficients}{Matrix with columns for the posterior mode, posterior mean, posterior standard
#' deviation, monte carlo error, and tail probabilities (posterior probability of observing a 
#' value for the coefficient as extreme as the prior mean)}
#' \item{Percentiles}{Matrix with estimated percentiles associated with the posterior density}
#' \item{pD}{Estimated effective number of parameters}
#' \item{deviance}{Vector with draws for the deviance}
#' \item{DIC}{Estimated DIC statistic}
#' \item{iters}{Average number of candidates per generated draws}
#' @details The \code{summary.glmb} function summarizes the output from the \code{glmb} function.
#' Key output includes mean residuals, information related to the prior, mean coefficients 
#' with associated stats, percentiles for the coefficients, as well as the effective number of
#' parameters and the DIC statistic. 
#' 
#' @seealso \code{\link{lmb}}, \code{\link{glmb}}, \code{\link{summary}}, \code{[stats]\link{summary.lm}},\code{[stats]\link{summary.glm}}. 
#' @example inst/examples/Ex_summary.glmb.R
#' @export
#' @method summary glmb


summary.glmb<-function(object,...){
  
  
  #est.disp <- FALSE
  
  #res<-residuals(object)
  #nobs=length(object$y)
  #n=nrow(object$coefficients)
  #df.r <- length(object$y)-object$pD
  
  dispersion=get_dispersion(object)
  #if(!is.null(object$DIC)) DIC=object$DIC
  #else DIC=NA
  #print("Mean Dispersion")
    mres<-colMeans(residuals(object))

  l1<-length(object$coef.means)
  n<-length(object$coefficients[,1])
  percentiles<-matrix(0,nrow=l1,ncol=7)
  se<-sqrt(diag(var(object$coefficients)))
#  if(object$family$family=="quasipoisson") se=se*sqrt(mean(dispersion))
  mc<-se/object$n
  mc<-se/n
  priorrank<-matrix(0,nrow=l1,ncol=1)
  pval1<-matrix(0,nrow=l1,ncol=1)
  pval2<-matrix(0,nrow=l1,ncol=1)
  
  for(i in 1:l1){
    percentiles[i,]<-quantile(object$coefficients[,i],probs=c(0.01,0.025,0.05,0.5,0.95,0.975,0.99))
    test<-append(object$coefficients[,i],object$Prior$mean[i])
    test2<-rank(test)
    priorrank[i,1]<-test2[n+1]
    priorrank[i,1]<-test2[n+1]
    pval1[i,1]<-priorrank[i,1]/(n+1)
    pval2[i,1]<-min(pval1[i,1],1-pval1[i,1])
    
  }

  ## Note: This restricts function to classes glmb and lmb
  ## Could break calls to this function from other classes
  
  if(!is.null(object$glm)){
  glmsummary<-summary(object$glm)
  se1<-sqrt(diag(glmsummary$cov.scaled))
  ml=as.numeric(object$glm$coefficients)
  } 

  
  
  if(!is.null(object$lm)){
    lm_summary=summary.lm(object$lm)
    se1=lm_summary$coefficients[,2]
    ml=lm_summary$coefficients[,1]
  }  
    
  Tab1<-cbind("Prior Mean"=object$Prior$mean,"Prior.sd"=as.numeric(sqrt(diag(object$Prior$Variance)))
              ,"Max Like."=ml,"Like.sd"=se1)
  TAB<-cbind("Post.Mode"=as.numeric(object$coef.mode),"Post.Mean"=object$coef.means,"Post.Sd"=se,"MC Error"=as.numeric(mc),"Pr(tail)"=as.numeric(pval2))
  TAB2<-cbind("1.0%"=percentiles[,1],"2.5%"=percentiles[,2],"5.0%"=percentiles[,3],Median=as.numeric(percentiles[,4]),"95.0%"=percentiles[,5],"97.5%"=as.numeric(percentiles[,6]),"99.0%"=as.numeric(percentiles[,7]))
  
  rownames(TAB2)<-rownames(TAB)
  
  res<-list(
            call=object$call,
            n=n,
            residuals=mres,
            coefficients1=Tab1,
            coefficients=TAB,
            Percentiles=TAB2,
            pD=object$pD,
            deviance=object$deviance,
            DIC=object$DIC,
  #          DIC=DIC,
            dispersion=mean(object$dispersion),
  #          dispersion=mean(dispersion),
            iters=mean(object$iters)
            )
  
  class(res)<-"summary.glmb"
  
  res
  
}

#' @rdname summary.glmb
#' @export
#' @method print summary.glmb

print.summary.glmb<-function(x,digits = max(3, getOption("digits") - 3),...){
  cat("Call\n")
  print(x$call)
  cat("\nExpected Deviance Residuals:\n")
  print(x$residuals,digits=digits)
  cat("\nPrior and Maximum Likelihood Estimates with Standard Deviations\n\n")
  printCoefmat(x$coefficients1,digits=digits)
  cat("\nBayesian Estimates Based on",x$n,"iid draws\n\n")
  printCoefmat(x$coefficients,digits=digits,P.values=TRUE,has.Pvalue=TRUE)
  cat("\nDistribution Percentiles\n\n")
  printCoefmat(x$Percentiles,digits=digits)
  cat("\nEffective Number of Parameters:",x$pD,"\n")
  cat("Expected Residual Deviance:",mean(x$deviance),"\n")
  cat("DIC:",x$DIC,"\n\n")
  cat("Expected Mean dispersion:",x$dispersion,"\n")
  cat("Sq.root of Expected Mean dispersion:",sqrt(x$dispersion),"\n\n")
  cat("Mean Likelihood Subgradient Candidates Per iid sample:",x$iters,"\n\n")
  
  
}


get_dispersion<-function(object){
  
  df.r <- length(object$y)-object$pD
  
  
  if(!is.null(object$dispersion)) object$dispersion
  else{
    
    if(object$family$family=="quasipoisson"){
      n=nrow(object$coefficients)
      disp_temp=rep(0,n)
      m=length(object$y)  
      k=ncol(object$x)    
      res_temp=matrix(0,nrow=n,ncol=m)
      fit_temp=object$x%*%t(object$coefficients)
      for(l in 1:n){
        if(is.null(object$fit$offset))         fit_temp[1:m,l]=exp(fit_temp[1:m,l])
        else fit_temp[1:m,l]=exp(object$fit$offset+fit_temp[1:m,l])
        
        res_temp[l,1:m]=(object$y-fit_temp[1:m,l])
        disp_temp[l]=(1/(m-k))*sum(res_temp[l,1:m]^2*object$prior.weights/fit_temp[1:m,l])
  
      }
    
     return(disp_temp) 
    }
    
  }
  
}

