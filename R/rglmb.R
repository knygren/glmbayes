rglmb<-function(n=1,y,x,mu,P,wt=1,dispersion=NULL,nu=NULL,V=NULL,family=gaussian(),offset2=rep(0,nobs),start=NULL,Gridtype=3)
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
  
  okfamilies <- c("gaussian","poisson","binomial","quasipoisson","quasibinomial","Gamma")
  if(family$family %in% okfamilies){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-c("log")		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-c("logit","probit","cloglog")		
    if(family$family=="Gamma") oklinks<-c("log")		
    if(family$link %in% oklinks){
      
      famfunc<-glmbfamfunc(family)
      f1<-famfunc$f1
      f2<-famfunc$f2
      f3<-famfunc$f3
      f5<-famfunc$f5
      f6<-famfunc$f6
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
  
  
  if(family$family=="poisson"||family$family=="binomial")dispersion2<-1
  else dispersion2<-dispersion
  
  if(family$family=="gaussian"){
    dispersion2=dispersion
    if(is.null(dispersion)){dispersion2=0}
    if(dispersion2>0){
      
      outlist<-glmbsim_Gauss_cpp(n=n,y=y,x=x,mu=mu,P=P,offset2=offset2,wt=wt,dispersion=dispersion,famfunc=famfunc,f1=f1,f2=f2,f3=f3,start=mu)
    }
    
    else{
      
      outlist=rnorm_gamma_reg(n=n,y=y,x=x,mu=mu,P=P,nu=nu,V=V,offset2=offset2,wt=wt)
      
    }
    
  }
  
  else{
    if(is.null(dispersion)){dispersion2=1}

    # f1, f2, and f3 passed here - Likely legacy of R code
    ## Can eliminate and replace with calling of corresponding c++ functions
    outlist<-glmbsim_NGauss_cpp(n=n,y=y,x=x,mu=mu,P=P,offset2=offset2,wt=wt,dispersion=dispersion2,
                                famfunc=famfunc,f1=f1,f2=f2,f3=f3,
                                start=start,family=family$family,link=family$link,Gridtype=Gridtype)
    
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










