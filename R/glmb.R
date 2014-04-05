glmb<-function(n,formula, family = binomial,dispersion=NULL,mu,Sigma,Gridtype=1, data, weights, subset2, 
    na.action, start = NULL, etastart, mustart, offset, control = list(...), 
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
    ...) UseMethod("glmb")



glmb.default<-function (n,formula, family = binomial,dispersion=NULL,mu,Sigma,Gridtype=1, data, weights, subset2, 
    na.action, start = NULL, etastart, mustart, offset, control = list(...), 
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
    ...) 
{
   # Note, renamed subset argument to subset2 as it caused conflict with subset function

    call <- match.call()
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
    m <- match(c("formula", "data", "subset2", "weights", "na.action", 
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
    
    # Initialize
	
    y<-fit$y	
    x<-fit$x
    b<-fit$coefficients	
    P<-solve(Sigma) 
    wtin<-fit$prior.weights	


      sim<-rglmb(n=n,y=y,x=x,mu=mu,P=P,wt=wtin,dispersion=dispersion,offset2=offset,family=family,start=b,Gridtype=Gridtype)
		
	dispersion2<-sim$dispersion
	famfunc<-sim$famfunc
	

	Prior<-list(mean=as.numeric(mu),Variance=Sigma)
	names(Prior$mean)<-colnames(fit$x)
	colnames(Prior$Variance)<-colnames(fit$x)
	rownames(Prior$Variance)<-colnames(fit$x)
#	DICinfo<-glmbdic(sim$coefficients,y=y,x=x,f1=famfunc$f1,f4=famfunc$f4,wt=wtin)
	DICinfo<-glmbdic(sim$coefficients,y=y,x=x,f1=famfunc$f1,f4=famfunc$f4,wt=wtin/dispersion2,dispersion2)



if (!is.null(offset)) {linear.predictors<-t(offset+x%*%t(sim$coefficients))}
if (is.null(offset)) {linear.predictors<-t(x%*%t(sim$coefficients))}
	linkinv<-fit$family$linkinv
	fitted.values<-linkinv(linear.predictors)

	outlist<-list(glm=fit,coefficients=sim$coefficients,coef.means=colMeans(sim$coefficients),
      coef.mode=sim$PostMode,Prior=Prior,
	fitted.values=fitted.values,
	linear.predictors=linear.predictors,
	deviance=DICinfo$Deviance,
	pD=DICinfo$pD,
	Dbar=DICinfo$Dbar,Dthetabar=DICinfo$Dthetabar,
	DIC=DICinfo$DIC,famfunc=famfunc,iters=sim$iters,dispersion=dispersion)

	outlist$call<-match.call()

	class(outlist)<-c(outlist$class,"glmb")
	outlist
}




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



summary.glmb<-function(object,...){

mres<-colMeans(residuals(object))


l1<-length(object$coef.means)
n<-length(object$coefficients[,1])
percentiles<-matrix(0,nrow=l1,ncol=7)
se<-sqrt(diag(var(object$coefficients)))
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

glmsummary<-summary(object$glm)
se1<-sqrt(diag(glmsummary$cov.scaled))


Tab1<-cbind("Prior Mean"=object$Prior$mean,"Prior.sd"=as.numeric(sqrt(diag(object$Prior$Variance)))
,"Max Like."=as.numeric(object$glm$coefficients),"Like.sd"=se1)
TAB<-cbind("Post.Mode"=as.numeric(object$coef.mode),"Post.Mean"=object$coef.means,"Post.Sd"=se,"MC Error"=as.numeric(mc),"Pr(tail)"=as.numeric(pval2))
TAB2<-cbind("1.0%"=percentiles[,1],"2.5%"=percentiles[,2],"5.0%"=percentiles[,3],Median=as.numeric(percentiles[,4]),"95.0%"=percentiles[,5],"97.5%"=as.numeric(percentiles[,6]),"99.0%"=as.numeric(percentiles[,7]))

rownames(TAB2)<-rownames(TAB)

res<-list(call=object$call,n=n,residuals=mres,coefficients1=Tab1,coefficients=TAB,Percentiles=TAB2,pD=object$pD,deviance=object$deviance,DIC=object$DIC,iters=mean(object$iters))

class(res)<-"summary.glmb"

res

}



print.summary.glmb<-function(x,...){
cat("Call\n")
print(x$call)
cat("\nExpected Deviance Residuals:\n")
print(x$residuals,digits=4)
cat("\nPrior and Maximum Likelihood Estimates with Standard Deviations\n\n")
printCoefmat(x$coefficients1,digits=4)
cat("\nBayesian Estimates Based on",x$n,"iid draws\n\n")
printCoefmat(x$coefficients,digits=4,P.values=TRUE,has.Pvalue=TRUE)
cat("\nDistribution Percentiles\n\n")
printCoefmat(x$Percentiles,digits=4)
        cat("\nEffective Number of Parameters:",x$pD,"\n")
        cat("Expected Residual Deviance:",mean(x$deviance),"\n")
        cat("DIC:",x$DIC,"\n\n")
        cat("Mean Likelihood Subgradient Candidates Per iid sample:",x$iters,"\n\n")


}

residuals.glmb<-function(object,...)
{
	y<-object$glm$y	
	n<-length(object$coefficients[,1])
	fitted.values<-object$fitted.values

	dev.residuals<-object$glm$family$dev.resids
	DevRes<-matrix(0,nrow=n,ncol=length(y))
#	yrep<-matrix(0,nrow=n,ncol=length(y))


	for(i in 1:n)
	{
	DevRes[i,]<-sign(y-fitted.values[i,])*sqrt(dev.residuals(y,fitted.values[i,],1))
#	if(fit$family$family=="poisson")yrep[i,]<-rpois(length(y),fitted.values[i,])

	}

	colnames(DevRes)<-names(y)
	DevRes
}


extractAIC.glmb<-function(fit,scale=NULL,k=NULL,...)
{
c(fit$pD,fit$DIC)

}

logLik.glmb<-function(object,...){

y<-object$glm$y
x<-object$glm$x
wt<-object$glm$prior.weights
dispersion<-object$dispersion


#alpha # to be updated 
f1<-object$famfunc$f1
n<-length(object$coefficients[,1])

logLikout<-matrix(0,nrow=n,ncol=1)
#f1temp(b=out$coefficients,y=out$y,x=out$x,alpha=0,wt=Claims/dispersion)
for(i in 1:n){
logLikout[i,1]<--f1(object$coefficients[i,],y=y,x=x,wt=wt/dispersion)
}
logLikout
}
