rglmb<-function(n=1,y,x,mu,P,wt=1,dispersion=NULL,nu=NULL,V=NULL,family=gaussian(),offset2=rep(0,nobs),start=NULL,Gridtype=2)UseMethod("rglmb")

rglmb.default_Old<-function(n=1,y,x,mu,P,wt=1,dispersion=NULL,family=gaussian(),offset2=rep(0,nobs),start=NULL,Gridtype=2
)
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
	wt2<-wt/dispersion2


	mu<-as.matrix(mu)
	Prior<-list(mean=as.numeric(mu),Precision=P)

	# Produce Maximum a posteriori estimates and get outputs)

	alpha<-x%*%mu+offset2  # Offset
  
  
  
	if(family$family=="gaussian")
		{
		x2<-x*sqrt(wt2)
		y1<-y-offset2
		y2<-y1*sqrt(wt2)

		k = length(mu)
	    	RA = chol(P)
    		W = rbind(x2, RA)
    		z = c(y2, as.vector(RA %*% mu))
    		IR = backsolve(chol(crossprod(W)), diag(k))
		b2<-crossprod(t(IR)) %*% crossprod(W, z) 

		out<-matrix(0,nrow=n,ncol=k)	
		draws<-matrix(1,n)
		LL<-matrix(1,n)

		for(i in 1:n){
    		out[i,]<-b2 + IR %*% rnorm(k)
		LL[i]<-f1(b=out[i,],y=y,x=x,alpha=offset2,wt=wt2)	
		}
		
#		outlist<-list(coefficients=out,PostMode=b2,Prior=Prior,n=n,iters=draws,famfunc=famfunc,Envelope=NULL,dispersion=dispersion2)
    outlist<-list(coefficients=out,PostMode=b2,Prior=Prior,iters=draws,famfunc=famfunc,Envelope=NULL,dispersion=dispersion2,loglike=LL)
		colnames(outlist$coefficients)<-colnames(x)
		

		}	
	
	else{

#	  start.timeE<-Sys.time()
	  
  
    
	opt1<-optim(par=start-mu, fn=f2, gr = f3,y=y,x=x,mu=mu-mu,P=P,
	alpha=alpha,wt=wt2,method="BFGS", hessian = TRUE)

	
  
	min1<-opt1$value
#	grad1<-f3(b1,y,x,mu-mu,P=P,alpha=alpha,wt=wt2)
	conver1<-opt1$convergence
	if(conver1>0) stop("Posterior Optimization failed")
	

  b1<-as.matrix(opt1$par)
  A1<-opt1$hessian	
	

	# Replace With solution from optim (may be more stable)

	A2<-A1
	b2<-b1	

	# Standardize Model Again to have diagonal variance-covariance matrix at 
	# posterior mode

	eS1 <- eigen(A2, symmetric = TRUE)
    	ev1 <- eS1$values
	
  if(length(eS1$values)==1){D1<-matrix(eS1$values,1)
        
  }
  else  D1<-diag(eS1$values)
  
  L2<-sqrt(D1)%*%t(eS1$vectors)
	L2Inv<-eS1$vectors%*%sqrt(solve(D1))
	b3<-L2%*%b2
	mu3<-L2%*%(mu-mu)
	x3<-x%*%L2Inv
	P3<-t(L2Inv)%*%P%*%L2Inv

	# Find diagonal matrix that has "smaller" precision than prior	
	# Follows Definition 3, and procedure on p. 1150 in Nygren 
	# Puts model into standard form 
	
  if(length(P3)==1) {
    P3Diag<-P3
    epsilon<-P3
    
  }
  else{
  
	P3Diag<-diag(diag(P3))
	epsilon<-P3Diag
  }
  
  
  scale<-1
	check<-0

	while(check==0){
		epsilon<-scale*P3Diag
		P4<-P3-epsilon				
		eStemp <- eigen(P4, symmetric = TRUE)
    		evtemp <- min(eStemp$values)

		if(evtemp>0) check<-1
		else scale<-scale/2
		}

	check2<-0
	scale2<-scale
	while(check2==0){
		scale<-scale+(scale2/10)
		epsilon_temp<-scale*P3Diag
		P4_temp<-P3-epsilon_temp
		eStemp <- eigen(P4_temp, symmetric = TRUE)
    		evtemp <- min(eStemp$values)

		if(evtemp>0){
						epsilon<-epsilon_temp
						P4<-P4_temp	
						}		
		else{check2<-1}
		}

	A3<-diag(nvars)-epsilon	

	# Standardize Again to Standard form

	eS2 <- eigen(epsilon, symmetric = TRUE)
    	ev2 <- eS2$value
  if(length(eS2$values)==1){D2<-eS2$values}
  else{D2<-diag(eS2$values)}
  
  
	L3<-sqrt(D2)%*%t(eS2$vectors)
	L3Inv<-eS2$vectors%*%sqrt(solve(D2))
	b4<-L3%*%b3
	mu4<-L3%*%mu3
	x4<-x3%*%L3Inv
	A4<-t(L3Inv)%*%A3%*%L3Inv
	P5<-t(L3Inv)%*%P4%*%L3Inv
	P6Temp<-P5+diag(nvars)
#	epsilontemp<-t(L3Inv)%*%epsilon%*%L3Inv



	opt3<-optim(par=b4, fn=f2, gr = f3,y=y,x=x4,mu=mu4-mu4,P=P6Temp,
	alpha=alpha,wt=wt2,method="BFGS", hessian = TRUE)


#	min4<-opt3$value
	b4temp<-opt3$par
#	conver3<-opt3$convergence
	Atemp<-opt3$hessian	
	atemp<-sqrt(diag(Atemp))


#  minlist<-list(min=min1,min4=min4)
#	valcheck<-list(b1=b1,b2=b2)
#	hcheck<-list(A1=A1,A2=A2)
#	gradcheck<-list(grad1=grad1)
#	ctest<-t(f3(b=b4,y=y,x=x4,mu=mu4,P=P5,alpha=alpha,wt=wt2))
	ctest2<-t(f3(b=b4temp,y=y,x=x4,mu=mu4,P=P5,alpha=alpha,wt=wt2))

  
	
	
	# Create enveloping function

	
#	Envelope<-glmbenvelope(bStar=b4temp,A=A4,f2=f2,f3=f3,f5=f5,f6=f6,y=y,x=x4,mu=mu4,P=P5,alpha=alpha,wt=wt2,Gridtype=Gridtype,n=n)
	
  if(n==1){
  Envelope<-glmbenvelope_c(bStar=b4temp, A=A4,y=y, x=x4,mu=mu4,P=P5,alpha=alpha,wt=wt2,family=family$family,link=family$link,Gridtype=Gridtype, n=n,sortgrid=FALSE)
  
  }
  if(n>1){
  Envelope<-glmbenvelope_c(bStar=b4temp, A=A4,y=y, x=x4,mu=mu4,P=P5,alpha=alpha,wt=wt2,family=family$family,link=family$link,Gridtype=Gridtype, n=n,sortgrid=TRUE)
  }

	
	
	qa1<-max(abs(b4-b4temp)*atemp)
	qa2<-max(abs(ctest2+t(b4temp))*atemp)
	if(qa1>0.01) {warning("Possible numeric imprecision during transformations type 1. \n Consider a stronger prior.")}
	if(qa2>0.01) {warning("Possible numeric imprecision during transformations type 2. \n Consider a stronger prior.")}

#	optlist<-list(b4=b4,b4temp=as.matrix(b4temp),ctest=as.matrix(t(ctest)),ctest2=as.matrix(t(ctest2)),qa1=qa1,qa2=qa2)

	# Simulate
# New glmbsim (Cpp)

#start.time<-Sys.time()

sim<-glmbsim_cpp(n=n,y=y,x=x4,mu=mu4,P=P5,alpha=alpha,wt=wt2,f2=f2,Envelope=Envelope,family=family$family,link=family$link)
#end.time<-Sys.time()
#time.taken<-end.time-start.time
#print("Time for new glmbsim")
#print(time.taken)

	# Undo Standardization

	# Warning-must verify this is correct (Drop mu could cause issues)

	out<-t(drop(mu)+L2Inv%*%L3Inv%*%t(sim$out))

  # Warning: Must verify this is correct 
  
	LL<-matrix(1,n)
	
	for(i in 1:n){
		  LL[i]<-f1(b=out[i,],y=y,x=x,alpha=offset2,wt=wt2)	
	}
	
  

#  outlist<-list(coefficients=out,PostMode=b2+mu,Prior=Prior,iters=sim$draws,famfunc=famfunc,Envelope=Envelope,dispersion=dispersion2,loglike=LL,b1=b1,A1=A1)
  outlist<-list(coefficients=out,PostMode=b2+mu,Prior=Prior,iters=sim$draws,famfunc=famfunc,Envelope=Envelope,dispersion=dispersion2,loglike=LL)
  colnames(outlist$coefficients)<-colnames(x)

	}

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


rglmb.default<-function(n=1,y,x,mu,P,wt=1,dispersion=NULL,nu=NULL,V=NULL,family=gaussian(),offset2=rep(0,nobs),start=NULL,Gridtype=2
)
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
    if(dispersion2>0){outlist<-glmbsim_Gauss_cpp(n=n,y=y,x=x,mu=mu,P=P,offset2=offset2,wt=wt,dispersion=dispersion,famfunc=famfunc,f1=f1,f2=f2,f3=f3,start=mu)
    }

    else{
         
      
         sim<-rmultireg(n=n,Y=matrix((y-offset2)*sqrt(wt),nrow=length(y)),X=x*sqrt(wt),Bbar=mu,A=P,nu=nu,V=V)
         draws<-matrix(1,n)
         LL<-matrix(1,n)
         
         for(i in 1:n){
            LL[i]<-f1(b=sim$B[i,],y=y,x=x,alpha=offset2,wt=wt/sim$Sigma[i])	
         }
         
         outlist<-list(coefficients=sim$B,PostMode=sim$BStar,
                       Prior=list(mean=as.numeric(mu),Precision=P),
                       iters=draws,famfunc=famfunc,Envelope=NULL,
                       dispersion=sim$Sigma,loglike=LL)
         
        }
    
  }
    
    
    
    
  else{
 #   dispersion2=dispersion
    if(is.null(dispersion)){dispersion2=1}
    outlist<-glmbsim_NGauss_cpp(n=n,y=y,x=x,mu=mu,P=P,offset2=offset2,wt=wt,dispersion=dispersion2,famfunc=famfunc,f1=f1,f2=f2,f3=f3,
    start=mu,family=family$family,link=family$link)
    
  }
  
  colnames(outlist$coefficients)<-colnames(x)
  
  outlist$call<-match.call()
  
  class(outlist)<-c(outlist$class,"rglmb")
  outlist
  
}






