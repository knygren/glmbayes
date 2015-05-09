
glmbfamfunc<-function(family) UseMethod("glmbfamfunc")

glmbfamfuncset<-function(family){

# need to add handling for offsets 

# f1-(negative)Log-Likelihood function
# f2-(negative)Log-Posterior density
# f3-(negative) Gradient for Log-Posterior density 
# f4-Deviance function (note multiplies the difference in two times log-likelihood by dispersion)

if(family$family=="gaussian")
	{
	f1<-function(b,y,x,alpha=0,wt=1){
			Xb<-alpha+x%*%b
			-sum(dnorm(y, mean=Xb,sd=sqrt(1/wt),log=TRUE))
		}
	f2<-function(b,y,x,mu,P,alpha=0,wt=1){
			Xb<-alpha+x%*%b
			-sum(dnorm(y, mean=Xb,sd=sqrt(1/wt),log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
		}
	f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      l2<-length(y)
      ltemp<-length(wt)
      yxb2<-NULL
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
      Ptemp<-diag(wt)      
      }
      xb<-alpha+x%*%b-y
			yXb2<-P%*%(b-mu)-t(x)%*%(Ptemp%*%xb)
      yxb2
		}
	f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
		dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dnorm(y, mean=y,sd=sqrt(1/wt),log=TRUE)))
				}

  # Edit This
	f5<-f2_gaussian
	f6<-f3_gaussian
	
  
  
	}

# Check if Poisson weight should be outside or inside function


if(family$family=="poisson"||family$family=="quasipoisson")
	{
	f1<-function(b,y,x,alpha=0,wt=1){
			lambda<-t(exp(alpha+x%*%b))
-sum(dpois(y, lambda,log=TRUE)*wt)
		}
	f2<-function(b,y,x,mu,P,alpha=0,wt=1){
			lambda<-t(exp(alpha+x%*%b))
			-sum(dpois(y, lambda,log=TRUE)*wt)+0.5*t((b-mu))%*%P%*%(b-mu)
		}
	f3<-function(b,y,x,mu,P,alpha=0,wt=1){
		-t(x)%*%((y-exp(alpha+x%*%b))*wt)+P%*%(b-mu)
		
	}
	f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
		dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dpois(y, y,log=TRUE)*wt))
				}
  f5<-f2_poisson
  f6<-f3_poisson

	}

if(family$family %in%  c("binomial","quasibinomial") && family$link=="logit")
	{
	f1<-function(b,y,x,alpha=0,wt=1){
			lambda<-t(exp(alpha+x%*%b))
			p<-lambda/(1+lambda)
-sum(dbinom(round(wt*y),round(wt), p,log=TRUE))
		}
	f2<-function(b,y,x,mu,P,alpha=0,wt=1){
			lambda<-t(exp(alpha+x%*%b))
			p<-lambda/(1+lambda)
-sum(dbinom(round(wt*y),round(wt),p,log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
		}
	f3<-function(b,y,x,mu,P,alpha=0,wt=1){
			p<-1/(1+t(exp(-alpha-x%*%b)))
		t(x)%*%((t(p)-y)*wt)+P%*%(b-mu)
		}
	f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
		dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dbinom(round(wt*y),round(wt),y,log=TRUE)))
		}
  f5<-f2_binomial_logit
  f6<-f3_binomial_logit
}

if(family$family %in%  c("binomial","quasibinomial") && family$link=="probit")
	{
	f1<-function(b,y,x,alpha=0,wt=1){
			p<-pnorm(alpha+x%*%b)
			-sum(dbinom(round(wt*y),round(wt), p,log=TRUE))
			}
	f2<-function(b,y,x,mu,P,alpha=0,wt=1){
			p<-pnorm(alpha+x%*%b)
-sum(dbinom(round(wt*y),round(wt),p,log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
		}


	f3<-function(b,y,x,mu,P,alpha=0,wt=1){
		p1<-pnorm(alpha+x%*%b)
		p2<-pnorm(-alpha-x%*%b)
-t(x)%*%as.matrix(((y*dnorm(alpha+x%*%b)/p1)-(1-y)*dnorm(alpha+x%*%b)/p2)*wt)+P%*%(b-mu)
		}

	f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
		dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dbinom(round(wt*y),round(wt),y,log=TRUE)))
		}
  f5<-f2_binomial_probit
  f6<-f3_binomial_probit

	}
if(family$family %in%  c("binomial","quasibinomial") && family$link=="cloglog")
	{
	f1<-function(b,y,x,alpha=0,wt=1){
			Xb<-alpha+x%*%b
			p<-1-exp(-exp(Xb))
			-sum(dbinom(round(wt*y),round(wt), p,log=TRUE))
			}
	f2<-function(b,y,x,mu,P,alpha=0,wt=1){
			Xb<-alpha+x%*%b
			p<-1-exp(-exp(Xb))
-sum(dbinom(round(wt*y),round(wt),p,log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
		}
	f3<-function(b,y,x,mu,P,alpha=0,wt=1){
			Xb<-alpha+x%*%b
			p1<-1-exp(-exp(Xb))
			p2<-exp(-exp(Xb))
			atemp<-exp(Xb-exp(Xb))
-t(x)%*%as.matrix(((y*atemp/p1)-(1-y)*atemp/p2)*wt)+P%*%(b-mu)
		}
	f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
		dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dbinom(round(wt*y),round(wt),y,log=TRUE)))
		}
  f5<-f2_binomial_cloglog
  f6<-f3_binomial_cloglog

}
if(family$family=="Gamma" && family$link=="log")
	{

	f1<-function(b,y,x,alpha=0,wt=1){
	mu<-t(exp(alpha+x%*%b))
	disp2<-1/wt
	-sum(dgamma(y,shape=1/disp2,scale=mu*disp2,log=TRUE))
	}

#  f1<-f1_gamma
	f2<-function(b,y,x,mu,P,alpha=0,wt=1){
	mu2<-t(exp(alpha+x%*%b))
	disp2<-1/wt
		-sum(dgamma(y,shape=1/disp2,scale=mu2*disp2,log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
	}

  
	f3<-function(b,y,x,mu,P,alpha=0,wt=1){
			mu2<-t(exp(alpha+x%*%b))
		t(x)%*%(t(1-y/mu2)*wt)+P%*%(b-mu)
		}

	f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
	disp2<-1/wt

	dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dgamma(y,shape=1/disp2,scale=y*disp2,log=TRUE)))

		}

  f5<-f2_gamma
  f6<-f3_gamma


	}	

	list(f1=f1,f2=f2,f3=f3,f4=f4,f5=f5,f6=f6)
	
}



glmbfamfunc.default<-function(family)
{

out<-glmbfamfuncset(family=family)
out$call<-match.call()

class(out)<-"glmbfamfunc"
out
}

print.glmbfamfunc<-function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nNegative Log-Likelihood Function:\nf1<-")
  print(x$f1)
  cat("\n(Negative Log-Posterior Function:\nf2<-")
  print(x$f2)
  cat("\nNegative Log-Posterior Gradient Function:\nf3<-")
  print(x$f3)
  cat("\nDeviance Function:\nf4<-")
  print(x$f4)
 	

}


