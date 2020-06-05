
rnorm_gamma_reg<-function(n,y,x,mu,P,nu,V,offset2,wt){

## Should add dimension checks here  
## Should do away with need for passing f1 since it is a known function
## Should move core part of rmultireg inside this code to avoid call

famfunc=glmbfamfunc(gaussian())  
f1=famfunc$f1
    
sim<-rmultireg(n=n,Y=matrix((y-offset2)*sqrt(wt),nrow=length(y)),X=x*sqrt(wt),Bbar=mu,A=P,nu=nu,V=V)


draws<-matrix(1,n)
LL<-matrix(1,n)

for(i in 1:n){
  
#  LL[i]<-f1(b=sim$B[i,],y=y,x=x,alpha=offset2,wt=wt/sim$Sigma[i])	
  LL[i]=-f1(b=sim$B[i,],y=y,x=x,alpha=offset2,wt=wt/sim$Sigma[i])	
}

outlist<-list(coefficients=sim$B,PostMode=sim$BStar,
              Prior=list(mean=as.numeric(mu),Precision=P),
              iters=draws,famfunc=famfunc,Envelope=NULL,
              dispersion=sim$Sigma,loglike=LL)

colnames(outlist$coefficients)<-colnames(x)

outlist$call<-match.call()

class(outlist)<-c(outlist$class,"rglmb")

return(outlist)

}