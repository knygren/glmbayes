data(menarche)

summary(menarche)
plot(Menarche/Total ~ Age, data=menarche)

Age2=menarche$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche$Menarche/menarche$Total
wt=menarche$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-(0.1*(log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(0.1*mu[2,1]/2)^2  # Allows slope to be up to 1 times as large as point estimate 

out<-rglmb(n = 1000, y=y, x=x, mu=mu, P=solve(V1), wt = wt, 
family = binomial(logit), Gridtype = 3) 
summary(out)

famfunc<-glmbfamfunc(binomial(logit))

f1<-famfunc$f1
f2<-famfunc$f2
f3<-famfunc$f3
f5<-famfunc$f5
f6<-famfunc$f6

dispersion2<-as.numeric(1.0)
start <- mu
offset2=rep(as.numeric(0.0),length(y))
P=solve(V1)
n=1000

## Appears that the type for some of these arguments are important/problematic

outlist<-glmbsim_NGauss_cpp(n=as.integer(n),y=as.vector(y),
x=as.matrix(x),mu=as.vector(mu),P=as.matrix(P),
offset2=as.vector(offset2),wt=as.vector(wt),dispersion=as.numeric(dispersion2),
famfunc=famfunc,f1=f1,f2=f2,f3=f3,start=as.vector(start),family="binomial",
link="logit",Gridtype=as.integer(3))



### This allows use of the rglmb summary function 
### add interface for glmbsim_NGauss_cpp later

#outlist$call<-match.call()
colnames(outlist$coefficients)<-colnames(x)
class(outlist)<-c(outlist$class,"rglmb")
#summary(outlist)

## This envelope seems fine

Env1=outlist$Envelope
#Env1

#####

### Adjust weight for dispersion

wt2=wt/dispersion2

################################     #####################

# Shift mean vector to offset so that adjusted model has 0 mean

alpha=x%*%as.vector(mu)+offset2
mu2=0*as.vector(mu)
P2=P
x2=x

############################     ####################################

#####  Optimization step to find posterior mode and associated Precision



parin=start-mu

## Corrected this call - sent wrong mu (mu instead of mu2)
opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x),mu=as.vector(mu2),
      P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
      method="BFGS",hessian=TRUE
      )

#opt_out
bstar=opt_out$par  ## Posterior mode for adjusted model
bstar
bstar+as.vector(mu)  # mode for actual model
A1=opt_out$hessian # Approximate Precision at mode


## Try new function here

Standard_Mod=glmb_Standardize_Model(y=as.vector(y), x=as.matrix(x),P=as.matrix(P),
                       bstar=as.matrix(bstar,ncol=1), A1=as.matrix(A1))



bstar2=Standard_Mod$bstar2  ## This matches
A=Standard_Mod$A
x2=Standard_Mod$x2
mu2=Standard_Mod$mu2
P2=Standard_Mod$P2
L2Inv=Standard_Mod$L2Inv
L3Inv=Standard_Mod$L3Inv

## QC Check - These should match

L2Inv%*%L3Inv%*%as.matrix(bstar2,ncol=1)  
bstar

Env2=glmbenvelope_c(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
    as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),as.vector(wt2),
    family="binomial",link="logit",Gridtype=as.integer(3), n=as.integer(n),sortgrid=TRUE)

## These now seem to match

Env1
Env2

