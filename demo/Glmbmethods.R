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

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 

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
famfunc=famfunc,f1=f1,f2=f2,f3=f3,start=as.vector(start),family="binomial",link="logit",
Gridtype=as.integer(3))

#####

### Adjust weight for dispersion

wt2=wt/dispersion2

# Shift mean vector to offset so that adjusted model has 0 mean

alpha=x%*%as.vector(mu)+offset2
mu2=0*as.vector(mu)

#####  Optimization step to find posterior mode and associated Precision

parin=start-mu

opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x),mu=as.vector(mu),
      P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
      method="BFGS",hessian=TRUE
      )

#opt_out
b2=opt_out$par  ## Posterior mode for adjusted model
b2+as.vector(mu)  # mode for actual model
A1=opt_out$hessian # Approximate Precision at mode

# Find eigenvalues and standardize to model with ~ Identity precision at posterior mode

A_eigen=eigen(A1)
D1=diag(A_eigen$values)
L2=sqrt(D1)%*%t(A_eigen$vectors)
L2Inv=A_eigen$vectors%*%sqrt(solve(D1))

## Apply transformation

#arma::mat b3=L2*b2;   
#arma::mat mu3=L2*mu2; // These are needed but will not be used to pass 

#arma::mat x3=x2*L2Inv;
#arma::mat P3=trans(L2Inv)*P2*L2Inv;






