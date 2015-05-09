## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

n<-1000
mysd<-1

mu<-matrix(0,5)
V0<-((mysd)^2)*diag(5)

glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(),x=TRUE)

glmb.D93<-glmb(n=n,counts ~ outcome + treatment, family = poisson(),mu=mu,Sigma=V0,Gridtype=1)

summary(glm.D93)

summary(glmb.D93)


y<-glm.D93$y
x<-glm.D93$x
b1<-glm.D93$coefficients
wt1<-glm.D93$prior.weights
dispersion<-1
wt2<-wt1/dispersion
alpha1<-rep(0,length(y))
mu<-matrix(0,5)
P_0<-1*diag(5)


# Initialize Simulation
n<-1000

betaout<-matrix(0,nrow=n,ncol=length(y))
alphaout<-matrix(0,nrow=n,ncol=length(mu))
xtemp<-diag(1) 
P<-1

P_0<-0.01*P_0

qc1<-rglmb_rand(n=100,y=y,x=x,mu=mu,P_0=P_0,P=P,wt=wt2,dispersion=dispersion,
                nu=NULL,V=NULL,family=poisson(log),offset2=alpha1,start=mu,Gridtype=3)

summary(qc1)


2*qc1$loglike

glmb.D93$Dbar


logLik(glm.D93)
logLik(glmb.D93)

logLik.glmb



#############################################################################
y=as.vector(18)
x=as.matrix(1)
mu=as.vector(0)
P=as.matrix(1)
wt=as.vector(1)



check<-rglmb(n = 1000, y, x, mu, P, wt = wt, dispersion=NULL,nu=NULL,
      V=NULL,family = poisson(), offset2 = rep(0, 1), start = NULL, Gridtype = 2)

summary(check)

famfunc<-glmbfamfunc(poisson(log))
f1<-famfunc$f1
f2<-famfunc$f2
f3<-famfunc$f3
f5<-famfunc$f5
f6<-famfunc$f6

glmbsim_NGauss_cpp(1,as.vector(t(y)),as.matrix(x),mu,as.matrix(P),as.vector(rep(0,1)),wt,dispersion,
                    famfunc,f1,f2,f3,mu,family="poisson",link="log",Gridtype=2)

glmbsim_NGauss2_cpp(1,as.vector(t(y)),as.matrix(x),mu,as.matrix(P),as.vector(rep(0,1)),wt,dispersion,
                   famfunc,f1,f2,f3,mu,family="poisson",link="log",Gridtype=2)

check<-rglmb(n = 1, y, x, mu, P, wt = wt, dispersion=NULL,nu=NULL,
             V=NULL,family = poisson(), offset2 = rep(0, 1), start = log(18), Gridtype = 2)


summary(check)

check[1]

rglmb


