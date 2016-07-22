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


#interupts(10000000000)


y<-glm.D93$y
x<-glm.D93$x
b1<-glm.D93$coefficients
wt1<-glm.D93$prior.weights
dispersion<-1
wt2<-wt1/dispersion
alpha1<-rep(0,length(y))
mu<-matrix(0,5)

P<-0.1   # Determines beta constant (Smaller values shrinks it  which helps convergence)
P<-1.0*P

m0<-1.5

m0/(1+m0) # Approximate prior weight

P*t(x)%*%x
P_0<-as.matrix(m0*P*t(x)%*%x)

P_0


qc1<-rglmb_rand(n=1000,y=y,x=x,mu=mu,P_0=P_0,P=P,wt=wt2,dispersion=dispersion,
                nu=NULL,V=NULL,family=poisson(log),offset2=alpha1,start=mu,Gridtype=3,
                epsilon_converge=0.01)


summary(qc1)



summary(qc1$randcoefficients)


library(coda)


mcmcout<-mcmc(qc1$coefficients)
plot(mcmcout)



mcmcout2<-mcmc(qc1$randcoefficients)
plot(mcmcout2)



effectiveSize(mcmcout2)
autocorr.plot(mcmcout2)



densplot(mcmcout2)




library(plotly)
library(coda)



tm <- seq(1, 100, by = 10)

plot_ly(x = c(1:100), y = qc1$coefficients[1:100,1], text = paste(tm, "Iteration"))

hist(qc1$coefficients[1:100,1],prob=TRUE)

lines(density(qc1$coefficients[1:100,1]))

hist(qc1$coefficients[1:100,2],prob=TRUE)

lines(density(qc1$coefficients[1:100,2]))

hist(qc1$coefficients[1:100,3],prob=TRUE)

lines(density(qc1$coefficients[1:100,3]))

hist(qc1$coefficients[1:100,4],prob=TRUE)

lines(density(qc1$coefficients[1:100,4]))


plot_ly(x = c(1:100), y = qc1$coefficients[1:100,2], text = paste(tm, "Iteration"))

plot_ly(x = c(1:100), y = qc1$coefficients[1:100,3], text = paste(tm, "Iteration"))


plot_ly(x = c(1:100), y = qc1$coefficients[1:100,4], text = paste(tm, "Iteration"))


plot_ly(x = c(1:3732), y = qc1$coefficients[1:3732,1], text = paste(tm, "Iteration"))

plot_ly(x = c(1:1000), y = qc1$coefficients[1:1000,1], text = paste(tm, "Iteration"))

hist(qc1$coefficients[1:1000,1],prob=TRUE)

lines(density(qc1$coefficients[1:1000,1]))

hist(qc1$coefficients[1:3732,1],prob=TRUE)

lines(density(qc1$coefficients[1:3732,1]))

hist(qc1$coefficients[1:1000,2],prob=TRUE)

lines(density(qc1$coefficients[1:1000,2]))

hist(qc1$coefficients[1:3732,2],prob=TRUE)

lines(density(qc1$coefficients[1:3732,2]))




library(ggplot2)

help(plot)

XTPX<-t(x)%*%(P*diag(9))%*%x

eigdecomp<-eigen(XTPX)

eigdecomp
values<-eigdecomp$values
vectors<-eigdecomp$vectors

values
vectors

values^(-0.5) * t(vectors)


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


