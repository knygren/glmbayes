## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group,x=TRUE,y=TRUE)

n<-1000
y<-lm.D9$y
x<-lm.D9$x
mu<-c(0,0)
P<-0.001*diag(2)
x
shape<-0.001
rate<-0.001

n<-11000
dispout<-matrix(0,n)
bout<-matrix(0,nrow=n,ncol=length(mu))

#initialize b

b<-as.matrix(mu)

# RUn Two-Block Gibbs Sampler



for(i in 1:11000){

dispout[i]<-rglmbdisp(1,y,x,b,alpha= rep(0, length(y)),shape,rate,family=gaussian())

outtemp<-rglmb(n = 1, y, x, mu, P, wt = 1, dispersion=dispout[i],family = gaussian(), offset2 = rep(0, length(y)), start = NULL, Gridtype = 2)

bout[i,]<-outtemp$coefficients
b<-t(as.matrix(outtemp$coefficients))

}

alpha<-rep(0, length(y))
mu1<-exp(alpha+x%*%b)


SS<-(y-mu1)/mu1-log(y/mu1)
SS

library(ars)

rglmbdisp(1,y,x,b,alpha= rep(0, length(y)),shape,rate,family=Gamma(log))



# MEan estimate for regression coefficients

colMeans(bout[1001:11000,])

# Standard deviation 

sqrt(diag(var(bout[1001:11000,])))

#  Comparison to lm Residual standard errors

sqrt(1/mean(1/dispout))

bout

sqrt(mean(disp))

summary(lm.D9)


library("R2OpenBUGS")

x<-lm.D9$x
y<-lm.D9$y
mu<-c(0,0)
P<-0.001*diag(2)


n<-length(y)
l1<-length(mu)
offset2<-0*c(1:n)

b<-matrix(1,l1)

alpha<-0
n1<-length(y)
y1<-as.matrix(y)-alpha
xb<-x%*%b
res<-y1-xb
SS<-res*res







data<-list(n=n,l1=l1,y=y,X=x,offset=offset2,mu=mu,P=P)

inits<-function(){list(beta=c(0,0),Pdata=1)}

parameters<-c("beta","Pdata","Sigma")
model.file <- system.file("model", "gaussian.txt", package="R2OpenBUGS")
height.sim <- bugs(data, inits, parameters, model.file,n.chains = 3, n.iter = 1000,working.directory = NULL)

attach.bugs(height.sim)

height.sim
colMeans(beta)

sqrt(diag(var(beta)))

summary(lm.D9)
#sqrt(mean(Sigma))
#mean(sqrt(Sigma))

sqrt(1/mean(Pdata))


