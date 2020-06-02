## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group,x=TRUE,y=TRUE)

lm_summary=summary(lm.D9)

lm_summary
lm.D9$coefficients
lm_summary$sigma


n<-1000
y<-lm.D9$y
x<-as.matrix(lm.D9$x)


mu<-c(0,0)
P<-0.1*diag(2)

#solve(P)
x
b=lm.D9$coefficients
b_old=lm.D9$coefficients

v0=lm_summary$sigma
Prior_Sd=(0.1*sqrt(v0))^2

v0/Prior_SD=shape/sqrt(shape)=sqrt(shape)
shape=(v0/Prior_Sd)^2
rate=shape/v0
shape
rate

dispout<-rglmbdisp(n=n,y,x,b,alpha= rep(0, length(y)),shape=shape,rate=rate,family=gaussian())
summary(dispout)

dispout$coefficients
class(dispout$coefficients)


n
y
x
mu
P
dispout
b


### This part works now

wt=1
wt2=0*c(1:length(y))+wt


outtemp<-rglmb(n = 1000, y, x, mu, P, wt = wt2, dispersion=dispout$coefficient[1],family = gaussian(), offset2 = rep(0, length(y)), start = b, Gridtype = 3)
outtemp<-rglmb(n = 1000, y, x, mu, P, wt = 1, dispersion=dispout$coefficient[1],family = gaussian(), offset2 = rep(0, length(y)), start = b, Gridtype = 3)

summary(outtemp)


help(rglmb)


# RUn Two-Block Gibbs Sampler
n=100

dispout=0*c(1:n)
bout=matrix(0,nrow=n,ncol=length(mu))

## Edited start




help(rglmb)
help(rglmbdisp)

## Seems to be triggered only some way into simulation - need to backtrack and then add
## error check to catch this earlier

# error:  chol(): decomposition failed

rep(0,n)
diag(diag(length(y)))

n=1000

rlm1<-function(n,y,x,mu,P,shape,rate,alpha,start,burnin=100){

dispout=rep(0,n)
bout=matrix(0,nrow=n,ncol=length(mu))

for(i in 1:burnin){
  
  out=rglmbdisp(1,y,x,start,alpha= rep(0, length(y)),shape=shape,rate=rate,family=gaussian())
  dispout[i]<-out$coefficient
  outtemp<-rglmb(n = 1, y, x, mu, P, wt = 1, dispersion=dispout[i],family = gaussian(), offset2 = rep(0, length(y)), start = start, Gridtype = 3)
  
  bout[i,]<-outtemp$coefficients
  b<-t(as.matrix(outtemp$coefficients))
  
}

  
for(i in 1:n){

out=rglmbdisp(1,y,x,start,alpha= rep(0, length(y)),shape=shape,rate=rate,family=gaussian())
dispout[i]<-out$coefficient
outtemp<-rglmb(n = 1, y, x, mu, P, wt = 1, dispersion=dispout[i],family = gaussian(), offset2 = rep(0, length(y)), start = start, Gridtype = 3)

bout[i,]<-outtemp$coefficients
b<-t(as.matrix(outtemp$coefficients))

      }

outlist1=list(coefficients=dispout,Prior=list(shape=shape,rate=rate))
outlist1$call<-match.call()

class(outlist1)<-c(outlist1$class,"rglmbdisp")

return(list(bout=bout,dispout=outlist1))
}


temp=rlm1(n,y,x,mu,P,shape,rate,alpha,start=b_old,burnin=100)  
summary(temp$dispout)
summary(temp$bout)








outlist=list(coefficients=dispout,Prior=list(shape=shape,rate=rate))
outlist$call=


summary(bout)
summary(dispout)



i=1
out=rglmbdisp(1,y,x,b,alpha= rep(0, length(y)),shape=shape,rate=rate,family=gaussian())
dispout[i]<-out$coefficients
outtemp<-rglmb(n = 1, y, x, mu, P, wt = 1, dispersion=dispout[i],family = gaussian(), offset2 = rep(0, length(y)), start = b_old, Gridtype = 3)
bout[i,]<-outtemp$coefficients
b<-t(as.matrix(outtemp$coefficients))

i=2
out=rglmbdisp(1,y,x,b,alpha= rep(0, length(y)),shape=shape,rate=rate,family=gaussian())
dispout[i]<-out$coefficients
outtemp<-rglmb(n = 1, y, x, mu, P, wt = 1, dispersion=dispout[i],family = gaussian(), offset2 = rep(0, length(y)), start = b_old, Gridtype = 3)
bout[i,]<-outtemp$coefficients
b<-t(as.matrix(outtemp$coefficients))


i=3
out=rglmbdisp(1,y,x,b,alpha= rep(0, length(y)),shape=shape,rate=rate,family=gaussian())
dispout[i]<-out$coefficients
outtemp<-rglmb(n = 1, y, x, mu, P, wt = 1, dispersion=dispout[i],family = gaussian(), offset2 = rep(0, length(y)), start = b_old, Gridtype = 3)
bout[i,]<-outtemp$coefficients
b<-t(as.matrix(outtemp$coefficients))


dispout

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


