
rm(list = ls())

data(carinsca)
carinsca$Merit <- ordered(carinsca$Merit)
carinsca$Class <- factor(carinsca$Class)
options(contrasts=c("contr.treatment","contr.treatment"))
attach(carinsca)
out <- glm(Claims/Insured~Merit+Class,family="poisson")

summary(out,cor=F)


out <- glm(Cost/Claims~Merit+Class,family=Gamma(link="log"),weights=Claims,x=TRUE)
outsum<-summary(out)

sqrt(diag(outsum$cov.unscaled))
sqrt(diag(outsum$cov.scaled))

y1<-out$y
x1<-out$x
b1<-out$coefficients
wt1<-out$prior.weights

alpha1<-rep(0,length(y1))


mu<-matrix(0,8)
P<-10*diag(8)


scale<-0.1275 # SAS estimate
dispersion<-1/scale



glm(Cost/Claims~Merit+Class,family=Gamma(link="log"),weights=Claims)


start.time<-Sys.time()

outb <- glmb(n=1000,Cost/Claims~Merit+Class,family=Gamma(link="log"),
             mu=mu,Sigma=0.1*diag(8),dispersion=dispersion, weights=Claims,x=TRUE,Gridtype=2)
end.time<-Sys.time()
time.taken<-end.time-start.time
print("Time for full simulation")
print(time.taken)


summary(outb)

# Get family functions

f1<-outb$famfunc$f1
f2<-outb$famfunc$f2

f1(b=b1,y=y1,x=x1,alpha=alpha1,wt=wt1)

b2<-matrix(b1,nrow=8,ncol=1)
f1_gamma(b=b2,y=y1,x=x1,alpha=alpha1,wt=wt1)


l2<-3^8
NegLL1<-matrix(0,nrow=l2,ncol=1)
b2<-matrix(b1,nrow=8,ncol=3^8)
b2[,1:10]

start.time<-Sys.time()

for(j in 1:l2){
NegLL1[j]<-f2(b=b2[,j],y=y1,x=x1,mu=mu,P=P,alpha=alpha1,wt=wt1)
}

end.time<-Sys.time()
end.time-start.time




start.time<-Sys.time()

NegLL2<-f2_gamma(b=b2,y=y1,x=x1,mu=mu,P=P,alpha=alpha1,wt=wt1)
end.time<-Sys.time()
end.time-start.time

NegLL1[1:10]

NegLL2[1:10]



as.matrix(y1)
as.matrix(alpha1)
as.matrix(wt1)




objects()


out2<-rglmb(n = 1000, out$y, out$x, mu=mu, P=P, wt = out$prior.weights, dispersion = dispersion, family = Gamma(link="log"), 
            offset2 = rep(0, 20), start = out$coefficients, Gridtype = 2) 

summary(out2)

a0<-0.01
m0<-0.01
b<-colMeans(out2$coefficients)

out3<-rglmbdisp(n=10000,y=out$y,x=out$x,b=b,alpha= rep(0, length(out$y)),wt=out$prior.weights,shape=m0,rate=a0,family=Gamma(link=log))

summary(out3)

f1<-out2$famfunc$f1


f4<-out2$famfunc$f4


# OpenBugs appears to use precision to calculate DIC

wt2<-out$prior.weights*mean(1/out3)


Dthetabar<-2*f1(b,y=out$y,x=out$x,alpha=rep(0,20),wt=wt2)

Dthetabar

Devtemp<-matrix(0,10)

for(j in 1:10000)
{
  wttemp<-out$prior.weights/out3[j]
  Devtemp[j]<-2*f1(b,y=out$y,x=out$x,alpha=rep(0,20),wt=wttemp)
}

Devtemp[1:10]

Dbar<-mean(Devtemp)

pD<-Dbar-Dthetabar
pD
DIC<-pD+Dbar
DIC

######  Joint Estimation - Two -Block Gibbs Sampler ###


y1<-out$y
x1<-out$x
b1<-out$coefficients
wt1<-out$prior.weights

alpha1<-rep(0,length(y1))


mu<-matrix(0,8)
P<-10*diag(8)


scale<-0.1275 # SAS estimate
dispersion<-1/scale

a0<-0.01
m0<-0.01

outbetas<-matrix(0,nrow=11000,ncol=8)
outdisp<-matrix(0,nrow=11000)

disp<-dispersion
start.time<-Sys.time()
for(i in 1:100)
{
  start.time1<-Sys.time()
  
  outtemp1<-rglmb(n = 1, out$y, out$x, mu=mu, P=P, wt = out$prior.weights, dispersion = disp, family = Gamma(link="log"), 
        offset2 = rep(0, 20), start = out$coefficients, Gridtype = 2) 
  end.time1<-Sys.time()
  time.taken1<-end.time1-start.time1
#  print("Time for rglmb")
#  print(time.taken1)
  
  outbetas[i,1:8]<-outtemp1$coefficients[1,1:8]
  b<-outbetas[i,1:8]  
  start.time2<-Sys.time()  
  outdisp[i,1]<-rglmbdisp(n=1,y=out$y,x=out$x,b=b,alpha= rep(0, length(out$y)),wt=out$prior.weights,shape=m0,rate=a0,family=Gamma(link=log))
  end.time2<-Sys.time()
  time.taken2<-end.time2-start.time2
#  print("Time for rglmbdisp")
#  print(time.taken2)
  
#  print(i) 
  disp<-outdisp[i,1]  
}  
end.time<-Sys.time()
time.taken<-end.time-start.time
print("Time for loop")
print(time.taken)



mean(outdisp[40:100,1])

colMeans(outbetas[40:100,])
