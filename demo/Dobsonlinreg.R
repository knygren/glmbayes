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
### This part works now

wt=1
wt2=0*c(1:length(y))+wt

outtemp<-rglmb(n = 1000, y, x, mu, P, wt = wt2, dispersion=dispout$coefficient[1],family = gaussian(), offset2 = rep(0, length(y)), start = b, Gridtype = 3)

outtemp<-rglmb(n = 1000, y, x, mu, P, wt = wt2,dispersion=dispout$coefficient[1],family = gaussian(), offset2 = rep(0, length(y)), start = b, Gridtype = 3)

summary(outtemp)


outtemp2<-rglmb(n = 1000, y, x, mu, P, wt = wt2, dispersion=NULL,
nu=3,V=3,family = gaussian(), offset2 = rep(0, length(y)), start = b, Gridtype = 3)
summary(outtemp2)

#outtemp2

outtemp3<-rnorm_gamma_reg(n=1000,y,x,mu,P,nu=3,V=3,offset2= rep(0, length(y)),wt=wt2)

#outtemp3
summary(outtemp3)



