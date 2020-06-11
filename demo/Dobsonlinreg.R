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

n<-10000
y<-lm.D9$y
x<-as.matrix(lm.D9$x)

### Set old regression and dispersion coefficients

b_old=lm.D9$coefficients
v_old=lm_summary$sigma^2

#### Set up for rglmb_dispersion

n0=0.1
shape=n0/2
rate=shape*v_old

rate/(shape-1)
rate/shape

#a1<-shape+n1/2
#b1<-rate+sum(SS)/2

#out<-1/rgamma(n,shape=a1,rate=b1)
## v0=sum(SS)/n1 ~ (2*rate+sum(SS))/(2*shape+n1)=[(2*rate+sum(SS))/2]/((2*shape+n1)/2) 
n1=length(y)
SS=v_old*n1

n1 # This is equal to 20


n0=2 # Prior observations
v_prior=v_old  # Prior point estimate for variance (the mean of (1/dispersion=1/v_prior))
wt0=(n0/n1)  

## set shape=0.01*(n1/2)
## set rate= 0.01*SS/2

shape=wt0*(n1/2)   ###  Shape is prior observations /2
rate=shape*v_prior  ### rate is essentiall prior SS - V in rmultireg should be this
rate/shape ## Should match v_prior (currently also v_old)

## We see that this currently matches 
### (test different v_prior with various prior observations below) 

dispout<-rglmb_dispersion(n=n,y,x,b_old,alpha= rep(0, length(y)),
shape=shape,rate=rate,family=gaussian())
mean(dispout$dispersion) 
v_prior
v_old
#summary(dispout)  ## Summary function not working - need to write

### Set up test rglmb regressions without and with prior for dispersion 

mu<-c(0,0)
mu=b_old  ### For testing purposes, set prior=b_old
P<-0.1*diag(2)
wt2=rep(1,length(y))

### Check

outtemp1<-glmb(n = 1000, weight ~ group, mu=mu, Sigma = solve(P), 
               dispersion=v_old,family = gaussian(),  start =b_old, Gridtype = 3)
## Could use a residuals function here -- For now, maybe run the glmb function

summary(outtemp1)
mean(colMeans(residuals(outtemp1)^2))
v_old
colMeans((outtemp1$coefficients))
b_old

outtemp2<-rglmb(n = 1000, y, x, mu=mu, P=P, wt = wt2, dispersion=v_old,
family = gaussian(), offset2 = rep(0, length(y)), start = b_old, Gridtype = 3)
summary(outtemp2)
colMeans((outtemp2$coefficients))
b_old


outtemp3<-rglmb(n = 1000, y, x, mu=mu, P=P, wt = wt2, shape=shape, V=rate,
                family = gaussian(), offset2 = rep(0, length(y)), start = b_old, Gridtype = 3)
summary(outtemp2)
colMeans((outtemp2$coefficients))
b_old


outtemp3<-glmb(n = 10000, weight ~ group, mu=mu, Sigma = solve(P), shape=shape, V=rate,
family = gaussian(),  start = b_old, Gridtype = 3)
## Could use a residuals function here -- For now, maybe run the glmb function

summary(outtemp3)
mean(colMeans(residuals(outtemp3)^2))
v_old
colMeans((outtemp3$coefficients))
b_old
mean(outtemp3$dispersion)  ## Seems slightly smaller --> Needs qc
v_old
rm(dispersion)
