## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)

ps=Prior_Setup(weight ~ group)
mu=ps$mu
V=ps$Sigma
mu[1,1]=mean(weight)

Prior_Check(weight ~ group,family =gaussian(),
           pfamily=dNormal(mu=mu,Sigma=V))

## Will move this step inside the Prior_Check
lm.D9 <- lm(weight ~ group,x=TRUE,y=TRUE)
disp_ML=sigma(lm.D9)^2
n_prior=2
shape=n_prior/2
rate= disp_ML*shape

lmb.D9=lmb(weight ~ group,dNormal_Gamma(mu,V/disp_ML,shape=shape,rate=rate))
summary(lmb.D9)


n<-10000
y<-lmb.D9$y
x<-as.matrix(lmb.D9$x)

prior_list=list(mu=mu,Sigma=V/disp_ML,shape=shape,rate=rate)

ng.D9=rnorm_gamma_reg(n=1000,y=y,x=x,
  prior_list=prior_list)

print(ng.D9)
#summary(ng.D9)  This failed - Need to investigate
summary(ng.D9$coefficients)

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
prior_list=list(shape=shape,rate=rate,beta=b_old)
dispout<-rglmb_dispersion(n=n,y,x,prior_list=prior_list,
offset= rep(0, length(y)),family=gaussian())

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
prior=list(mu=mu,Sigma=solve(P),dispersion=v_old)
outtemp1<-glmb(n = 1000, weight ~ group, family = gaussian(),
pfamily=dNormal(mu=mu,Sigma=solve(P),dispersion=v_old))
## Could use a residuals function here -- For now, maybe run the glmb function

summary(outtemp1)
mean(colMeans(residuals(outtemp1)^2))
v_old
colMeans((outtemp1$coefficients))
b_old
prior=list(mu=mu,P=P,dispersion=v_old)
outtemp2<-rglmb(n = 1000, y, x, family = gaussian(),
pfamily=dNormal(mu=mu,Sigma=solve(P),dispersion=v_old),
offset = rep(0, length(y)), weights = wt2)
summary(outtemp2)
colMeans((outtemp2$coefficients))
b_old

prior=list(mu=mu,Sigma=solve(P),shape=shape,rate=rate)
outtemp3<-glmb(n = 10000, weight ~ group,family = gaussian(),  
dNormal_Gamma(mu=mu,Sigma=solve(P),shape=shape,rate=rate))

summary(outtemp3)
mean(colMeans(residuals(outtemp3)^2))
v_old
colMeans((outtemp3$coefficients))
b_old
mean(outtemp3$dispersion)  ## Seems slightly smaller --> Needs qc
v_old
