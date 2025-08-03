pkgname <- "glmbayes"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('glmbayes')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("AMI")
### * AMI

flush(stderr()); flush(stdout())

### Name: AMI
### Title: amitriptyline overdose data
### Aliases: AMI
### Keywords: datasets, Bayesian Binomial Regression

### ** Examples

data(menarche2)



cleanEx()
nameEx("EnvelopeBuild")
### * EnvelopeBuild

flush(stderr()); flush(stdout())

### Name: EnvelopeBuild
### Title: Builds Envelope function for simulation
### Aliases: EnvelopeBuild
### Keywords: internal

### ** Examples

data(menarche2)

summary(menarche2)
plot(Menarche/Total ~ Age, data=menarche2)

Age2=menarche2$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche2$Menarche/menarche2$Total
wt=menarche2$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 1 times as large as point estimate 

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


###### Adjust weight for dispersion

wt2=wt/dispersion2

######################### Shift mean vector to offset so that adjusted model has 0 mean

alpha=x%*%as.vector(mu)+offset2
mu2=0*as.vector(mu)
P2=P
x2=x


#####  Optimization step to find posterior mode and associated Precision

parin=start-mu

opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x),mu=as.vector(mu2),
              P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
              method="BFGS",hessian=TRUE
)

bstar=opt_out$par  ## Posterior mode for adjusted model
bstar
bstar+as.vector(mu)  # mode for actual model
A1=opt_out$hessian # Approximate Precision at mode

## Standardize Model

Standard_Mod=glmb_Standardize_Model(y=as.vector(y), x=as.matrix(x),P=as.matrix(P),
                                    bstar=as.matrix(bstar,ncol=1), A1=as.matrix(A1))

bstar2=Standard_Mod$bstar2  
A=Standard_Mod$A
x2=Standard_Mod$x2
mu2=Standard_Mod$mu2
P2=Standard_Mod$P2
L2Inv=Standard_Mod$L2Inv
L3Inv=Standard_Mod$L3Inv

Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),as.vector(wt2),
family="binomial",link="logit",Gridtype=as.integer(3), n=as.integer(n),
sortgrid=TRUE)

## These now seem to match

Env2



cleanEx()
nameEx("EnvelopeOpt")
### * EnvelopeOpt

flush(stderr()); flush(stdout())

### Name: EnvelopeOpt
### Title: Optimizes Envelope function for simulation
### Aliases: EnvelopeOpt
### Keywords: internal

### ** Examples

data(menarche2)

summary(menarche2)
plot(Menarche/Total ~ Age, data=menarche2)

Age2=menarche2$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche2$Menarche/menarche2$Total
wt=menarche2$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 1 times as large as point estimate 

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



###### Adjust weight for dispersion

wt2=wt/dispersion2

######################### Shift mean vector to offset so that adjusted model has 0 mean

alpha=x%*%as.vector(mu)+offset2
mu2=0*as.vector(mu)
P2=P
x2=x


#####  Optimization step to find posterior mode and associated Precision

parin=start-mu

opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x),mu=as.vector(mu2),
              P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
              method="BFGS",hessian=TRUE
)

bstar=opt_out$par  ## Posterior mode for adjusted model
bstar
bstar+as.vector(mu)  # mode for actual model
A1=opt_out$hessian # Approximate Precision at mode

## Standardize Model

Standard_Mod=glmb_Standardize_Model(y=as.vector(y), x=as.matrix(x),
P=as.matrix(P),bstar=as.matrix(bstar,ncol=1), A1=as.matrix(A1))

bstar2=Standard_Mod$bstar2  
A=Standard_Mod$A
x2=Standard_Mod$x2
mu2=Standard_Mod$mu2
P2=Standard_Mod$P2
L2Inv=Standard_Mod$L2Inv
L3Inv=Standard_Mod$L3Inv

Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),as.vector(wt2),
family="binomial",link="logit",Gridtype=as.integer(3), 
n=as.integer(n),sortgrid=TRUE)

## These now seem to match

Env2



cleanEx()
nameEx("EnvelopeSort")
### * EnvelopeSort

flush(stderr()); flush(stdout())

### Name: EnvelopeSort
### Title: Sorts Envelope function for simulation
### Aliases: EnvelopeSort
### Keywords: internal

### ** Examples

data(menarche2)

summary(menarche2)
plot(Menarche/Total ~ Age, data=menarche2)

Age2=menarche2$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche2$Menarche/menarche2$Total
wt=menarche2$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 1 times as large as point estimate 

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

#outlist<-rnnorm_reg_cpp(n=as.integer(n),y=as.vector(y),x=as.matrix(x),
#mu=as.vector(mu),P=as.matrix(P),offset2=as.vector(offset2),
#wt=as.vector(wt),dispersion=as.numeric(dispersion2),famfunc=famfunc,
#f1=f1,f2=f2,f3=f3,start=as.vector(start),family="binomial",
#link="logit",Gridtype=as.integer(3))

### This allows use of the rglmb summary function 
### add interface for glmbsim_NGauss_cpp later

#outlist$call<-match.call()
#colnames(outlist$coefficients)<-colnames(x)
#class(outlist)<-c(outlist$class,"rglmb")
#summary(outlist)
#Env1=outlist$Envelope


###### Adjust weight for dispersion

wt2=wt/dispersion2

######################### Shift mean vector to offset so that adjusted model has 0 mean

alpha=x%*%as.vector(mu)+offset2
mu2=0*as.vector(mu)
P2=P
x2=x


#####  Optimization step to find posterior mode and associated Precision

parin=start-mu

opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x),mu=as.vector(mu2),
              P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
              method="BFGS",hessian=TRUE
)

bstar=opt_out$par  ## Posterior mode for adjusted model
bstar
bstar+as.vector(mu)  # mode for actual model
A1=opt_out$hessian # Approximate Precision at mode

## Standardize Model

Standard_Mod=glmb_Standardize_Model(y=as.vector(y), x=as.matrix(x),
P=as.matrix(P),bstar=as.matrix(bstar,ncol=1), A1=as.matrix(A1))

bstar2=Standard_Mod$bstar2  
A=Standard_Mod$A
x2=Standard_Mod$x2
mu2=Standard_Mod$mu2
P2=Standard_Mod$P2
L2Inv=Standard_Mod$L2Inv
L3Inv=Standard_Mod$L3Inv

Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),as.vector(wt2),
family="binomial",link="logit",Gridtype=as.integer(3), 
n=as.integer(n),sortgrid=TRUE)

## These now seem to match

#Env1
Env2





cleanEx()
nameEx("Inv_f3_gaussian")
### * Inv_f3_gaussian

flush(stderr()); flush(stdout())

### Name: Inv_f3_gaussian
### Title: Derives the inverse of the gradient vector.
### Aliases: Inv_f3_gaussian

### ** Examples

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb extractAIC-------------------------------------------------------------
extractAIC(glmb.D93)



cleanEx()
nameEx("Neg_logLik")
### * Neg_logLik

flush(stderr()); flush(stdout())

### Name: Neg_logLik
### Title: Negative Log-Likelihood for a Generalized Linear Model
### Aliases: Neg_logLik Neg_logLik2
### Keywords: internal

### ** Examples

1+1
10+1



cleanEx()
nameEx("Normal_ct")
### * Normal_ct

flush(stderr()); flush(stdout())

### Name: Normal_ct
### Title: The Central Normal Distribution
### Aliases: Normal_ct pnorm_ct rnorm_ct
### Keywords: internal

### ** Examples

pnorm_ct(0.2,0.4)
exp(pnorm_ct(0.2,0.4))
pnorm_ct(0.2,0.4,log.p=FALSE)
log(pnorm_ct(0.2,0.4,log.p=FALSE))
## Example where difference between two pnorm calls fail
## but call to pnorm_ct works
pnorm(0.5)-pnorm(0.4999999999999999)
pnorm_ct(0.4999999999999999,0.5,log.p=FALSE)



cleanEx()
nameEx("Prior_Check")
### * Prior_Check

flush(stderr()); flush(stdout())

### Name: Prior_Check
### Title: Checks for Prior-data conflicts
### Aliases: Prior_Check

### ** Examples

## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
## Step 1: Set up Prior
ps=Prior_Setup(counts ~ outcome + treatment)
mu=ps$mu
V=ps$Sigma
# Step2A: Check the Prior
Prior_Check(counts ~ outcome + treatment,family = poisson(),
            pfamily=dNormal(mu=mu,Sigma=V))
# Step2B: Update and Re-Check the Prior
mu[1,1]=log(mean(counts))
Prior_Check(counts ~ outcome + treatment,family = poisson(),
            pfamily=dNormal(mu=mu,Sigma=V))



cleanEx()
nameEx("Prior_Setup")
### * Prior_Setup

flush(stderr()); flush(stdout())

### Name: Prior_Setup
### Title: Setup Prior Objects
### Aliases: Prior_Setup

### ** Examples

## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
## Step 1: Set up Prior
ps=Prior_Setup(counts ~ outcome + treatment)
mu=ps$mu
V=ps$Sigma
print(ps)



cleanEx()
nameEx("Set_Grid")
### * Set_Grid

flush(stderr()); flush(stdout())

### Name: Set_Grid
### Title: Calculate Log-densities for Grid Components
### Aliases: Set_Grid
### Keywords: internal

### ** Examples

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb extractAIC-------------------------------------------------------------
extractAIC(glmb.D93)



cleanEx()
nameEx("anova.glmb")
### * anova.glmb

flush(stderr()); flush(stdout())

### Name: anova.glmb
### Title: Analysis of Deviance for Bayesian Generalized Linear Model Fits
### Aliases: anova.glmb

### ** Examples

data(menarche2)
## ----Analysis Setup-----------------------------------------------------------
## Number of variables in model
Age=menarche2$Age
nvars=2
## Reference Ages for setting of priors and Age_Difference
ref_age1=13  # user can modify this
ref_age2=15  ## user can modify this
## Define variables used later in analysis
Age2=menarche2$Age-ref_age1
Age_Diff=ref_age2-ref_age1
mu1=as.matrix(c(0,1.098612),ncol=1)
V1<-1*diag(nvars)
V1[1,1]=0.18687882
V1[2,2]=0.10576217
V1[1,2]=-0.03389182
V1[2,1]=-0.03389182
Menarche_Model_Data=data.frame(Age=menarche2$Age,Total=menarche2$Total,
                               Menarche=menarche2$Menarche,Age2)
glmb.out1<-glmb(n=1000,cbind(Menarche, Total-Menarche) ~Age2,family=binomial(logit),
                pfamily=dNormal(mu=mu1,Sigma=V1),data=Menarche_Model_Data)

# Prediction from original model
pred1=predict(glmb.out1,type="response")

## Get Original Residuals, their means, and credible bounds
res_out=residuals(glmb.out1)
colMeans(res_out)

## Set up to simulate new data and residuals
res_mean=colMeans(res_out)
res_low1=apply(res_out,2,FUN=quantile,probs=c(0.025))
res_high1=apply(res_out,2,FUN=quantile,probs=c(0.975))

## Simulate new data and get residuals for simulated data

ysim1=simulate(glmb.out1,nsim=1,seed=NULL,pred=pred1,family="binomial",
               prior.weights=weights(glmb.out1))


res_ysim_out1=residuals(glmb.out1,ysim=ysim1)
res_low=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.025))
res_high=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.975))

# Plot Credible Interval bounds for Deviance Residuals

plot(res_mean~Age,ylim=c(-2.5,2.5),
main="Credible Interval Bound for Menarche - Logit Model Deviance Residuals",
xlab = "Age", ylab = "Avg. Dev. Res")
lines(Age, 0*res_mean,lty=1)
lines(Age, res_low,lty=1)
lines(Age, res_high,lty=1)
lines(Age, res_low1,lty=2)
lines(Age, res_high1,lty=2)




cleanEx()
nameEx("carinsca")
### * carinsca

flush(stderr()); flush(stdout())

### Name: carinsca
### Title: Canadian Automobile Insurance Claims for 1957-1958
### Aliases: carinsca
### Keywords: datasets, Bayesian Poisson Regression, Bayesian Gamma
###   Regression

### ** Examples

data(carinsca)
## maybe str(carinsca) ; plot(carinsca) ...



cleanEx()
nameEx("case.names.glmb")
### * case.names.glmb

flush(stderr()); flush(stdout())

### Name: case.names.glmb
### Title: Case and Variable Names of Fitted Models
### Aliases: case.names.glmb variable.names.glmb

### ** Examples

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))

## ----glmb vcov----------------------------------------------------------------
case.names(glmb.D93)




cleanEx()
nameEx("confint.glmb")
### * confint.glmb

flush(stderr()); flush(stdout())

### Name: confint.glmb
### Title: Credible Intervals for Model Parameters
### Aliases: confint.glmb

### ** Examples


## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb confint-------------------------------------------------------------
confint(glmb.D93)



cleanEx()
nameEx("deviance.rglmb")
### * deviance.rglmb

flush(stderr()); flush(stdout())

### Name: deviance.rglmb
### Title: Model Deviance
### Aliases: deviance.rglmb

### ** Examples

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb extractAIC-------------------------------------------------------------
extractAIC(glmb.D93)



cleanEx()
nameEx("dummy.coef.glmb")
### * dummy.coef.glmb

flush(stderr()); flush(stdout())

### Name: dummy.coef.glmb
### Title: Extract Coefficients in Original Coding
### Aliases: dummy.coef.glmb print.dummy_coef.glmb

### ** Examples

set.seed(333)
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
mysd<-1
mu<-matrix(0,5)
mu[1]=log(mean(counts))
V0<-((mysd)^2)*diag(5)
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment, 
family = poisson(),pfamily=dNormal(mu=mu,Sigma=V0))
summary(glmb.D93)


myd=dummy.coef(glmb.D93)
print(myd)



cleanEx()
nameEx("extractAIC.glmb")
### * extractAIC.glmb

flush(stderr()); flush(stdout())

### Name: extractAIC.glmb
### Title: Extract DIC from a Fitted Bayesian Model
### Aliases: extractAIC.glmb extractDIC

### ** Examples

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb extractAIC-------------------------------------------------------------
extractAIC(glmb.D93)



cleanEx()
nameEx("extractAIC.rglmb")
### * extractAIC.rglmb

flush(stderr()); flush(stdout())

### Name: extractAIC.rglmb
### Title: Extract DIC from a Fitted Bayesian Model
### Aliases: extractAIC.rglmb

### ** Examples

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb extractAIC-------------------------------------------------------------
extractAIC(glmb.D93)



cleanEx()
nameEx("f2_gaussian_vector")
### * f2_gaussian_vector

flush(stderr()); flush(stdout())

### Name: f2_gaussian_vector
### Title: Derives the inverse of the gradient vector.
### Aliases: f2_gaussian_vector

### ** Examples

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb extractAIC-------------------------------------------------------------
extractAIC(glmb.D93)



cleanEx()
nameEx("formula.summary.rglmb")
### * formula.summary.rglmb

flush(stderr()); flush(stdout())

### Name: formula.summary.rglmb
### Title: Extract Log-Likelihood
### Aliases: formula.summary.rglmb

### ** Examples

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb logLik-------------------------------------------------------------
colMeans(logLik(glmb.D93))



cleanEx()
nameEx("glmb")
### * glmb

flush(stderr()); flush(stdout())

### Name: glmb
### Title: Fitting Bayesian Generalized Linear Models
### Aliases: glmb print.glmb

### ** Examples

set.seed(333)
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
mysd<-1
mu<-matrix(0,5)
mu[1]=log(mean(counts))
V0<-((mysd)^2)*diag(5)
glmb.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu,Sigma=V0))
summary(glmb.D93)

# Menarche Binomial Data Example 
data(menarche2)
Age2=menarche2$Age-13

#Priors Are Derived in Vignette
nvars=2
## Logit Prior and Model
mu1=as.matrix(c(0,1.098612),ncol=1)
V1<-1*diag(nvars)
V1[1,1]=0.18687882
V1[2,2]=0.10576217
V1[1,2]=-0.03389182
V1[2,1]=-0.03389182
glmb.out1<-glmb(cbind(Menarche, Total-Menarche) ~ Age2,
family=binomial(logit),pfamily=dNormal(mu=mu1,Sigma=V1), data=menarche2)
summary(glmb.out1)

## Probit Prior and Model
mu2=as.matrix(c(0,0.6407758),ncol=1)
V2<-1*diag(nvars)
V2[1,1]=0.07158369
V2[2,2]=0.03453205
V2[1,2]=-0.01512075
V2[2,1]=-0.01512075
glmb.out2<-glmb(cbind(Menarche, Total-Menarche) ~ Age2,
family=binomial(probit),pfamily=dNormal(mu=mu2,Sigma=V2), data=menarche2)
summary(glmb.out2)

## clog-log Prior and Model
mu2=as.matrix(c(-0.3665129,0.6002727),ncol=1)
V2<-1*diag(nvars)
V2[1,1]=0.11491322
V2[2,2]=0.03365986
V2[1,2]=-0.03502783
V2[2,1]=-0.03502783
glmb.out3<-glmb(cbind(Menarche, Total-Menarche) ~ Age2,
                family=binomial(cloglog),pfamily=dNormal(mu=mu2,Sigma=V2), data=menarche2)
summary(glmb.out3)

DIC_Out=rbind(extractAIC(glmb.out1),extractAIC(glmb.out2),extractAIC(glmb.out3))
rownames(DIC_Out)=c("logit","probit","clog-log")
colnames(DIC_Out)=c("pD","DIC")
DIC_Out



cleanEx()
nameEx("glmb.influence.measures")
### * glmb.influence.measures

flush(stderr()); flush(stdout())

### Name: glmb.influence.measures
### Title: Bayesian Regression Diagnostics
### Aliases: glmb.influence.measures rstandard.glmb rstudent.glmb
###   glmb.dffits dfbetas.glmb glmb.covratio cooks.distance.glmb

### ** Examples

set.seed(333)
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
mysd<-1
mu<-matrix(0,5)
mu[1]=log(mean(counts))
V0<-((mysd)^2)*diag(5)
glmb.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu,Sigma=V0))

## Start setup here [First output from earlier optim optimization]
betastar=glmb.D93$coef.mode  # Posterior mode from optim
x=glmb.D93$x
y=glmb.D93$y
#offset=glmb.D93$offset   # not present in the output --> For now set to 0 vector
offset2=0*y   # Should return this from lower level functions
weights2=glmb.D93$prior.weights

## Check influence measures for original model
fit=glmb.wfit(x,y,weights2,offset2,family=poisson(),Bbar=mu,P=solve(V0),betastar)
influence.measures(fit)

print(fit)
print(glmb.D93$coef.mode)

### Now try a strong prior with poorly chosen intercept
mu1=0*mu
V1=0.1*V0
glmb2.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu1,Sigma=V1))

Bbar2=mu1  # Prior mean
betastar2=glmb2.D93$coef.mode  # Posterior mode from optim
fit2=glmb.wfit(x,y,weights2,offset2,family=poisson(),Bbar2,P=solve(V1),betastar2)

influence.measures(fit2)

print(fit2)
print(glmb2.D93$coef.mode)



influence(glmb2.D93)
influence(glmb2.D93,do.coef=TRUE)
influence(glmb2.D93,do.coef=FALSE)

#glmb2.D93$qr=glmb2.D93$fit$qr
#influence.measures(glmb2.D93,influence(glmb2.D93))
#influence.measures(glmb2.D93$fit,influence(glmb2.D93))
#influence.measures(glmb2.D93$fit,influence(glmb2.D93,do.coef=FALSE))
glmb.influence.measures(glmb2.D93,influence(glmb2.D93))


## Do not need method functions
dfbeta(glmb2.D93,influence(glmb2.D93))
dfbeta(glmb2.D93$fit,influence(glmb2.D93))
hatvalues(glmb2.D93,influence(glmb2.D93))
hatvalues(glmb2.D93$fit,influence(glmb2.D93))

# Methods (implmented)

rstandard(glmb2.D93,infl=influence(glmb2.D93))
rstandard(glmb2.D93)


# Need Metods

## Needs a method function
#dfbetas(glmb2.D93,influence(glmb2.D93))
dfbetas(glmb2.D93)
dfbetas(glmb2.D93,influence(glmb2.D93,do.coef=TRUE))
dfbetas(glmb2.D93$fit,influence(glmb2.D93))




# Needs a method function
cooks.distance(glmb2.D93)
cooks.distance(glmb2.D93$fit,influence(glmb2.D93))



# Needs a method function (now works)
#rstandard(glmb2.D93$fit,influence(glmb2.D93))


# Needs a method function
rstudent(glmb2.D93)
rstudent(glmb2.D93$fit,influence(glmb2.D93))


# Not a method, separate function
glmb.dffits(glmb2.D93)
dffits(glmb2.D93$fit,influence(glmb2.D93)) 

# Not a method - separate function
glmb.covratio(glmb2.D93)  
covratio(glmb2.D93$fit,influence(glmb2.D93))  



# Needs a method function but requires a different approach
# The influence function actually stores this measure
## This seems to require more work to create a method function
## Should 
#hat(glmb2.D93,influence(glmb2.D93))
#hat(glmb2.D93$fit,influence(glmb2.D93))

infl=influence(glmb2.D93)
hat2=infl$hat               
hat2               



cleanEx()
nameEx("glmb.wfit")
### * glmb.wfit

flush(stderr()); flush(stdout())

### Name: glmb.wfit
### Title: Fitter Function for Bayesian Generalized Linear Models
### Aliases: glmb.wfit

### ** Examples

set.seed(333)
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
mysd<-1
mu<-matrix(0,5)
mu[1]=log(mean(counts))
V0<-((mysd)^2)*diag(5)
glmb.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu,Sigma=V0))

## Start setup here [First output from earlier optim optimization]
betastar=glmb.D93$coef.mode  # Posterior mode from optim
x=glmb.D93$x
y=glmb.D93$y
#offset=glmb.D93$offset   # not present in the output --> For now set to 0 vector
offset2=0*y   # Should return this from lower level functions
weights2=glmb.D93$prior.weights

## Check influence measures for original model
fit=glmb.wfit(x,y,weights2,offset2,family=poisson(),Bbar=mu,P=solve(V0),betastar)
influence.measures(fit)

print(fit)
print(glmb.D93$coef.mode)

### Now try a strong prior with poorly chosen intercept
mu1=0*mu
V1=0.1*V0
glmb2.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu1,Sigma=V1))

Bbar2=mu1  # Prior mean
betastar2=glmb2.D93$coef.mode  # Posterior mode from optim
fit2=glmb.wfit(x,y,weights2,offset2,family=poisson(),Bbar2,P=solve(V1),betastar2)

influence.measures(fit2)

print(fit2)
print(glmb2.D93$coef.mode)



influence(glmb2.D93)
influence(glmb2.D93,do.coef=TRUE)
influence(glmb2.D93,do.coef=FALSE)

#glmb2.D93$qr=glmb2.D93$fit$qr
#influence.measures(glmb2.D93,influence(glmb2.D93))
#influence.measures(glmb2.D93$fit,influence(glmb2.D93))
#influence.measures(glmb2.D93$fit,influence(glmb2.D93,do.coef=FALSE))
glmb.influence.measures(glmb2.D93,influence(glmb2.D93))


## Do not need method functions
dfbeta(glmb2.D93,influence(glmb2.D93))
dfbeta(glmb2.D93$fit,influence(glmb2.D93))
hatvalues(glmb2.D93,influence(glmb2.D93))
hatvalues(glmb2.D93$fit,influence(glmb2.D93))

# Methods (implmented)

rstandard(glmb2.D93,infl=influence(glmb2.D93))
rstandard(glmb2.D93)


# Need Metods

## Needs a method function
#dfbetas(glmb2.D93,influence(glmb2.D93))
dfbetas(glmb2.D93)
dfbetas(glmb2.D93,influence(glmb2.D93,do.coef=TRUE))
dfbetas(glmb2.D93$fit,influence(glmb2.D93))




# Needs a method function
cooks.distance(glmb2.D93)
cooks.distance(glmb2.D93$fit,influence(glmb2.D93))



# Needs a method function (now works)
#rstandard(glmb2.D93$fit,influence(glmb2.D93))


# Needs a method function
rstudent(glmb2.D93)
rstudent(glmb2.D93$fit,influence(glmb2.D93))


# Not a method, separate function
glmb.dffits(glmb2.D93)
dffits(glmb2.D93$fit,influence(glmb2.D93)) 

# Not a method - separate function
glmb.covratio(glmb2.D93)  
covratio(glmb2.D93$fit,influence(glmb2.D93))  



# Needs a method function but requires a different approach
# The influence function actually stores this measure
## This seems to require more work to create a method function
## Should 
#hat(glmb2.D93,influence(glmb2.D93))
#hat(glmb2.D93$fit,influence(glmb2.D93))

infl=influence(glmb2.D93)
hat2=infl$hat               
hat2               



cleanEx()
nameEx("glmb_Standardize_Model")
### * glmb_Standardize_Model

flush(stderr()); flush(stdout())

### Name: glmb_Standardize_Model
### Title: Standardize A Non-Gaussian Model
### Aliases: glmb_Standardize_Model
### Keywords: internal

### ** Examples

data(menarche2)

summary(menarche2)
plot(Menarche/Total ~ Age, data=menarche2)

Age2=menarche2$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche2$Menarche/menarche2$Total
wt=menarche2$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 1 times as large as point estimate 

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



###### Adjust weight for dispersion

wt2=wt/dispersion2

######################### Shift mean vector to offset so that adjusted model has 0 mean

alpha=x%*%as.vector(mu)+offset2
mu2=0*as.vector(mu)
P2=P
x2=x


#####  Optimization step to find posterior mode and associated Precision

parin=start-mu

opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x),mu=as.vector(mu2),
              P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
              method="BFGS",hessian=TRUE
)

bstar=opt_out$par  ## Posterior mode for adjusted model
bstar
bstar+as.vector(mu)  # mode for actual model
A1=opt_out$hessian # Approximate Precision at mode

## Standardize Model

Standard_Mod=glmb_Standardize_Model(y=as.vector(y), x=as.matrix(x),P=as.matrix(P),
                                    bstar=as.matrix(bstar,ncol=1), A1=as.matrix(A1))

bstar2=Standard_Mod$bstar2  
A=Standard_Mod$A
x2=Standard_Mod$x2
mu2=Standard_Mod$mu2
P2=Standard_Mod$P2
L2Inv=Standard_Mod$L2Inv
L3Inv=Standard_Mod$L3Inv



cleanEx()
nameEx("glmbayes")
### * glmbayes

flush(stderr()); flush(stdout())

### Name: glmbayes
### Title: Bayesian Generalized Linear Models (iid Samples)
### Aliases: glmbayes PackageOverview
### Keywords: internal

### ** Examples

## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

## Call to glm
glm.D93 <- glm(counts ~ outcome + treatment, 
               family = poisson())

## Using glmb
## Step 1: Set up Prior
ps=Prior_Setup(counts ~ outcome + treatment)
mu=ps$mu
V=ps$Sigma
# Step2A: Check the Prior
Prior_Check(counts ~ outcome + treatment,family = poisson(),
            pfamily=dNormal(mu=mu,Sigma=V))
# Step2B: Update and Re-Check the Prior
mu[1,1]=log(mean(counts))
Prior_Check(counts ~ outcome + treatment,family = poisson(),
            pfamily=dNormal(mu=mu,Sigma=V))
# Step 3: Call the glmb function
glmb.D93<-glmb(counts ~ outcome + treatment, family=poisson(), 
               pfamily=dNormal(mu=mu,Sigma=V))

## ----Printed_Views------------------------------------------------------------
## Printed view of the output from the glm function 
print(glm.D93)
## Printed view of the output from the glmb function 
print(glmb.D93)

## ----Methods---------------------------------------------------------------
## Methods for class "lm"
methods(class="lm")

## Methods for class "glm"
methods(class="glm")

## Methods for class "glmb"
methods(class="glmb")

## ----summary--------------------------------------------------------------
## summary output for the "glm" class
summary(glm.D93)

## summary output for the "glm" class
summary(glmb.D93)

## ----fitted outputs-------------------------------------------------------
## fitted outputs for the glm function
fitted(glm.D93)

## ----glmb fitted outputs------------------------------------------------------
## mean of fitted outputs for the glm function
colMeans(fitted(glmb.D93))

## ----predictions----------------------------------------------------------
## predictions for the glm function
predict(glm.D93)

## predictions for the glmb function
colMeans(predict(glmb.D93)) 

## ----residuals------------------------------------------------------------
## residuals for the glm function
residuals(glm.D93)

## residuals for the glmb function
colMeans(residuals(glmb.D93))

## ----vcov-----------------------------------------------------------------
## vcov for the glm function
vcov(glm.D93)

## vcov for the glm function
vcov(glmb.D93)

## ----confint--------------------------------------------------------------
## confint for the glm function
confint(glm.D93)

## confint for the glm function
confint(glmb.D93)

## ----AIC/DIC------------------------------------------------------------------
## AIC for the glm function (equivalent degrees of freedom and the AIC)
extractAIC(glm.D93)

## DIC for the glmb function (estimated effective number of parameters and the DIC)
extractAIC(glmb.D93)

## ----Deviance-------------------------------------------------------------
## Deviance for the glm function
deviance(glm.D93)

## Deviance for the glmb function
mean(deviance(glmb.D93))

## ----logLik---------------------------------------------------------------
## Deviance for the glm function
logLik(glm.D93)

## Deviance for the glmb function
mean(logLik(glmb.D93))

## ----Model Frame----------------------------------------------------------
## Model Frame for the glm function
model.frame(glm.D93)

## Model Frame for the glmb function
model.frame(glmb.D93$glm)

## ----formula--------------------------------------------------------------
## formula for the glm function
formula(glm.D93)

## ----formula-------------------------------------------------------------
## formula for the glmb function
formula(glmb.D93)

## ----family--------------------------------------------------------------
## family for the glm function
family(glm.D93)

## family for the glmb function
family(glmb.D93$glm)

## ----nobs-----------------------------------------------------------------
## nobs for the glm function
nobs(glm.D93)

## nobs for the glmb function
nobs(glmb.D93)

## ----show-----------------------------------------------------------------
## show for the glm function
show(glm.D93)

## show for the glmb function
show(glmb.D93)



cleanEx()
nameEx("glmbfamfunc")
### * glmbfamfunc

flush(stderr()); flush(stdout())

### Name: glmbfamfunc
### Title: Return family functions used during simulation and post
###   processing
### Aliases: glmbfamfunc print.glmbfamfunc
### Keywords: internal

### ** Examples

famfunc<-glmbfamfunc(binomial(logit))

print(famfunc)

f1<-famfunc$f1
f2<-famfunc$f2
f3<-famfunc$f3
f5<-famfunc$f5
f6<-famfunc$f6



cleanEx()
nameEx("influence.glmb")
### * influence.glmb

flush(stderr()); flush(stdout())

### Name: influence.glmb
### Title: Bayesian Regression Diagnostics
### Aliases: influence.glmb

### ** Examples

set.seed(333)
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
mysd<-1
mu<-matrix(0,5)
mu[1]=log(mean(counts))
V0<-((mysd)^2)*diag(5)
glmb.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu,Sigma=V0))

## Start setup here [First output from earlier optim optimization]
betastar=glmb.D93$coef.mode  # Posterior mode from optim
x=glmb.D93$x
y=glmb.D93$y
#offset=glmb.D93$offset   # not present in the output --> For now set to 0 vector
offset2=0*y   # Should return this from lower level functions
weights2=glmb.D93$prior.weights

## Check influence measures for original model
fit=glmb.wfit(x,y,weights2,offset2,family=poisson(),Bbar=mu,P=solve(V0),betastar)
influence.measures(fit)

print(fit)
print(glmb.D93$coef.mode)

### Now try a strong prior with poorly chosen intercept
mu1=0*mu
V1=0.1*V0
glmb2.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu1,Sigma=V1))

Bbar2=mu1  # Prior mean
betastar2=glmb2.D93$coef.mode  # Posterior mode from optim
fit2=glmb.wfit(x,y,weights2,offset2,family=poisson(),Bbar2,P=solve(V1),betastar2)

influence.measures(fit2)

print(fit2)
print(glmb2.D93$coef.mode)



influence(glmb2.D93)
influence(glmb2.D93,do.coef=TRUE)
influence(glmb2.D93,do.coef=FALSE)

#glmb2.D93$qr=glmb2.D93$fit$qr
#influence.measures(glmb2.D93,influence(glmb2.D93))
#influence.measures(glmb2.D93$fit,influence(glmb2.D93))
#influence.measures(glmb2.D93$fit,influence(glmb2.D93,do.coef=FALSE))
glmb.influence.measures(glmb2.D93,influence(glmb2.D93))


## Do not need method functions
dfbeta(glmb2.D93,influence(glmb2.D93))
dfbeta(glmb2.D93$fit,influence(glmb2.D93))
hatvalues(glmb2.D93,influence(glmb2.D93))
hatvalues(glmb2.D93$fit,influence(glmb2.D93))

# Methods (implmented)

rstandard(glmb2.D93,infl=influence(glmb2.D93))
rstandard(glmb2.D93)


# Need Metods

## Needs a method function
#dfbetas(glmb2.D93,influence(glmb2.D93))
dfbetas(glmb2.D93)
dfbetas(glmb2.D93,influence(glmb2.D93,do.coef=TRUE))
dfbetas(glmb2.D93$fit,influence(glmb2.D93))




# Needs a method function
cooks.distance(glmb2.D93)
cooks.distance(glmb2.D93$fit,influence(glmb2.D93))



# Needs a method function (now works)
#rstandard(glmb2.D93$fit,influence(glmb2.D93))


# Needs a method function
rstudent(glmb2.D93)
rstudent(glmb2.D93$fit,influence(glmb2.D93))


# Not a method, separate function
glmb.dffits(glmb2.D93)
dffits(glmb2.D93$fit,influence(glmb2.D93)) 

# Not a method - separate function
glmb.covratio(glmb2.D93)  
covratio(glmb2.D93$fit,influence(glmb2.D93))  



# Needs a method function but requires a different approach
# The influence function actually stores this measure
## This seems to require more work to create a method function
## Should 
#hat(glmb2.D93,influence(glmb2.D93))
#hat(glmb2.D93$fit,influence(glmb2.D93))

infl=influence(glmb2.D93)
hat2=infl$hat               
hat2               



cleanEx()
nameEx("lmb")
### * lmb

flush(stderr()); flush(stdout())

### Name: lmb
### Title: Fitting Bayesian Linear Models
### Aliases: lmb print.lmb

### ** Examples

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

lm.D9 <- lm(weight ~ group,x=TRUE,y=TRUE)
disp_ML=sigma(lm.D9)^2
n_prior=2
shape=n_prior/2
rate= disp_ML*shape

# Conjugate Normal_Gamma Prior 
lmb.D9=lmb(weight ~ group,dNormal_Gamma(mu,V/disp_ML,shape=shape,rate=rate))
summary(lmb.D9)

# Independent_Normal_Gamma_Prior
lmb.D9_v2=lmb(weight ~ group,dIndependent_Normal_Gamma(mu,V,shape=shape,rate=rate))
summary(lmb.D9_v2)



cleanEx()
nameEx("logLik.glmb")
### * logLik.glmb

flush(stderr()); flush(stdout())

### Name: logLik.glmb
### Title: Extract Log-Likelihood
### Aliases: logLik.glmb

### ** Examples

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb logLik-------------------------------------------------------------
colMeans(logLik(glmb.D93))



cleanEx()
nameEx("menarche")
### * menarche

flush(stderr()); flush(stdout())

### Name: menarche2
### Title: Age Of Menarche In Warsaw
### Aliases: menarche2
### Keywords: datasets, Bayesian Binomial Regression

### ** Examples

data(menarche2)



cleanEx()
nameEx("pfamily")
### * pfamily

flush(stderr()); flush(stdout())

### Name: pfamily
### Title: Prior Family Objects for Bayesian Models
### Aliases: pfamily dNormal dGamma dNormal_Gamma dIndependent_Normal_Gamma
###   print.pfamily

### ** Examples

mu=c(0,0)
Sigma=diag(2)

npf<-dNormal(mu,Sigma)  # Normal pfamily
str(dNormal(mu,Sigma))

## Example where # Normal pfamily is used

## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
mysd<-1
mu<-matrix(0,5)
mu[1]=log(mean(counts))
V0<-((mysd)^2)*diag(5)
glmb.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu,Sigma=V0))
summary(glmb.D93)



cleanEx()
nameEx("predict.glmb")
### * predict.glmb

flush(stderr()); flush(stdout())

### Name: predict.glmb
### Title: Predict Method for Bayesian GLM Fits
### Aliases: predict.glmb

### ** Examples

data(menarche2)
## ----Analysis Setup-----------------------------------------------------------
## Number of variables in model
Age=menarche2$Age
nvars=2
## Reference Ages for setting of priors and Age_Difference
ref_age1=13  # user can modify this
ref_age2=15  ## user can modify this
## Define variables used later in analysis
Age2=Age-ref_age1
Age_Diff=ref_age2-ref_age1
mu1=as.matrix(c(0,1.098612),ncol=1)
V1<-1*diag(nvars)
V1[1,1]=0.18687882
V1[2,2]=0.10576217
V1[1,2]=-0.03389182
V1[2,1]=-0.03389182
Menarche_Model_Data=data.frame(Age=menarche2$Age,Total=menarche2$Total,
                               Menarche=menarche2$Menarche,Age2)
glmb.out1<-glmb(n=1000,cbind(Menarche, Total-Menarche) ~Age2,family=binomial(logit),
                pfamily=dNormal(mu=mu1,Sigma=V1),data=Menarche_Model_Data)

# Prediction from original model
pred0=predict(glmb.out1,type="link")
colMeans(pred0)

pred1=predict(glmb.out1,type="response")
colMeans(pred1)


## Generate new data
Age_New <- seq(8, 20, 0.25)
Age2_New=Age_New-13
mod_Object=glmb.out1
obs_size=median(menarche2$Total) ## Counts for sim from Binomial
olddata=data.frame(Age=menarche2$Age,
                   Menarche=menarche2$Menarche,Total=menarche2$Total,Age2=Age2)

newdata=data.frame(Age=Age_New,Age2=Age2_New)

# Simulate for newdata

pred_menarche=predict(mod_Object,newdata=newdata,olddata=olddata,type="response")
pred_m=colMeans(pred_menarche)

n=nrow(mod_Object$coefficients)
pred_y=matrix(0,nrow=n,ncol=length(Age_New))
for(i in 1:n){
  pred_y[i,1:length(Age_New)]=rbinom(length(Age_New),size=obs_size,
  prob=pred_menarche[i,1:length(Age_New)
                                                                                      ])
}

# Produce various predictions
pred_y_m=colMeans(pred_y/obs_size)
quant1_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.025))
quant2_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.975))
quant1_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.025))
quant2_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.975))

#Plot Predictions for newdata
plot(Menarche/Total ~ Age, data=menarche2,
main="Percentage of girls Who have had their first period")
lines(Age_New, pred_m,lty=1)
lines(Age_New, quant1_m,lty=2)
lines(Age_New, quant2_m,lty=2)
lines(Age_New, quant1_m_y,lty=2)
lines(Age_New, quant2_m_y,lty=2)




cleanEx()
nameEx("rGamma_reg")
### * rGamma_reg

flush(stderr()); flush(stdout())

### Name: rGamma_reg
### Title: The Bayesian Generalized Linear Model Dispersion Distribution
### Aliases: rGamma_reg print.rGamma_reg summary.rGamma_reg

### ** Examples

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

prior_list=list(beta=b_old,shape=shape,rate=rate)

dispout<-rGamma_reg(n=n,y,x,prior_list=prior_list,
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
## Could use a residuals function here -- For now, maybe run the glmb function

summary(outtemp3)
mean(colMeans(residuals(outtemp3)^2))
v_old
colMeans((outtemp3$coefficients))
b_old
mean(outtemp3$dispersion)  ## Seems slightly smaller --> Needs qc
v_old



cleanEx()
nameEx("rNormal_Gamma_reg")
### * rNormal_Gamma_reg

flush(stderr()); flush(stdout())

### Name: rNormal_Gamma_reg
### Title: The Bayesian Normal-Gamma Regression Distribution
### Aliases: rNormal_Gamma_reg

### ** Examples

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
ngamma.D9=rNormal_Gamma_reg(n=1000,y=y,x=x,
  prior_list=prior_list)
summary(ngamma.D9$coefficients)




cleanEx()
nameEx("rNormal_reg")
### * rNormal_reg

flush(stderr()); flush(stdout())

### Name: rNormal_reg
### Title: The Bayesian Generalized Linear Model with Normal Prior
###   Distribution
### Aliases: rNormal_reg

### ** Examples

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
ngamma.D9=rNormal_Gamma_reg(n=1000,y=y,x=x,
  prior_list=prior_list)
summary(ngamma.D9$coefficients)




cleanEx()
nameEx("rNormal_reg.wfit")
### * rNormal_reg.wfit

flush(stderr()); flush(stdout())

### Name: rNormal_reg.wfit
### Title: Fitter Function for Bayesian Linear Models
### Aliases: rNormal_reg.wfit

### ** Examples


## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb confint-------------------------------------------------------------
confint(glmb.D93)



cleanEx()
nameEx("residuals.glmb")
### * residuals.glmb

flush(stderr()); flush(stdout())

### Name: residuals.glmb
### Title: Accessing Bayesian Generalized Linear Model Fits
### Aliases: residuals.glmb residuals.lmb

### ** Examples

data(menarche2)
## ----Analysis Setup-----------------------------------------------------------
## Number of variables in model
Age=menarche2$Age
nvars=2
## Reference Ages for setting of priors and Age_Difference
ref_age1=13  # user can modify this
ref_age2=15  ## user can modify this
## Define variables used later in analysis
Age2=menarche2$Age-ref_age1
Age_Diff=ref_age2-ref_age1
mu1=as.matrix(c(0,1.098612),ncol=1)
V1<-1*diag(nvars)
V1[1,1]=0.18687882
V1[2,2]=0.10576217
V1[1,2]=-0.03389182
V1[2,1]=-0.03389182
Menarche_Model_Data=data.frame(Age=menarche2$Age,Total=menarche2$Total,
                               Menarche=menarche2$Menarche,Age2)
glmb.out1<-glmb(n=1000,cbind(Menarche, Total-Menarche) ~Age2,family=binomial(logit),
                pfamily=dNormal(mu=mu1,Sigma=V1),data=Menarche_Model_Data)

# Prediction from original model
pred1=predict(glmb.out1,type="response")

## Get Original Residuals, their means, and credible bounds
res_out=residuals(glmb.out1)
colMeans(res_out)

## Set up to simulate new data and residuals
res_mean=colMeans(res_out)
res_low1=apply(res_out,2,FUN=quantile,probs=c(0.025))
res_high1=apply(res_out,2,FUN=quantile,probs=c(0.975))

## Simulate new data and get residuals for simulated data

ysim1=simulate(glmb.out1,nsim=1,seed=NULL,pred=pred1,family="binomial",
               prior.weights=weights(glmb.out1))


res_ysim_out1=residuals(glmb.out1,ysim=ysim1)
res_low=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.025))
res_high=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.975))

# Plot Credible Interval bounds for Deviance Residuals

plot(res_mean~Age,ylim=c(-2.5,2.5),
main="Credible Interval Bound for Menarche - Logit Model Deviance Residuals",
xlab = "Age", ylab = "Avg. Dev. Res")
lines(Age, 0*res_mean,lty=1)
lines(Age, res_low,lty=1)
lines(Age, res_high,lty=1)
lines(Age, res_low1,lty=2)
lines(Age, res_high1,lty=2)




cleanEx()
nameEx("rglmb")
### * rglmb

flush(stderr()); flush(stdout())

### Name: rglmb
### Title: The Bayesian Generalized Linear Model Distribution
### Aliases: rglmb print.rglmb

### ** Examples

data(menarche2)

summary(menarche2)
plot(Menarche/Total ~ Age, data=menarche2)

Age2=menarche2$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche2$Menarche/menarche2$Total
wt=menarche2$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 

out<-rglmb(n = 1000, y=y, x=x, pfamily=dNormal(mu=mu,Sigma=V1), weights = wt, 
           family = binomial(logit)) 
summary(out)

# Add mean(out$iters to rglmb summary function)
mean(out$iters)



cleanEx()
nameEx("rindep_norm_gamma_reg_std_R")
### * rindep_norm_gamma_reg_std_R

flush(stderr()); flush(stdout())

### Name: rindep_norm_gamma_reg_std_R
### Title: The Standard Bayesian Indendepent Normal-Gamma Regression
###   Distribution
### Aliases: rindep_norm_gamma_reg_std_R

### ** Examples


## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group,y=TRUE,x=TRUE)

p_setup=Prior_Setup(lm.D9)
dispersion=sigma(lm.D9)^2

## Initialize dispersion
p=1/dispersion

mu=p_setup$mu
Sigma_prior=p_setup$Sigma
y=lm.D9$y
x=lm.D9$x

Prior_Check(weight ~ group,gaussian(),dNormal(mu=mu,Sigma=Sigma_prior))

## This simple "standardization" for the prior appears to 
## correct the issue with the prior

## "Standardize" the prior for the intercept by setting the
## mean equal to the mean of the dependent variable
## and setting the variance equal to the variance of the dependent variable

mu[1,1]=mean(y)
Sigma_prior[1,1]=var(y)

## For factors, set the "standad prior to mean=0 and variance driven by the range
## of the dependent variable 

mu[2,1]=0
Sigma_prior[2,2]=((max(y)-min(y))/1.96)^2

#Prior_Check(lm.D9,mu,Sigma_prior)
Prior_Check(weight ~ group,gaussian(),dNormal(mu=mu,Sigma=Sigma_prior))


#prior=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion)
glmb.D9=glmb(weight~group, family=gaussian(),dNormal(mu=mu,Sigma=Sigma_prior,dispersion=dispersion))
post_mode=glmb.D9$coef.mode

sum_out1=summary(glmb.D9)

## Try mean one standard deviation away to see how it works
#mu[1,1]=mu[1,1]-2*sum_out1$coefficients[1,3]
#mu[2,1]=mu[2,1]+2*sum_out1$coefficients[1,3]

#prior=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion)

# Temporarily lower the prior variance
Sigma_prior=1*Sigma_prior

glmb.D9_v2=glmb(n=1000,weight~group, family=gaussian(),
dNormal(mu=mu,Sigma=Sigma_prior,dispersion=dispersion))

n_prior=2
shape=n_prior/2
rate= dispersion*shape


summary(glmb.D9)
summary(glmb.D9_v2)


#lm_out=lm(y ~ x-1) # returns same model
RSS=sum(residuals(lm.D9)^2)

n_prior=4
n_data=length(y)
shape=(n_prior/2)
rate=n_prior*RSS/n_data



prior_list=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion,
                 shape=shape,rate=rate,Precision=solve(Sigma_prior))


set.seed(360)

 ptm <- proc.time()
 sim2=rindependent_norm_gamma_reg(n=1000,y,x,prior_list=prior_list,
offset=NULL,weights=1,max_disp_perc=0.99)
 proc.time()-ptm

 
 
 
 
 



cleanEx()
nameEx("rindependent_norm_gamma_reg")
### * rindependent_norm_gamma_reg

flush(stderr()); flush(stdout())

### Name: rindependent_norm_gamma_reg
### Title: The Bayesian Indendepent Normal-Gamma Regression Distribution
### Aliases: rindependent_norm_gamma_reg rindependent_norm_gamma_reg_v3
###   rindependent_norm_gamma_reg_v2 rindependent_norm_gamma_reg_v4

### ** Examples


## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group,y=TRUE,x=TRUE)

p_setup=Prior_Setup(lm.D9)
dispersion=sigma(lm.D9)^2

## Initialize dispersion
p=1/dispersion

mu=p_setup$mu
Sigma_prior=p_setup$Sigma
y=lm.D9$y
x=lm.D9$x

Prior_Check(weight ~ group,gaussian(),dNormal(mu=mu,Sigma=Sigma_prior))

## This simple "standardization" for the prior appears to 
## correct the issue with the prior

## "Standardize" the prior for the intercept by setting the
## mean equal to the mean of the dependent variable
## and setting the variance equal to the variance of the dependent variable

mu[1,1]=mean(y)
Sigma_prior[1,1]=var(y)

## For factors, set the "standad prior to mean=0 and variance driven by the range
## of the dependent variable 

mu[2,1]=0
Sigma_prior[2,2]=((max(y)-min(y))/1.96)^2

#Prior_Check(lm.D9,mu,Sigma_prior)
Prior_Check(weight ~ group,gaussian(),dNormal(mu=mu,Sigma=Sigma_prior))


#prior=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion)
glmb.D9=glmb(weight~group, family=gaussian(),dNormal(mu=mu,Sigma=Sigma_prior,dispersion=dispersion))
post_mode=glmb.D9$coef.mode

sum_out1=summary(glmb.D9)

## Try mean one standard deviation away to see how it works
#mu[1,1]=mu[1,1]-2*sum_out1$coefficients[1,3]
#mu[2,1]=mu[2,1]+2*sum_out1$coefficients[1,3]

#prior=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion)

# Temporarily lower the prior variance
Sigma_prior=1*Sigma_prior

glmb.D9_v2=glmb(n=1000,weight~group, family=gaussian(),
dNormal(mu=mu,Sigma=Sigma_prior,dispersion=dispersion))

n_prior=2
shape=n_prior/2
rate= dispersion*shape


summary(glmb.D9)
summary(glmb.D9_v2)


#lm_out=lm(y ~ x-1) # returns same model
RSS=sum(residuals(lm.D9)^2)

n_prior=4
n_data=length(y)
shape=(n_prior/2)
rate=n_prior*RSS/n_data



prior_list=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion,
                 shape=shape,rate=rate,Precision=solve(Sigma_prior))


set.seed(360)

 ptm <- proc.time()
 sim2=rindependent_norm_gamma_reg(n=1000,y,x,prior_list=prior_list,
offset=NULL,weights=1,max_disp_perc=0.99)
 proc.time()-ptm

 
 
 
 
 



cleanEx()
nameEx("rlmb")
### * rlmb

flush(stderr()); flush(stdout())

### Name: rlmb
### Title: The Bayesian Linear Model Distribution
### Aliases: rlmb rlmb.print print.rlmb

### ** Examples

## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)

ps=Prior_Setup(weight ~ group)
x=ps$x
mu=ps$mu
V=ps$Sigma
mu[1,1]=mean(weight)

Prior_Check(weight ~ group,family =gaussian(),
            pfamily=dNormal(mu=mu,Sigma=V))
lm.D9 <- lm(weight ~ group,x=TRUE,y=TRUE)
y=lm.D9$y

## Dispersion for maximum likelihood estimate
disp_ML=sigma(lm.D9)^2
n_prior=2
shape=n_prior/2
rate= disp_ML*shape

# Two-Block Gibbs sampler
set.seed(180)
dispersion2=disp_ML
# Run 1000 burn-in iterations
for(i in 1:1000){
  out1=rlmb(n = 1, y=y, x=x, pfamily=dNormal(mu=mu,Sigma=V,dispersion=dispersion2))
  out2=rlmb(n = 1, y=y, x=x, pfamily=dGamma(shape=shape,rate=rate,beta=out1$coefficients[1,]))
  dispersions=out2$dispersion
}


# Run 1000 additional iterations and store output
beta_out<-matrix(0,nrow=1000, ncol=2)
disp_out=rep(0,1000)
for(i in 1:1000){
  out1=rlmb(n = 1, y=y, x=x, pfamily=dNormal(mu=mu,Sigma=V,dispersion=dispersion2))
  out2=rlmb(n = 1, y=y, x=x, pfamily=dGamma(shape=shape,rate=rate,beta=out1$coefficients[1,]))
  dispersions=out2$dispersion
  beta_out[i,1:2]=out1$coefficients[1,1:2]
  disp_out[i]=out2$dispersion
}

colMeans(beta_out)
mean(disp_out)

# Same model using Independent_Normal_Gamma_Prior
lmb.D9_v2=lmb(weight ~ group,dIndependent_Normal_Gamma(mu,V,shape=shape,rate=rate))
summary(lmb.D9_v2)





cleanEx()
nameEx("rnnorm_reg_std")
### * rnnorm_reg_std

flush(stderr()); flush(stdout())

### Name: rnnorm_reg_std
### Title: The Bayesian Generalized Linear Model Distribution in Standard
###   Form
### Aliases: rnnorm_reg_std
### Keywords: internal

### ** Examples

data(menarche2)

summary(menarche2)
plot(Menarche/Total ~ Age, data=menarche2)

Age2=menarche2$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche2$Menarche/menarche2$Total
wt=menarche2$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 1 times as large as point estimate 

#out<-rglmb(n = 1000, y=y, x=x, mu=mu, P=solve(V1), wt = wt, 
#           family = binomial(logit), Gridtype = 3) 
#summary(out)

famfunc<-glmbfamfunc(binomial(logit))

f1<-famfunc$f1
f2<-famfunc$f2  # Used in optim and glmbsim_cpp
f3<-famfunc$f3  # Used in optim
f5<-famfunc$f5
f6<-famfunc$f6

dispersion2<-as.numeric(1.0)
start <- mu
offset2=rep(as.numeric(0.0),length(y))
P=solve(V1)
n=1000


###### Adjust weight for dispersion

wt2=wt/dispersion2

######################### Shift mean vector to offset so that adjusted model has 0 mean

alpha=x%*%as.vector(mu)+offset2
mu2=0*as.vector(mu)
P2=P
x2=x


#####  Optimization step to find posterior mode and associated Precision

parin=start-mu

opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x),mu=as.vector(mu2),
              P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
              method="BFGS",hessian=TRUE
)

bstar=opt_out$par  ## Posterior mode for adjusted model
bstar
bstar+as.vector(mu)  # mode for actual model
A1=opt_out$hessian # Approximate Precision at mode

## Standardize Model

Standard_Mod=glmb_Standardize_Model(y=as.vector(y), x=as.matrix(x),P=as.matrix(P),
                                    bstar=as.matrix(bstar,ncol=1), A1=as.matrix(A1))

bstar2=Standard_Mod$bstar2  
A=Standard_Mod$A
x2=Standard_Mod$x2
mu2=Standard_Mod$mu2
P2=Standard_Mod$P2
L2Inv=Standard_Mod$L2Inv
L3Inv=Standard_Mod$L3Inv

Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
                   as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),as.vector(wt2),
                   family="binomial",link="logit",Gridtype=as.integer(3), 
                   n=as.integer(n),sortgrid=TRUE)

## These now seem to match

Env2

#int n, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, 
#NumericVector alpha, NumericVector wt, Function f2,
#Rcpp::List Envelope, Rcpp::CharacterVector family, 
#Rcpp::CharacterVector link, int progbar

### Note: getting the types correct here is important but potentially difficult for users
### May be better to call an R function wrapper that checks (and converts when possible) 
##  to correct types

sim=rnnorm_reg_std(n=as.integer(n),y=as.vector(y),x=as.matrix(x2),mu=as.matrix(mu2,ncol=1),
                   P=as.matrix(P2),alpha=as.vector(alpha),wt=as.vector(wt2),
                   f2=f2,Envelope=Env2,family="binomial",link="logit",as.integer(0))

out=L2Inv%*%L3Inv%*%t(sim$out)

for(i in 1:n){
  out[,i]=out[,i]+mu
}

summary(t(out))
mean(sim$draws)



cleanEx()
nameEx("setlogP")
### * setlogP

flush(stderr()); flush(stdout())

### Name: setlogP
### Title: Calculate constants used during sampling from Likelihood
###   subgradient densities
### Aliases: setlogP
### Keywords: internal

### ** Examples

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb extractAIC-------------------------------------------------------------
extractAIC(glmb.D93)



cleanEx()
nameEx("simulate.glmb")
### * simulate.glmb

flush(stderr()); flush(stdout())

### Name: simulate.glmb
### Title: Simulate Responses
### Aliases: simulate.glmb

### ** Examples

data(menarche2)
## ----Analysis Setup-----------------------------------------------------------
## Number of variables in model
Age=menarche2$Age
nvars=2
## Reference Ages for setting of priors and Age_Difference
ref_age1=13  # user can modify this
ref_age2=15  ## user can modify this
## Define variables used later in analysis
Age2=menarche2$Age-ref_age1
Age_Diff=ref_age2-ref_age1
mu1=as.matrix(c(0,1.098612),ncol=1)
V1<-1*diag(nvars)
V1[1,1]=0.18687882
V1[2,2]=0.10576217
V1[1,2]=-0.03389182
V1[2,1]=-0.03389182
Menarche_Model_Data=data.frame(Age=menarche2$Age,Total=menarche2$Total,
                               Menarche=menarche2$Menarche,Age2)
glmb.out1<-glmb(n=1000,cbind(Menarche, Total-Menarche) ~Age2,family=binomial(logit),
                pfamily=dNormal(mu=mu1,Sigma=V1),data=Menarche_Model_Data)

# Prediction from original model
pred1=predict(glmb.out1,type="response")

## Get Original Residuals, their means, and credible bounds
res_out=residuals(glmb.out1)
colMeans(res_out)

## Set up to simulate new data and residuals
res_mean=colMeans(res_out)
res_low1=apply(res_out,2,FUN=quantile,probs=c(0.025))
res_high1=apply(res_out,2,FUN=quantile,probs=c(0.975))

## Simulate new data and get residuals for simulated data

ysim1=simulate(glmb.out1,nsim=1,seed=NULL,pred=pred1,family="binomial",
               prior.weights=weights(glmb.out1))


res_ysim_out1=residuals(glmb.out1,ysim=ysim1)
res_low=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.025))
res_high=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.975))

# Plot Credible Interval bounds for Deviance Residuals

plot(res_mean~Age,ylim=c(-2.5,2.5),
main="Credible Interval Bound for Menarche - Logit Model Deviance Residuals",
xlab = "Age", ylab = "Avg. Dev. Res")
lines(Age, 0*res_mean,lty=1)
lines(Age, res_low,lty=1)
lines(Age, res_high,lty=1)
lines(Age, res_low1,lty=2)
lines(Age, res_high1,lty=2)




cleanEx()
nameEx("summary.glmb")
### * summary.glmb

flush(stderr()); flush(stdout())

### Name: summary.glmb
### Title: Summarizing Bayesian Generalized Linear Model Fits
### Aliases: summary.glmb print.summary.glmb

### ** Examples

###########################  Example for glmb function:
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

mysd<-1
mu<-matrix(0,5)
mu[1,1]=log(mean(counts))
V0<-((mysd)^2)*diag(5)
glmb.D93<-glmb(counts ~ outcome + treatment, family = poisson(),
pfamily=dNormal(mu=mu,Sigma=V0))
summary(glmb.D93)

###########################  Example for lmb function

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

## May move this step inside the Prior_Check function
lm.D9 <- lm(weight ~ group,x=TRUE,y=TRUE)
disp_ML=sigma(lm.D9)^2
n_prior=2
shape=n_prior/2
rate= disp_ML*shape

lmb.D9=lmb(weight ~ group,dNormal_Gamma(mu,V/disp_ML,shape=shape,rate=rate))
summary(lmb.D9)



cleanEx()
nameEx("summary.rglmb")
### * summary.rglmb

flush(stderr()); flush(stdout())

### Name: summary.rglmb
### Title: Summarizing Bayesian Generalized Linear Model Distribution
###   Functions
### Aliases: summary.rglmb print.summary.rglmb

### ** Examples

data(menarche2)

summary(menarche2)
plot(Menarche/Total ~ Age, data=menarche2)

Age2=menarche2$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche2$Menarche/menarche2$Total
wt=menarche2$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 

out<-rglmb(n = 1000, y=y, x=x, pfamily=dNormal(mu=mu,Sigma=V1), weights = wt, 
           family = binomial(logit)) 
summary(out)

print(summary(out),digits=4)

mean(out$iters)



cleanEx()
nameEx("vcov.glmb")
### * vcov.glmb

flush(stderr()); flush(stdout())

### Name: vcov.glmb
### Title: Calculate Variance-Covariance Matrux for a Fitted Model Object
### Aliases: vcov.glmb

### ** Examples

## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))

## ----glmb vcov----------------------------------------------------------------
vcov(glmb.D93)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
