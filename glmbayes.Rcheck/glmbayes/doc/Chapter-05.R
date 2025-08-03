## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup,results = "hide"---------------------------------------------------
library(glmbayes)

## ----lm and lmb,results = "hide"----------------------------------------------
# Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
##Classical model-call
lm.D9 <- lm(weight ~ group)

# Bayesian Prior Setup and Call-Prior Should Generally be Customized
p_info=Prior_Setup(weight~group)
mu1=p_info$mu    
Sigma1=p_info$Sigma

#lmb.D9 <- lmb(weight ~ group,pfamily=dNormal_Gamma(mu1,Sigma1,shape=4,rate=0.1))

## ----glm and glmb,results = "hide"--------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
##Classical model-call
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())

# Bayesian Prior Setup and Call-Prior Should Generally be Customized
p_info2=Prior_Setup(counts ~ outcome + treatment)
mu2=p_info2$mu    
Sigma2=p_info2$Sigma

glmb.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu2,Sigma=Sigma2))


## ----Calling dNormal,results = "hide"-----------------------------------------
dNormal(mu=mu2,Sigma=Sigma2,dispersion=NULL)

## ----Calling dNormal_Gamma,results = "hide"-----------------------------------

dNormal_Gamma(mu=mu1,Sigma=Sigma1,shape=4,rate=0.1)

## ----Calling dGamma,results = "hide"------------------------------------------
b=lm.D9$coefficients
dGamma(shape=4,rate=0.1,beta=b)

## ----Calling Prior_Setup,results = "hide"-------------------------------------
Prior_Setup(weight ~ group)

## ----Calling Prior_Check,results = "hide"-------------------------------------
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())


Prior_Check(counts ~ outcome + treatment, family = poisson(),dNormal(mu2,Sigma2))

## ----Calling pfamily,results = "hide"-----------------------------------------
pfamily(glmb.D93)

## ----lm methods---------------------------------------------------------------
methods(class="lm")

## ----glm methods--------------------------------------------------------------
methods(class="glm")

## ----glmb methods-------------------------------------------------------------
methods(class="glmb")

## ----lmb methods--------------------------------------------------------------
methods(class="lmb")

## ----calling lmb,results = "hide"---------------------------------------------
## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
## Classical call
lm.D9 <- lm(weight ~ group,x=TRUE,y=TRUE)
lm_summary=summary(lm.D9)
## Bayesian Call
n<-10000
y<-lm.D9$y
#x<-as.matrix(lm.D9$x)

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

mu<-c(0,0)
mu=b_old  ### For testing purposes, set prior=b_old
P<-0.1*diag(2)
# Bayesian Prior Setup and Call-Prior Should Generally be Customized
p_info=Prior_Setup(weight~group)
x=p_info$x
outtemp4=rlmb(n=1000,y=y,x=x,pfamily=dNormal_Gamma(mu=mu,Sigma=solve(P),shape=shape,rate=rate))

#summary(outtemp4)


## ----calling rglmb------------------------------------------------------------
data(menarche2)
Age2=menarche2$Age-13

y=menarche2$Menarche/menarche2$Total
wt=menarche2$Total

## Use model.frame and model.matrix to derive x
#mf=model.frame(formula)
#x=model.matrix(formula,mf)

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

## Modify Prior_Setup so it can take a model matrix as well as an model object
mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
V1<-1*diag(as.numeric(2.0))
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 

out<-rglmb(n = 1000, y=y, x=x, pfamily=dNormal(mu=mu,Sigma=V1), weights = wt, 
           family = binomial(logit)) 

summary(out)

## ----rglmb methods------------------------------------------------------------
methods(class="rglmb")

## ----summary.rglmb methods----------------------------------------------------
methods(class="summary.rglmb")

## ----rlmb methods-------------------------------------------------------------
methods(class="rlmb")

