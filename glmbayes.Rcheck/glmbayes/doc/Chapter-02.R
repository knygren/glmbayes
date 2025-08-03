## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup,echo = FALSE-------------------------------------------------------
library(glmbayes)

## ----menarche data------------------------------------------------------------
## Load menarche data
data(menarche2)
head(menarche2, 5)


## ----Analysis Setup-----------------------------------------------------------
## Number of variables in model
Age=menarche2$Age
nvars=2
set.seed(333)

## Reference Ages for setting of priors and Age_Difference
ref_age1=13  # user can modify this
ref_age2=15  ## user can modify this

## Define variables used later in analysis
Age2=Age-ref_age1
Age_Diff=ref_age2-ref_age1


## ----Prior Info---------------------------------------------------------------

## Point estimates at reference ages
m1=0.5  
m2=0.9

## Lower bound of prior credible intervals for point estimates
m1_lower=0.3
m2_lower=0.7

## Assumed correlation between the two (on link scale)
m_corr=0.4

## ----Logit: set up link function info and initialize prior matrices-----------

## Set up link function and initialize prior mean and Variance-Covariance matrices
bi_logit <- binomial(link="logit")
mu1<-matrix(0,nrow=nvars,ncol=1)
rownames(mu1)=c("Intercept","Age2")
colnames(mu1)=c("Prior Mean")
V1<-1*diag(nvars)
rownames(V1)=c("Intercept","Age2")
colnames(V1)=c("Intercept","Age2")

## ----Logit:set prior means----------------------------------------------------
## Prior mean for intercept is set to point estimate 
## at reference age1 (on logit scale)
mu1[1,1]=bi_logit$linkfun(m1)

## Prior mean for slope is set to difference in point estimates
## on logit scale divided by Age_Diff

mu1[2,1]=(bi_logit$linkfun(m2) -bi_logit$linkfun(m1))/Age_Diff 
print(mu1)



## ----Logit:set prior Variance Covariance matrix-------------------------------
## Implied standard deviations for point estimates on logit scale

sd_m1= (bi_logit$linkfun(m1) -bi_logit$linkfun(m1_lower))/1.96
sd_m2= (bi_logit$linkfun(m2) -bi_logit$linkfun(m2_lower))/1.96

## Also compute implied estimate for upper bound of confidence intervals

m1_upper=bi_logit$linkinv(bi_logit$linkfun(m1)+sd_m1*1.96)
m2_upper=bi_logit$linkinv(bi_logit$linkfun(m2)+sd_m2*1.96)
print("m1_upper is:")
m1_upper
print("m2_upper is:")
m2_upper


## Implied Standard deviation for slope (using variance formula for difference between two variables)
a=(1/Age_Diff)
sd_slope=sqrt((a*sd_m1)^2+(a*sd_m2)^2-2*a*a*(sd_m1*sd_m2*m_corr))

#Cov(m1,slope)=cov(m1, a*(m2-m1)) =a*E[(m1-E[m1])((m2-m1)-E[m2-m1])]
#   =a*E[(m1-E[m1])(m2-E[m2])]- a* E[(m1-E[m1])(m1-E[m1])]
##   =a*Cov[m1,m2] - a*Var[m1]
##  =a*sd_m1*sd_m2*m_corr-a* sd_m1*sd_m1
cov_V1=a*sd_m1*sd_m2*m_corr-a* sd_m1*sd_m1

# Set covariance matrix
V1[1,1]=sd_m1^2
V1[2,2]=sd_slope^2
V1[1,2]=cov_V1
V1[2,1]=V1[1,2]
print("V1 is:")
print(V1)

## ----Run Logit,results = "hide"-----------------------------------------------
Menarche_Model_Data=data.frame(Age=menarche2$Age,Total=as.integer(menarche2$Total),
Menarche=as.integer(menarche2$Menarche),Age2)

glmb.out1<-glmb(n=1000,cbind(Menarche, as.integer(Total-Menarche)) ~Age2,family=binomial(logit),
pfamily=dNormal(mu=mu1,Sigma=V1),data=Menarche_Model_Data)


## ----Print Logit--------------------------------------------------------------

# Print model output
print(glmb.out1)

# Print prior mean as comparison
print(t(mu1))


## ----Summary Logit------------------------------------------------------------
summary(glmb.out1)

