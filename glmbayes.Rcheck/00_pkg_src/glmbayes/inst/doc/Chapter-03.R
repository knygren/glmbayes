## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup,echo = FALSE-------------------------------------------------------
library(glmbayes)

## ----menarche data,results = "hide",echo = FALSE------------------------------
## Load menarche data
data(menarche2)
head(menarche2, 5)


## ----Analysis Setup,results = "hide",echo = FALSE-----------------------------
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


## ----Prior Info,results = "hide",echo = FALSE---------------------------------

## Point estimates at reference ages
m1=0.5  
m2=0.9

## Lower bound of prior credible intervals for point estimates
m1_lower=0.3
m2_lower=0.7

## Assumed correlation between the two (on link scale)
m_corr=0.4

## ----Logit: set up link function info and initialize prior matrices,results = "hide",echo = FALSE----

## Set up link function and initialize prior mean and Variance-Covariance matrices
bi_logit <- binomial(link="logit")
mu1<-matrix(0,nrow=nvars,ncol=1)
rownames(mu1)=c("Intercept","Age2")
colnames(mu1)=c("Prior Mean")
V1<-1*diag(nvars)
rownames(V1)=c("Intercept","Age2")
colnames(V1)=c("Intercept","Age2")

## ----Logit:set prior means,results = "hide",echo = FALSE----------------------
## Prior mean for intercept is set to point estimate 
## at reference age1 (on logit scale)
mu1[1,1]=bi_logit$linkfun(m1)

## Prior mean for slope is set to difference in point estimates
## on logit scale divided by Age_Diff

mu1[2,1]=(bi_logit$linkfun(m2) -bi_logit$linkfun(m1))/Age_Diff 
print(mu1)


## ----Logit:set prior Variance Covariance matrix,results = "hide",echo = FALSE----
## Implied standard deviations for point estimates on logit scale

sd_m1= (bi_logit$linkfun(m1) -bi_logit$linkfun(m1_lower))/1.96
sd_m2= (bi_logit$linkfun(m2) -bi_logit$linkfun(m2_lower))/1.96

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
print(V1)

## ----Run Logit,results = "hide",echo = FALSE----------------------------------
Menarche_Model_Data=data.frame(Age=menarche2$Age,Total=menarche2$Total,Menarche=menarche2$Menarche,Age2)
prior1=list(mu=mu1,Sigma=V1)
#glmb.out1<-glmb(n=10000,cbind(Menarche, Total-Menarche) ~ #Age2,family=binomial(logit),mu=mu1,Sigma=V1,data=Menarche_Model_Data)
glmb.out1<-glmb(n=10000,cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),pfamily=dNormal(mu=mu1,Sigma=V1),data=Menarche_Model_Data)


## ----Print Logit,results = "hide",echo = FALSE--------------------------------

# Print model output
print(glmb.out1)

# Print prior mean as comparison
print(t(mu1))


## ----Summary Logit,results = "hide",echo = FALSE------------------------------
summary(glmb.out1)

## ----Dobson Poisson Data------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

## ----Dobson Check_Prior-------------------------------------------------------
#glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(),x=TRUE)


Setup.D93=Prior_Setup(counts ~ outcome + treatment)

## ----Dobson Print Check_Prior-------------------------------------------------
print(Setup.D93)

## ----Venables Check_Prior-----------------------------------------------------
## example from Venables and Ripley (2002, pp. 190-2.)
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive = 20-numdead)

Setup.budworm=Prior_Setup(SF ~ sex*ldose)

## ----Venables Print Check_Prior-----------------------------------------------
print(Setup.budworm)

## ----Check Old Prior----------------------------------------------------------
Prior_Check(cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),pfamily=dNormal(mu1,V1),data=Menarche_Model_Data)



## ----Check Old Age------------------------------------------------------------
Prior_Check(cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit),pfamily=dNormal(mu1,V1)
,data=Menarche_Model_Data)



## ----Run wrong Age,results = "hide"-------------------------------------------

glmb.out2<-glmb(n=10000,cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit),pfamily=dNormal(mu=mu1,Sigma=V1),data=Menarche_Model_Data)

## ----Summary Logit - Wrong Age------------------------------------------------
summary(glmb.out2)

## ----Check 0 mean-------------------------------------------------------------
mu2=mu1
mu2[2,1]=0

pc=Prior_Check(cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),pfamily=dNormal(mu2,V1)
,data=Menarche_Model_Data)

print(pc)


## ----Run wrong mean,results = "hide"------------------------------------------
glmb.out3<-glmb(n=10000,cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),pfamily=dNormal(mu=mu2,Sigma=V1),data=Menarche_Model_Data)

## ----Summary Logit - 0 Slope--------------------------------------------------
summary(glmb.out3)

## ----Run wrong mean-adjusted Variance,results = "hide"------------------------
V2=V1
V2[2,2]=(pc[2,1]*sd_slope)^2
#glmb.out4<-glmb(n=10000,cbind(Menarche, Total-Menarche) ~ #Age2,family=binomial(logit),mu=mu2,Sigma=V2,data=Menarche_Model_Data)
glmb.out4<-glmb(n=10000,cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),pfamily=dNormal(mu=mu2,Sigma=V2),data=Menarche_Model_Data)

## ----Summary Logit - 0 Slope - Adjusted Variance------------------------------
summary(glmb.out4)

## ----DIC Comparison-----------------------------------------------------------
DIC_Out=rbind(extractAIC(glmb.out1),
extractAIC(glmb.out2),
extractAIC(glmb.out3),
extractAIC(glmb.out4))

colnames(DIC_Out)=c("pD","DIC")
rownames(DIC_Out)=c("Original Specification","Regular Age","Slope Mean=0 - No Var Adjustment",
                     "Slope Mean=0 - Var Adjustment")
print(DIC_Out)


