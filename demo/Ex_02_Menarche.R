## Part 1: Code From Vignette 
## Chapter 2: Specifying Multivariate Normal Priors and Interpreting Model Outputs

## ----menarche data------------------------------------------------------------
## Load menarche data
data(menarche2)
head(menarche2, 5)

readline("press any key to continue")

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

readline("press any key to continue")

## ----Prior Info---------------------------------------------------------------

## Point estimates at reference ages
m1=0.5  
m2=0.9

## Lower bound of prior credible intervals for point estimates
m1_lower=0.3
m2_lower=0.7

## Assumed correlation between the two (on link scale)
m_corr=0.4

readline("press any key to continue")

## ----Logit: set up link function info and initialize prior matrices-----------

## Set up link function and initialize prior mean and Variance-Covariance matrices
bi_logit <- binomial(link="logit")
mu1<-matrix(0,nrow=nvars,ncol=1)
rownames(mu1)=c("Intercept","Age2")
colnames(mu1)=c("Prior Mean")
V1<-1*diag(nvars)
rownames(V1)=c("Intercept","Age2")
colnames(V1)=c("Intercept","Age2")

readline("press any key to continue")

## ----Logit:set prior means----------------------------------------------------
## Prior mean for intercept is set to point estimate 
## at reference age1 (on logit scale)
mu1[1,1]=bi_logit$linkfun(m1)

## Prior mean for slope is set to difference in point estimates
## on logit scale divided by Age_Diff

mu1[2,1]=(bi_logit$linkfun(m2) -bi_logit$linkfun(m1))/Age_Diff 
print(mu1)

readline("press any key to continue")

## ----Logit:set prior Variance Covariance matrix-------------------------------
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

readline("press any key to continue")

## ----Run Logit,results = "hide"-----------------------------------------------
Menarche_Model_Data=data.frame(Age=menarche2$Age,Total=menarche2$Total,
                               Menarche=menarche2$Menarche,Age2)

glmb.out1<-glmb(n=1000,cbind(Menarche, Total-Menarche) ~Age2,family=binomial(logit),
                pfamily=dNormal(mu=mu1,Sigma=V1),data=Menarche_Model_Data)

readline("press any key to continue")
## ----Print Logit--------------------------------------------------------------

# Print model output
print(glmb.out1)

# Print prior mean as comparison
print(t(mu1))

readline("press any key to continue")

## ----Summary Logit------------------------------------------------------------
summary(glmb.out1)

readline("press any key to continue")

## Part 2: Code From Vignette 
## Chapter 3: Using the Prior_Setup and Prior_Check Functions

## ----Check Prior From Part 1----------------------------------------------------------
Prior_Check(cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),
            pfamily=dNormal(mu1,V1),data=Menarche_Model_Data)


## ----Check Prior When Regular Age Used ------------------------------------------------------------
Prior_Check(cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit),
            pfamily=dNormal(mu1,V1),data=Menarche_Model_Data)

## ----Run Regular Age Model-------------------------------------------

glmb.out2<-glmb(n=1000,cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit),pfamily=dNormal(mu=mu1,Sigma=V1),data=Menarche_Model_Data)
summary(glmb.out2)

readline("press any key to continue")

## ----Check 0 mean Prior -------------------------------------------------------------

mu2=mu1
mu2[2,1]=0
pc=Prior_Check(cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),pfamily=dNormal(mu2,V1)
               ,data=Menarche_Model_Data)
print(pc)

readline("press any key to continue")

## ----Run With 0 Prior Mean ------------------------------------------
glmb.out3<-glmb(n=1000,cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),pfamily=dNormal(mu=mu2,Sigma=V1),data=Menarche_Model_Data)
summary(glmb.out3)

readline("press any key to continue")

## ----Adjust Variance and Run With "Wrong Mean"  ------------------------

V2=V1
V2[2,2]=(pc[2,1]*sd_slope)^2
glmb.out4<-glmb(n=1000,cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),pfamily=dNormal(mu=mu2,Sigma=V2),data=Menarche_Model_Data)
summary(glmb.out4)

readline("press any key to continue")

## ----DIC Comparison-----------------------------------------------------------
DIC_Info=rbind(extractAIC(glmb.out1),
               extractAIC(glmb.out2),
               extractAIC(glmb.out3),
               extractAIC(glmb.out4))

colnames(DIC_Info)=c("pD","DIC")
rownames(DIC_Info)=c("Original Specification","Regular Age","Slope Mean=0 - No Var Adjustment",
                     "Slope Mean=0 - Var Adjustment")
print(DIC_Info)

readline("press any key to continue")


## Predicting for original data
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
  pred_y[i,1:length(Age_New)]=rbinom(length(Age_New),size=obs_size,prob=pred_menarche[i,1:length(Age_New)
                                                                                      ])
}


pred_y_m=colMeans(pred_y/obs_size)

quant1_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.025))
quant2_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.975))
quant1_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.025))
quant2_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.975))


plot(Menarche/Total ~ Age, data=menarche2,main="Percentage of girls Who have had their first period")
lines(Age_New, pred_m,lty=1)
lines(Age_New, quant1_m,lty=2)
lines(Age_New, quant2_m,lty=2)
lines(Age_New, quant1_m_y,lty=2)
lines(Age_New, quant2_m_y,lty=2)

readline("press any key to continue")

## Using Deviance Residuals

## Get Original Residuals, their means, and credible bounds
res_out=residuals(glmb.out1)
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

plot(res_mean~Age,ylim=c(-2.5,2.5),main="Credible Interval Bound for Menarche - Logit Model Deviance Residuals",xlab = "Age", ylab = "Avg. Dev. Res")
lines(Age, 0*res_mean,lty=1)
lines(Age, res_low,lty=1)
lines(Age, res_high,lty=1)
lines(Age, res_low1,lty=2)
lines(Age, res_high1,lty=2)

readline("press any key to continue")

# Look a the raw data

menarche2

readline("press any key to continue")

##########################################   Probit Model  ########################


## ----Probit: set up link function info and initialize prior matrices-----------

## Set up link function and initialize prior mean and Variance-Covariance matrices
bi_probit <- binomial(link="probit")
mu2<-matrix(0,nrow=nvars,ncol=1)
rownames(mu2)=c("Intercept","Age2")
colnames(mu2)=c("Prior Mean")
V2<-1*diag(nvars)
rownames(V2)=c("Intercept","Age2")
colnames(V2)=c("Intercept","Age2")

readline("press any key to continue")

## ----Probit:set prior means----------------------------------------------------
## Prior mean for intercept is set to point estimate 
## at reference age1 (on logit scale)
mu2[1,1]=bi_probit$linkfun(m1)

## Prior mean for slope is set to difference in point estimates
## on logit scale divided by Age_Diff

mu2[2,1]=(bi_probit$linkfun(m2) -bi_probit$linkfun(m1))/Age_Diff 
print(mu2)

readline("press any key to continue")

## ----Probit:set prior Variance Covariance matrix-------------------------------
## Implied standard deviations for point estimates on logit scale

sd_m1= (bi_probit$linkfun(m1) -bi_probit$linkfun(m1_lower))/1.96
sd_m2= (bi_probit$linkfun(m2) -bi_probit$linkfun(m2_lower))/1.96

## Implied Standard deviation for slope (using variance formula for difference between two variables)
a=(1/Age_Diff)
sd_slope=sqrt((a*sd_m1)^2+(a*sd_m2)^2-2*a*a*(sd_m1*sd_m2*m_corr))

#Cov(m1,slope)=cov(m1, a*(m2-m1)) =a*E[(m1-E[m1])((m2-m1)-E[m2-m1])]
#   =a*E[(m1-E[m1])(m2-E[m2])]- a* E[(m1-E[m1])(m1-E[m1])]
##   =a*Cov[m1,m2] - a*Var[m1]
##  =a*sd_m1*sd_m2*m_corr-a* sd_m1*sd_m1
cov_V2=a*sd_m1*sd_m2*m_corr-a* sd_m1*sd_m1

# Set covariance matrix
V2[1,1]=sd_m1^2
V2[2,2]=sd_slope^2
V2[1,2]=cov_V2
V2[2,1]=V2[1,2]
print(V2)

readline("press any key to continue")

## ----Run Probit,results = "hide"-----------------------------------------------
Menarche_Model_Data=data.frame(Age=menarche2$Age,Total=menarche2$Total,
                               Menarche=menarche2$Menarche,Age2)

glmb.out2<-glmb(n=1000,cbind(Menarche, Total-Menarche) ~Age2,family=binomial(probit),
                pfamily=dNormal(mu=mu2,Sigma=V2),data=Menarche_Model_Data)

print(glmb.out2)
summary(glmb.out2)

readline("press any key to continue")


####################################    The Clog-Log Model - Specification 1   #####################


## ----cloglog: set up link function info and initialize prior matrices-----------

## Set up link function and initialize prior mean and Variance-Covariance matrices
bi_cloglog <- binomial(link="cloglog")
mu3<-matrix(0,nrow=nvars,ncol=1)
rownames(mu3)=c("Intercept","Age2")
colnames(mu3)=c("Prior Mean")
V3<-1*diag(nvars)
rownames(V3)=c("Intercept","Age2")
colnames(V3)=c("Intercept","Age2")

readline("press any key to continue")

## ----cloglog:set prior means----------------------------------------------------
## Prior mean for intercept is set to point estimate 
## at reference age1 (on logit scale)
mu3[1,1]=bi_cloglog$linkfun(m1)

## Prior mean for slope is set to difference in point estimates
## on logit scale divided by Age_Diff

mu3[2,1]=(bi_cloglog$linkfun(m2) -bi_cloglog$linkfun(m1))/Age_Diff 
print(mu3)

readline("press any key to continue")

## ----Probit:set prior Variance Covariance matrix-------------------------------
## Implied standard deviations for point estimates on logit scale

sd_m1= (bi_cloglog$linkfun(m1) -bi_cloglog$linkfun(m1_lower))/1.96
sd_m2= (bi_cloglog$linkfun(m2) -bi_cloglog$linkfun(m2_lower))/1.96

## Implied Standard deviation for slope (using variance formula for difference between two variables)
a=(1/Age_Diff)
sd_slope=sqrt((a*sd_m1)^2+(a*sd_m2)^2-2*a*a*(sd_m1*sd_m2*m_corr))

#Cov(m1,slope)=cov(m1, a*(m2-m1)) =a*E[(m1-E[m1])((m2-m1)-E[m2-m1])]
#   =a*E[(m1-E[m1])(m2-E[m2])]- a* E[(m1-E[m1])(m1-E[m1])]
##   =a*Cov[m1,m2] - a*Var[m1]
##  =a*sd_m1*sd_m2*m_corr-a* sd_m1*sd_m1
cov_V3=a*sd_m1*sd_m2*m_corr-a* sd_m1*sd_m1

# Set covariance matrix
V3[1,1]=sd_m1^2
V3[2,2]=sd_slope^2
V3[1,2]=cov_V3
V3[2,1]=V3[1,2]
print(V3)

readline("press any key to continue")

## ----Run cloglog,results = "hide"-----------------------------------------------
Menarche_Model_Data=data.frame(Age=menarche2$Age,Total=menarche2$Total,
                               Menarche=menarche2$Menarche,Age2)

glmb.out3<-glmb(n=1000,cbind(Menarche, Total-Menarche) ~Age2,family=binomial(cloglog),
                pfamily=dNormal(mu=mu3,Sigma=V3),data=Menarche_Model_Data)

print(glmb.out3)
summary(glmb.out3)

readline("press any key to continue")


####################################    The Clog-Log Model - Specification 2   #####################


## ----cloglog: set up link function info and initialize prior matrices-----------

## Set up link function and initialize prior mean and Variance-Covariance matrices
bi_cloglog <- binomial(link="cloglog")
mu4<-matrix(0,nrow=nvars,ncol=1)
rownames(mu4)=c("Intercept","Age2")
colnames(mu4)=c("Prior Mean")
V4<-1*diag(nvars)
rownames(V4)=c("Intercept","Age2")
colnames(V4)=c("Intercept","Age2")

readline("press any key to continue")

## ----cloglog:set prior means----------------------------------------------------
## Prior mean for intercept is set to point estimate 
## at reference age1 (on logit scale)
mu4[1,1]=bi_cloglog$linkfun((1-m1))

## Prior mean for slope is set to difference in point estimates
## on logit scale divided by Age_Diff

mu4[2,1]=(bi_cloglog$linkfun((1-m2))-bi_cloglog$linkfun((1-m1)))/Age_Diff 
print(mu4)

readline("press any key to continue")

## ----Probit:set prior Variance Covariance matrix-------------------------------
## Implied standard deviations for point estimates on logit scale

sd_m1= (bi_cloglog$linkfun((1-m1_lower))-bi_cloglog$linkfun((1-m1)))/1.96
sd_m2= (bi_cloglog$linkfun((1-m2_lower))-bi_cloglog$linkfun((1-m2)))/1.96

## Implied Standard deviation for slope (using variance formula for difference between two variables)
a=(1/Age_Diff)
sd_slope=sqrt((a*sd_m1)^2+(a*sd_m2)^2-2*a*a*(sd_m1*sd_m2*m_corr))

#Cov(m1,slope)=cov(m1, a*(m2-m1)) =a*E[(m1-E[m1])((m2-m1)-E[m2-m1])]
#   =a*E[(m1-E[m1])(m2-E[m2])]- a* E[(m1-E[m1])(m1-E[m1])]
##   =a*Cov[m1,m2] - a*Var[m1]
##  =a*sd_m1*sd_m2*m_corr-a* sd_m1*sd_m1
cov_V4=a*sd_m1*sd_m2*m_corr-a* sd_m1*sd_m1

# Set covariance matrix
V4[1,1]=sd_m1^2
V4[2,2]=sd_slope^2
V4[1,2]=cov_V4
V4[2,1]=V4[1,2]
print(V4)

readline("press any key to continue")

## ----Run cloglog,results = "hide"-----------------------------------------------
Menarche_Model_Data=data.frame(Age=menarche2$Age,Total=menarche2$Total,
                               Menarche=menarche2$Menarche,Age2)

glmb.out4<-glmb(n=1000,cbind(Total-Menarche,Menarche) ~Age2,family=binomial(cloglog),
                pfamily=dNormal(mu=mu4,Sigma=V4),data=Menarche_Model_Data)

print(glmb.out4)
summary(glmb.out4)

readline("press any key to continue")

DIC_Out=rbind(
  extractAIC(glmb.out1),
  extractAIC(glmb.out2),
  extractAIC(glmb.out3),
  extractAIC(glmb.out4))

rownames(DIC_Out)=c("Logit","Probit","Clog-Log","Reverse Clog-Log")


DIC_Out

readline("press any key to continue")

glm.out1=glm(cbind(Menarche, Total-Menarche) ~Age2,family=binomial(logit),data=Menarche_Model_Data)
glm.out2=glm(cbind(Menarche, Total-Menarche) ~Age2,family=binomial(probit),data=Menarche_Model_Data)
glm.out3=glm(cbind(Menarche, Total-Menarche) ~Age2,family=binomial(cloglog),data=Menarche_Model_Data)
glm.out4=glm(cbind(Total-Menarche,Menarche ) ~Age2,family=binomial(cloglog),data=Menarche_Model_Data)

AIC_Out=rbind(extractAIC(glm.out1),
              extractAIC(glm.out2),
              extractAIC(glm.out3),
              extractAIC(glm.out4))

rownames(AIC_Out)=c("Logit","Probit","Clog-Log","Reverse Clog-Log")
colnames(AIC_Out)=c("edf","AIC")
AIC_Out

readline("press any key to continue")


################################################################################


## Predict and Simulate for newdata

mod_Object=glmb.out2
pred1=predict(glmb.out2,type="response")
pred1_m=colMeans(pred1)
n=nrow(mod_Object$coefficients)
pred_menarche=predict(mod_Object,newdata=newdata,olddata=olddata,type="response")
pred_y=matrix(0,nrow=n,ncol=length(Age_New))
for(i in 1:n){
  pred_y[i,1:length(Age_New)]=rbinom(length(Age_New),size=obs_size,
                                     prob=pred_menarche[i,1:length(Age_New)])
}


## Produce mean and quantile information

pred_m=colMeans(pred_menarche)
pred_y_m=colMeans(pred_y/obs_size)
quant1_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.025))
quant2_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.975))
quant1_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.025))
quant2_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.975))


plot(Menarche/Total ~ Age, data=menarche2,main="Percentage of girls Who have had their first period")
lines(Age_New, pred_m,lty=1)
lines(Age_New, quant1_m,lty=2)
lines(Age_New, quant2_m,lty=2)
lines(Age_New, quant1_m_y,lty=2)
lines(Age_New, quant2_m_y,lty=2)

readline("press any key to continue")


## Get Residuals for model, their means, and credible bounds
res_out=residuals(glmb.out2)
res_mean=colMeans(res_out)
sqrt(mean(res_mean^2))
res_low1=apply(res_out,2,FUN=quantile,probs=c(0.025))
res_high1=apply(res_out,2,FUN=quantile,probs=c(0.975))

## Simulate new data and get residuals for simulated data

ysim1=simulate(glmb.out2,nsim=1,seed=NULL,pred=pred1,family="binomial",
               prior.weights=weights(glmb.out2))
res_ysim_out1=residuals(glmb.out2,ysim=ysim1)
res_low=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.025))
res_high=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.975))


# Plot Credible Interval bounds for Deviance Residuals

plot(res_mean~Age,ylim=c(-2.5,2.5),
     main="Credible Interval Bounds for Probit Model Deviance Residuals",
     xlab = "Age", ylab = "Dev. Residuals")
lines(Age, 0*res_mean,lty=1)
lines(Age, res_low,lty=1)
lines(Age, res_high,lty=1)
lines(Age, res_low1,lty=2)
lines(Age, res_high1,lty=2)

readline("press any key to continue")

### clog-log Model

## Predict and Simulate for newdata
mod_Object=glmb.out3
pred1=predict(glmb.out3,type="response")
pred1_m=colMeans(pred1)
n=nrow(mod_Object$coefficients)
pred_menarche=predict(mod_Object,newdata=newdata,olddata=olddata,type="response")
pred_y=matrix(0,nrow=n,ncol=length(Age_New))
for(i in 1:n){
  pred_y[i,1:length(Age_New)]=rbinom(length(Age_New),size=obs_size,
                                     prob=pred_menarche[i,1:length(Age_New)])
}


## Produce mean and quantile information

pred_m=colMeans(pred_menarche)
pred_y_m=colMeans(pred_y/obs_size)
quant1_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.025))
quant2_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.975))
quant1_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.025))
quant2_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.975))


plot(Menarche/Total ~ Age, data=menarche2,main="Percentage of girls Who have had their first period")
lines(Age_New, pred_m,lty=1)
lines(Age_New, quant1_m,lty=2)
lines(Age_New, quant2_m,lty=2)
lines(Age_New, quant1_m_y,lty=2)
lines(Age_New, quant2_m_y,lty=2)

readline("press any key to continue")


## Get Residuals for model, their means, and credible bounds
res_out=residuals(glmb.out3)
res_mean=colMeans(res_out)
sqrt(mean(res_mean^2))
res_low1=apply(res_out,2,FUN=quantile,probs=c(0.025))
res_high1=apply(res_out,2,FUN=quantile,probs=c(0.975))

## Simulate new data and get residuals for simulated data

ysim1=simulate(glmb.out3,nsim=1,seed=NULL,pred=pred1,family="binomial",
               prior.weights=weights(glmb.out3))
res_ysim_out1=residuals(glmb.out3,ysim=ysim1)
res_low=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.025))
res_high=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.975))

# Plot Credible Interval bounds for Deviance Residuals

plot(res_mean~Age,ylim=c(-4.0,4.0),
     main="Credible Interval Bounds for cloglog Model Deviance Residuals",
     xlab = "Age", ylab = "Dev. Residuals")
lines(Age, 0*res_mean,lty=1)
lines(Age, res_low,lty=1)
lines(Age, res_high,lty=1)
lines(Age, res_low1,lty=2)
lines(Age, res_high1,lty=2)


readline("press any key to continue")

### Reverse clog-log Model

## Predict and Simulate for newdata
mod_Object=glmb.out4
pred1=predict(glmb.out4,type="response")
pred1_m=colMeans(pred1)
n=nrow(mod_Object$coefficients)
pred_menarche=predict(mod_Object,newdata=newdata,olddata=olddata,type="response")
pred_y=matrix(0,nrow=n,ncol=length(Age_New))
for(i in 1:n){
  pred_y[i,1:length(Age_New)]=rbinom(length(Age_New),size=obs_size,
                                     prob=pred_menarche[i,1:length(Age_New)])
}


## Produce mean and quantile information

pred_m=colMeans(pred_menarche)
pred_y_m=colMeans(pred_y/obs_size)
quant1_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.025))
quant2_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.975))
quant1_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.025))
quant2_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.975))


plot(Menarche/Total ~ Age, data=menarche2,main="Percentage of girls Who have had their first period")
lines(Age_New, 1-pred_m,lty=1)
lines(Age_New, 1-quant1_m,lty=2)
lines(Age_New, 1-quant2_m,lty=2)
lines(Age_New, 1-quant1_m_y,lty=2)
lines(Age_New, 1-quant2_m_y,lty=2)

readline("press any key to continue")

## Get Residuals for model, their means, and credible bounds
res_out=residuals(glmb.out4)
res_mean=colMeans(res_out)
sqrt(mean(res_mean^2))
res_low1=apply(res_out,2,FUN=quantile,probs=c(0.025))
res_high1=apply(res_out,2,FUN=quantile,probs=c(0.975))

## Simulate new data and get residuals for simulated data

ysim1=simulate(glmb.out4,nsim=1,seed=NULL,pred=pred1,family="binomial",
               prior.weights=weights(glmb.out4))
res_ysim_out1=residuals(glmb.out4,ysim=ysim1)
res_low=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.025))
res_high=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.975))


# Plot Credible Interval bounds for Deviance Residuals

plot(-res_mean~Age,ylim=c(-3.0,3.0),
     main="Credible Interval Bounds for Reverse clog-log Model Deviance Residuals",
     xlab = "Age", ylab = "Dev. Residuals")
lines(Age, 0*res_mean,lty=1)
lines(Age, -res_low,lty=1)
lines(Age, -res_high,lty=1)
lines(Age, -res_low1,lty=2)
lines(Age, -res_high1,lty=2)



