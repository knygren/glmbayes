rm(list=ls())

data(menarche2)
summary(menarche2)
#plot(Menarche/Total ~ Age, data=menarche2)

#################################################### Logit Model ############################################

glm.out1<-glm(cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit), data=menarche2)
summary(glm.out1)

residuals(glm.out1)

# Change age variable to be difference from 13 so that prior point estimate 
# of 0 (corresponding to 50% of population) having undergone menarche by age 13
# makes sense 

Age2=menarche2$Age-13

n<-1000

# Prior below implies point estimate that 10% of population at age 10
# and 90% by age 16 has undergone Menarche

mu<-matrix(0,nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(2)

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]<-10  ## May want to modify this
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 


glm.out2<-glm(cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit), data=menarche2)
summary(glm.out2)
residuals(glm.out2)

glm.out2$model
get_all_vars(formula(glm.out2$formula),menarche2)  ## Variables used to define alpha is in here     
all.vars(formula(glm.out2))

## Including prior.weights did not seem to change mean residuals much
glmb.out1<-glmb(n=n,cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),mu=mu,Sigma=V1,Gridtype=1, data=menarche2)




summary(glmb.out1)
olddata=data.frame(Age=menarche2$Age,
                   Menarche=menarche2$Menarche,Total=menarche2$Total,Age2=Age2)

pred1=predict(glmb.out1,type="response")


glmb.sum.out=summary(glmb.out1)
glmb.res=glmb.sum.out$residuals

# ~1.261 candidates per iid sample 

Age_New <- seq(8, 20, 0.25)
Age2_New=Age_New-13
mod_Object=glmb.out1
obs_size=median(menarche2$Total) ## Counts for sim from Binomial

newdata=data.frame(Age=Age_New,Age2=Age2_New)

mod_Object$glm$model
get_all_vars(formula(mod_Object$glm$formula),olddata)  ## Variables used to define alpha is in here     

pred_menarche=predict(mod_Object,newdata=newdata,olddata=olddata,type="response")

colMeans(pred_menarche)

pred_y=matrix(0,nrow=n,ncol=length(Age_New))
for(i in 1:n){
  pred_y[i,1:length(Age_New)]=rbinom(length(Age_New),size=obs_size,prob=pred_menarche[i,1:length(Age_New)
])
}



pred_m=colMeans(predict(mod_Object,newdata=newdata,olddata=olddata,type="response"))
pred_y_m=colMeans(pred_y/obs_size)

quant1_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.025))
quant2_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.975))
quant1_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.025))
quant2_m_y=apply(pred_y/obs_size,2,FUN=quantile,probs=c(0.975))


plot(Menarche/Total ~ Age, data=menarche2,main="Percentage of girls Who have had their first period")

## Plot for Fit on response scale
#plot(c(9,18), c(0,1), type = "n", xlab = "Age", ylab = "prob")
#text(menarche2$Age, menarche2$Menarche/menarche2$Total,label=c("M") )
#text(menarche2$Age, menarche2$Menarche/menarche2$Total)
lines(Age_New, pred_m,lty=1)
lines(Age_New, quant1_m,lty=2)
lines(Age_New, quant2_m,lty=2)
lines(Age_New, quant1_m_y,lty=2)
lines(Age_New, quant2_m_y,lty=2)

## Plot for Deviance Residuals

wts=glmb.out1$glm$prior.weights

res_out=residuals(glmb.out1)
res_mean=colMeans(res_out)

## These are not the intervals we want
## 

res_low1=apply(res_out,2,FUN=quantile,probs=c(0.025))
res_high1=apply(res_out,2,FUN=quantile,probs=c(0.975))

# Plot without ylim first and then set to make graph symmetric
plot(res_mean~Age2,main="Deviance Residuals for Menarche Model",xlab = "Age", ylab = "Avg. Dev. Res")
plot(res_mean~Age2,ylim=c(-2.5,2.5),main="Deviance Residuals for Menarche Model",xlab = "Age", ylab = "Avg. Dev. Res")
lines(Age2, 0*res_mean,lty=1)
lines(Age2, res_low1,lty=2)
lines(Age2, res_high1,lty=2)


ysim_temp=simulate.glmb(glmb.out1,nsim=1,seed=NULL,pred=pred1,family="binomial",
              wt=weights(glmb.out1$glm))

#colMeans(ysim_temp)



ysim1=glmb_simulate(pred1,family="binomial",wt=weights(glmb.out1$glm))

## ysim should replace y in calculations and not pred
res_ysim_out1=residuals(glmb.out1,ysim=ysim1)

colMeans(res_ysim_out1)
res_mean
res_low=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.025))
res_high=apply(res_ysim_out1,2,FUN=quantile,probs=c(0.975))

plot(res_mean~Age2,main="Deviance Residuals for Menarche Model",xlab = "Age", ylab = "Avg. Dev. Res")
plot(res_mean~Age2,ylim=c(-2.5,2.5),main="Deviance Residuals for Menarche Model",xlab = "Age", ylab = "Avg. Dev. Res")
lines(Age2, 0*res_mean,lty=1)
lines(Age2, res_low,lty=1)
lines(Age2, res_high,lty=1)


# Credible Interval bounds for Deviance Residuals

plot(res_mean~Age2,ylim=c(-2.5,2.5),main="Credible Interval Bound for Menarche - Logit Model Deviance Residuals",xlab = "Age", ylab = "Avg. Dev. Res")
lines(Age2, res_low,lty=1)
lines(Age2, res_high,lty=1)
lines(Age2, res_low1,lty=2)
lines(Age2, res_high1,lty=2)


ysim1[1:10,]

## These should be similar 
colMeans(ysim1)
colMeans(pred1)


wt=weights(glmb.out1$glm)
colMeans(ysim1/wt)

nvars=ncol(pred1)
nsims=nrow(pred1)
y_temp<-matrix(0,nrow=nrow(pred1),ncol=ncol(pred1))
wt=weights(glmb.out1$glm)
pred=pred1

i=1
 y_temp[i,1:nvars]=rbinom(n=nvars,size=round(wt),prob=pred[i,1:nvars])              





pred1[1,1:25]

ysim[1:10,]

methods(class="glmb")
methods(class="glm")

## Add this function - trivial but makes it more consistent
glmb.out1$glm$prior.weights
weights(glmb.out1$glm)



#################################################### Probit Model ############################################

mu<-matrix(0,nrow=2,ncol=1)
mu[2,1]=2*(0.9-0.5)/3 ## Implied Slope at age 13 should be 0.5*mu[1,1]*(0.9-0.5)/3 
V1<-1*diag(2)


# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]=((0.9-0.5)/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 

  
#Menarche_Prior_Error_Checks
glm.out2<-glm(cbind(Menarche,Total-Menarche)~Age2,family=binomial(probit), data=menarche2)
summary(glm.out2)

#  Note: Slope coefficient here is a bit more than half of the value of that in the logit model 

#Like_std=summary(glm.out2)$coefficients[,2]
#Menarche_Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V1)),glm.out2$coefficients,Like_std)

glmb.out2<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age2,family=binomial(probit),mu=mu,Sigma=V1,Gridtype=1,data=menarche2)
summary(glmb.out2)
olddata=data.frame(Age=menarche2$Age,Menarche=menarche2$Menarche,Total=menarche2$Total,Age2=Age2)
# ~1.251 candidates per iid sample


mod_Object=glmb.out2

pred_menarche=predict(mod_Object,newdata=newdata,olddata=olddata,type="response")
pred_m=colMeans(predict(mod_Object,newdata=newdata,olddata=olddata,type="response"))
quant1_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.025))
quant2_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.975))

## Plot for Fit on response scale
plot(c(8,20), c(0,1), type = "n", xlab = "Age", ylab = "prob")
text(menarche2$Age, menarche2$Menarche/menarche2$Total)
#text(menarche2$Age, menarche2$Menarche/menarche2$Total,label=c("M") )
lines(Age_New, pred_m,lty=1)
lines(Age_New, quant1_m,lty=2)
lines(Age_New, quant2_m,lty=2)

## Plot for Deviance Residuals

glmb.sum.out=summary(mod_Object)
glmb.res=glmb.sum.out$residuals

plot(c(9,18), c(-0.5,0.5), type = "n", xlab = "Age", ylab = "Avg. Dev. Res")
text(menarche2$Age, glmb.res,label=c("R") )
lines(Age_New, 0*pred_m,lty=1)


#################################################### Clog_Log Model ############################################

# Need to improve the prior specification here 

mu<-matrix(0,nrow=2,ncol=1)
mu[1,1]=log(-log(1-0.5))  # Derived prior
1-exp(-exp((mu[1,1])))    # implies point estimate at age 13 is 0.5
mu[2,1]=(log(-log(1-0.9))-mu[1,1])/3 # implies 90% menarche at age 16


V1<-1*diag(2)
V1[1,1]=((log(-log(1-0.9))-mu[1,1])/2)^2  ## Corresponding standard deviation
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 


glm.out3 <-glm(cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(cloglog), data=menarche2)
summary(glm.out3)

#1-exp(-exp(-0.59602))

#Like_std=summary(glm.out3)$coefficients[,2]
#Menarche_Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V1)),glm.out3$coefficients,Like_std)


glmb.out3<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age2,family=binomial(cloglog),mu=mu,Sigma=V1,Gridtype=1,data=menarche2)
summary(glmb.out3)

# Candidates per iid samples ~ 1.255



mod_Object=glmb.out3

pred_menarche=predict(mod_Object,newdata=newdata,olddata=olddata,type="response")
pred_m=colMeans(predict(mod_Object,newdata=newdata,olddata=olddata,type="response"))
quant1_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.025))
quant2_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.975))

## Plot for Fit on response scale
plot(c(8,20), c(0,1), type = "n", xlab = "Age", ylab = "prob")
text(menarche2$Age, menarche2$Menarche/menarche2$Total,label=c("M") )
lines(Age_New, pred_m,lty=1)
lines(Age_New, quant1_m,lty=2)
lines(Age_New, quant2_m,lty=2)

## Plot for Deviance Residuals

glmb.sum.out=summary(mod_Object)
glmb.res=glmb.sum.out$residuals

plot(c(9,18), c(-0.5,0.5), type = "n", xlab = "Age", ylab = "Avg. Dev. Res")
text(menarche2$Age, glmb.res,label=c("R") )
lines(Age_New, 0*pred_m,lty=1)


#################################################### Clog_Log Opposite Model ############################################

# Need to improve the prior specification here 

mu<-matrix(0,nrow=2,ncol=1)
mu[1,1]=log(-log(1-0.5))  # Derived prior
1-exp(-exp((mu[1,1])))    # implies point estimate at age 13 is 0.5
mu[2,1]=(log(-log(1-0.1))-mu[1,1])/3 # implies 90% menarche at age 16 [This was reversed!]


V1<-1*diag(2)
V1[1,1]=((log(-log(1-0.9))-mu[1,1])/2)^2  ## Corresponding standard deviation
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 


glm.out4 <-glm(cbind(Total-Menarche,Menarche) ~ Age2,family=binomial(cloglog), data=menarche2)
summary(glm.out4)

#1-exp(-exp(-0.59602))

#Like_std=summary(glm.out3)$coefficients[,2]
#Menarche_Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V1)),glm.out3$coefficients,Like_std)


glmb.out4<-glmb(n=n,cbind(Total-Menarche,Menarche)~Age2,family=binomial(cloglog),mu=mu,Sigma=V1,Gridtype=1,data=menarche2)
summary(glmb.out4)

# Candidates per iid samples ~ 1.255



mod_Object=glmb.out4

pred_menarche=predict(mod_Object,newdata=newdata,olddata=olddata,type="response")
pred_m=colMeans(predict(mod_Object,newdata=newdata,olddata=olddata,type="response"))
quant1_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.025))
quant2_m=apply(pred_menarche,2,FUN=quantile,probs=c(0.975))

## Plot for Fit on response scale
plot(c(8,20), c(0,1), type = "n", xlab = "Age", ylab = "prob")
text(menarche2$Age, menarche2$Menarche/menarche2$Total,label=c("M") )
lines(Age_New, 1-pred_m,lty=1)
lines(Age_New, 1-quant1_m,lty=2)
lines(Age_New, 1-quant2_m,lty=2)

## Plot for Deviance Residuals

glmb.sum.out=summary(mod_Object)
glmb.res=glmb.sum.out$residuals

## Reverse the residuals
plot(c(9,18), c(-0.5,0.5), type = "n", xlab = "Age", ylab = "Avg. Dev. Res")
text(menarche2$Age, -glmb.res,label=c("R") )
lines(Age_New, 0*pred_m,lty=1)

DIC_Info=rbind(
extractAIC(glmb.out1),
extractAIC(glmb.out2),
extractAIC(glmb.out3),
extractAIC(glmb.out4))

rownames(DIC_Info)=c("Logit","Probit","Clog-Log","Reverse Clog-Log")

AIC_Info=rbind(
  extractAIC(glm.out1),
  extractAIC(glm.out2),
  extractAIC(glm.out3),
  extractAIC(glm.out4))
