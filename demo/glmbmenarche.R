library(MASS)
data(menarche)

summary(menarche)
plot(Menarche/Total ~ Age, data=menarche)


glm.out1<-glm(cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit), data=menarche)
summary(glm.out1)

# Change age variable to be difference from 13 so that prior point estimate 
# of 0 (corresponding to 50% of population) having undergone menarche by age 13
# makes sense 

Age2=menarche$Age-13


n<-1000

# Prior Assume that 10% of population with menarche at age 10
# and 90% by age 16
mu<-matrix(0,nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(2)

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]<-10  


glm.out2<-glm(cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit), data=menarche)
summary(glm.out2)

Like_std=summary(glm.out2)$coefficients[,2]

Menarche_Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V1)),glm.out2$coefficients,Like_std)

# Note that maximum likelihood estimate for slope is consistent with prior
# but prior point estimate for slope seems to be low compared to
# maximum likelihood estimate - keep model as is for now.

Menarche_Prior_Error_Checks
mu[2,1]=1.5
glm.out2$coefficients

# This still gets stuck somewhere - Need to do extensive QC on where it ends up stuck
Menarche_Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V1)),glm.out2$coefficients,Like_std)
glmb.out1<-glmb(n=n,cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),mu=mu,Sigma=V1,Gridtype=1, data=menarche)

summary(glmb.out2)
summary(glmb.out1)

# Now this hangs...
mu<-matrix(0,nrow=2,ncol=1)
mu[2,1]=0

V1<-1*diag(2)


glm.out2<-glm(cbind(Menarche,Total-Menarche)~Age2,family=binomial(probit), data=menarche)

summary(glm.out2)

Like_std=summary(glm.out2)$coefficients[,2]
Menarche_Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V1)),glm.out2$coefficients,Like_std)


glmb.out2<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age2,family=binomial(probit),mu=mu,Sigma=V1,Gridtype=1,data=menarche)


summary(glm.out2)
summary(glmb.out2)

# This runs fine

mu<-matrix(0,nrow=2,ncol=1)
mu[2,1]=0

V1<-1*diag(2)


glm.out3 <-glm(cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(cloglog), data=menarche)

summary(glm.out3)

Like_std=summary(glm.out3)$coefficients[,2]
Menarche_Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V1)),glm.out3$coefficients,Like_std)

glmb.out3<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age2,family=binomial(cloglog),mu=mu,Sigma=V1,Gridtype=1,data=menarche)


summary(glm.out3)
summary(glmb.out3)


