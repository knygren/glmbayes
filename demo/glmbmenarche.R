library(MASS)
data(menarche)

summary(menarche)
plot(Menarche/Total ~ Age, data=menarche)

#################################################### Logit Model ############################################

glm.out1<-glm(cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit), data=menarche)
summary(glm.out1)

# Change age variable to be difference from 13 so that prior point estimate 
# of 0 (corresponding to 50% of population) having undergone menarche by age 13
# makes sense 

Age2=menarche$Age-13


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


glm.out2<-glm(cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit), data=menarche)
summary(glm.out2)

## Compare prior to Maximum likelihood estimates 

Like_std=summary(glm.out2)$coefficients[,2]
Menarche_Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V1)),glm.out2$coefficients,Like_std)

# Note that maximum likelihood estimate for intercept seems consistent with prior
# but prior point estimate for slope seems to be low compared to
# maximum likelihood estimate. Modify slope a bit

Menarche_Prior_Error_Checks

## Weird Cbars are printed during this call (need to modify code-remove printing)
## May be because Gridtype =1  

glmb.out1<-glmb(n=n,cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(logit),mu=mu,Sigma=V1,Gridtype=1, data=menarche)
summary(glmb.out1)

# ~1.261 candidates per iid sample 


#################################################### Probit Model ############################################

mu<-matrix(0,nrow=2,ncol=1)
mu[2,1]=2*(0.9-0.5)/3 ## Implied Slope at age 13 should be 0.5*mu[1,1]*(0.9-0.5)/3 
V1<-1*diag(2)


# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]=((0.9-0.5)/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 

  
Menarche_Prior_Error_Checks
glm.out2<-glm(cbind(Menarche,Total-Menarche)~Age2,family=binomial(probit), data=menarche)
summary(glm.out2)

#  Note: Slope coefficient here is a bit more than half of the value of that in the logit model 

Like_std=summary(glm.out2)$coefficients[,2]
Menarche_Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V1)),glm.out2$coefficients,Like_std)

glmb.out2<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age2,family=binomial(probit),mu=mu,Sigma=V1,Gridtype=1,data=menarche)
summary(glmb.out2)

# ~1.251 candidates per iid sample


#################################################### Clog_Log Model ############################################

# Need to improve the prior specification here 

mu<-matrix(0,nrow=2,ncol=1)
mu[1,1]=log(-log(1-0.5))  # Derived prior
1-exp(-exp((mu[1,1])))    # implies point estimate at age 13 is 0.5
mu[2,1]=(log(-log(1-0.9))-mu[1,1])/3 # implies 90% menarche at age 16


V1<-1*diag(2)
V1[1,1]=((log(-log(1-0.9))-mu[1,1])/2)^2  ## Corresponding standard deviation
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 


glm.out3 <-glm(cbind(Menarche, Total-Menarche) ~ Age2,family=binomial(cloglog), data=menarche)
summary(glm.out3)

#1-exp(-exp(-0.59602))

Like_std=summary(glm.out3)$coefficients[,2]
Menarche_Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V1)),glm.out3$coefficients,Like_std)


glmb.out3<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age2,family=binomial(cloglog),mu=mu,Sigma=V1,Gridtype=1,data=menarche)
summary(glmb.out3)

# Candidates per iid samples ~ 1.255

