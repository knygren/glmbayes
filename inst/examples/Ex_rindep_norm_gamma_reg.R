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
Sigma_prior=0.1*Sigma_prior

glmb.D9_v2=glmb(n=1000,weight~group, family=gaussian(),dNormal(mu=mu,Sigma=Sigma_prior,dispersion=dispersion))

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


# For now, pass Precision as well as it is needed by the rglmb summary function
prior_list=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion,
shape=shape,rate=rate,Precision=solve(Sigma_prior))

## Using new function

ptm <- proc.time()
sim1=rindependent_norm_gamma_reg(n=1000,y,x,prior_list=prior_list,offset=NULL,weights=1)
 proc.time()-ptm

 # Try a model where prior mean is equal to maximum likelihood estimate

#mu_temp=lm.D9$coefficients
mu_temp=mu

 ## Temporarily make prior dispersion very large 
#disp_temp=100*dispersion
prior_list2=list(mu=mu_temp,Sigma=Sigma_prior,dispersion=dispersion,
                 shape=shape,rate=rate,Precision=solve(Sigma_prior))
 ptm <- proc.time()
 sim2=rindependent_norm_gamma_reg_v2(n=1000,y,x,prior_list=prior_list2,offset=NULL,weights=1)
 proc.time()-ptm
 
summary(sim1)
summary(sim2)

mean(sim1$iters)
mean(sim2$iters)

mu=mu_temp

mean(sim1$dispersion) 
mean(1/sim1$dispersion)  # should be ~1.8 ## Estimate with only base terms --> 2.07 =(shape/rate)

disp_temp=rate/shape

glmb.D9_v3=glmb(n=1000,weight ~ group,family=gaussian(),dNormal_Gamma(mu,Sigma_prior/disp_temp,
            shape=shape,rate=rate))
summary(glmb.D9_v3)

cov(sim1$coefficients)
cov(glmb.D9_v3$coefficients)


## Now, replicate using two-block Gibbs sampling


## Sample for beta given dispersion

## Initialize dispersion

dispersion2=dispersion

# Loop through two-block-Gibbs sampler

#### Burn-in iterations

for(i in 1:1000){
  prior2=list(mu=mu,Sigma=Sigma_prior, dispersion=dispersion2)
  glmb_out1=glmb(n=1,y~x-1,family=gaussian(),dNormal(mu=mu,Sigma=Sigma_prior,dispersion=dispersion))
  
  b_old=glmb_out1$coefficients[1,]
#  beta_out[i,1:2]=glmb_out1$coefficients[1,1:2]
  
  ## sample for dispersion given beta [how did this run without passing these?]

  #shape=(n_prior/2)
  #rate=n_prior*RSS/n_data
  prior3=list(shape=shape, rate=rate,beta=b_old)
  
  #disp_out1<-rglmb_dispersion(n=1,y,x,prior_list=prior3,
  #offset= rep(0, length(y)),family=gaussian())
  
  disp_out1<-rGamma_reg(n=1,y,x,prior_list=prior3,offset= rep(0, length(y)),weights=1,family=gaussian())
  dispersion2=disp_out1$dispersion
  
#  disp_out[i,1]=disp_out1$dispersion
}


n_out=10000

disp_out<-matrix(0,nrow=n_out,ncol=1)
beta_out<-matrix(0,nrow=n_out,ncol=2)


## Post burn-in run
for(i in 1:n_out){

#prior2=list(mu=mu,Sigma=Sigma_prior, dispersion=dispersion2)
#glmb_out1=glmb(n=1,y~x-1,family=gaussian(),prior=prior2)
prior2=list(mu=mu,Sigma=Sigma_prior, dispersion=dispersion2)
glmb_out1=glmb(n=1,y~x-1,family=gaussian(),dNormal(mu=mu,Sigma=Sigma_prior,dispersion=dispersion))



b_old=glmb_out1$coefficients[1,]
beta_out[i,1:2]=glmb_out1$coefficients[1,1:2]

## sample for dispersion given beta
prior3=list(shape=shape, rate=rate,beta=b_old)

#disp_out1<-rglmb_dispersion(n=1,y,x,prior_list=prior3,
#offset= rep(0, length(y)),family=gaussian())
disp_out1<-rGamma_reg(n=1,y,x,prior_list=prior3,offset= rep(0, length(y)),weights=1,family=gaussian())
dispersion2=disp_out1$dispersion


dispersion2=disp_out1$dispersion

disp_out[i,1]=disp_out1$dispersion
}



## Look at dispersion

## Maximum likelihood
dispersion

# New sampler (when prior=maximum likelihood estimate seemed to match Two-Block Gibbs)

mean(sim1$dispersion)
var(sim1$dispersion)

#mean(sim2$dispersion)
#var(sim2$dispersion)

#mean(sim3$dispersion)
#var(sim3$dispersion)
mean(glmb.D9_v3$dispersion)
## Two-block Gibbs

mean(disp_out)
var(disp_out)

## Look at data precision [Grid currently seems to produce estimates that are too high]
## Likely need the penalty term for the higher RSS term in rejection method

# maximum likelihood

1/dispersion

# new sampler

mean(1/sim1$dispersion)
var(1/sim1$dispersion)


#mean(1/sim2$dispersion)
#var(1/sim2$dispersion)

#mean(1/sim3$dispersion)
#var(1/sim3$dispersion)

##  Two-block Gibbs

mean(1/disp_out)
var(1/disp_out)


## Look at regression parameters

## maximum likelihood

lm.D9$coefficients

## Posterior mode and mean at maximum likelihood precision

t(glmb.D9_v2$coef.mode)
t(glmb.D9_v2$coef.means)

##  mean for new simulation [could be a bit too close to the assumed posterior mode]
## Estimates seem to match those from sampler with fixed dispersion
## instead of those with variable dispersion from two-block Givvs

colMeans(sim1$coefficients)
#colMeans(sim2$coefficients)
#colMeans(sim3$coefficients)

## mean for two-block Gibbs
colMeans(beta_out)

## Look at covariance
## Variance with standardized prior for intercepts seems
## slighly higher than in the model from two-block Gibbs
## This might be because iterations for two-block Gibbs were reduced 
## and variance should be bounded from above by true density

cov(glmb.D9_v2$coefficients)  
cov(glmb.D9_v3$coefficients)  
cov(sim1$coefficients)  
#cov(sim2$coefficients)  
#cov(sim3$coefficients)  
## With 10,000 iterations, this now seems much closer to the output from the function....
cov(beta_out)

## Look at implied precision -- Results closer to two-block than
## to those from fixed dispersion - could difference be due to chance?

solve(cov(glmb.D9_v2$coefficients) ) 
solve(cov(sim1$coefficients)  )
solve(cov(beta_out))

cor(beta_out[,1],disp_out[,1])  # -0.4168 
cor(beta_out[,2],disp_out[,1])  # 0.2615

sd1=sqrt(diag(cov(glmb.D9$coefficients)))
sd2=sqrt(diag(cov(sim1$coefficients)))
sd3=sqrt(diag(cov(beta_out)))

sd2
sd3

## Standard deviation now slightly higher than two-block Gibbs

rbind(colMeans(sim1$coefficients),sd2)
rbind(colMeans(beta_out),sd3)

## t-tests are a bit inconclusive [Really need a multivariate test here]

t.test(sim1$coefficients[,1],beta_out[,1])
t.test(sim1$coefficients[,2],beta_out[,2])

#####################################################################3

mean(sim1$test_out)
mean(sim1$iters)  
1/mean(sim1$iters)  ## 4.9% acceptance rate for current example - prior not too far
