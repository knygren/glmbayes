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

Prior_Check(lm.D9,mu,Sigma_prior)

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

Prior_Check(lm.D9,mu,Sigma_prior)


prior=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion)
glmb.D9=glmb(n=1000,weight~group, family=gaussian(),prior=prior)
post_mode=glmb.D9$coef.mode

sum_out1=summary(glmb.D9)

## Try mean one standard deviation away to see how it works
#mu[1,1]=mu[1,1]-2*sum_out1$coefficients[1,3]
#mu[2,1]=mu[2,1]+2*sum_out1$coefficients[1,3]

prior=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion)
glmb.D9_v2=glmb(n=1000,weight~group, family=gaussian(),prior=prior)

summary(glmb.D9)
summary(glmb.D9_v2)


#lm_out=lm(y ~ x-1) # returns same model
RSS=sum(residuals(lm.D9)^2)

n_prior=4
n_data=length(y)
shape=(n_prior/2)
rate=n_prior*RSS/n_data

# For now, pass Precision as well as it is needed by the rglmb summary function
prior=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion,
shape=shape,rate=rate,Precision=solve(Sigma_prior))

## Using new function

#prior ## Currently standard normal so can't be cause of problem

## Try prior mean moved by 1 posterior sd from maximum likelihood estimate
## Took a couple of minutes to run 1000 iterations
## Seemed to simulate too close to the prior mean and not adequately penalize
## for being away from maximum likelihood estimate

#> t(glmb.D9$coef.means) Column means for fixed dispersion (from maximum likelihood)
#(Intercept)   groupTrt
#[1,]     4.7923 -0.1084

#> colMeans(beta_out) Column means for two-block Gibbs
#[1]  4.73585593 -0.07061893

#attr(prior,"Prior Type")=c("Normal")
#if(attr(prior,"Prior Type")=="Normal") print("Prior Type was Normal")
#attr(prior,"Prior Type")=c("Gamma")
#if(attr(prior,"Prior Type")=="Normal") print("Prior Type was Normal")
#if(attr(prior,"Prior Type")=="Gamma") print("Prior Type was Gamma")
#attr(prior,"Prior Type")=c("Normal-Gamma")
#if(attr(prior,"Prior Type")=="Normal") print("Prior Type was Normal")
#if(attr(prior,"Prior Type")=="Gamma") print("Prior Type was Gamma")
#if(attr(prior,"Prior Type")=="Normal-Gamma") print("Prior Type was Normal-Gamma")
#attr(prior,"Prior Type")=c("Independent-Normal-Gamma")
#print(paste("Prior Type was:",attr(prior,"Prior Type")))

sim1=rindep_norm_gamma_reg(n=1000,y,x,prior,offset2=NULL,wt=1)

summary(sim1)

## Note, last ratio is a constant for test 1 determined by distance between
## prior and maximum likelihood estimate
## for priors one sd away, it was -0.885 ## Note, much, much lower when approximate
## posterior mode used
## results seem to match Block Gibbs quite closely now that estimate 
## for joint posterior mode used to center candidates for beta
## 


qc1=cbind(sim1$dispersion,1/sim1$dispersion,sim1$test_out,
log(sim1$test_out[,2])*sim1$dispersion)

#qc1[1:30,]
#sim1$test_out[1:30,]


## Now, replicate using two-block Gibbs sampling

disp_out<-matrix(0,nrow=10000,ncol=1)
beta_out<-matrix(0,nrow=10000,ncol=2)

## Sample for beta given dispersion

## Initialize dispersion

dispersion2=dispersion

# Loop through two-block-Gibbs sampler

#### Burn-in iterations

for(i in 1:1000){
  prior2=list(mu=mu,Sigma=Sigma_prior, dispersion=dispersion2)
  glmb_out1=glmb(n=1,y~x-1,family=gaussian(),prior=prior2)
  
  b_old=glmb_out1$coefficients[1,]
#  beta_out[i,1:2]=glmb_out1$coefficients[1,1:2]
  
  ## sample for dispersion given beta [how did this run without passing these?]

  #shape=(n_prior/2)
  #rate=n_prior*RSS/n_data
  prior3=list(shape=shape, rate=rate,beta=b_old)
  
  disp_out1<-rglmb_dispersion(n=1,y,x,prior_list=prior3,
  offset= rep(0, length(y)),family=gaussian())
  dispersion2=disp_out1$dispersion
  
#  disp_out[i,1]=disp_out1$dispersion
}


## Post burn-in run
for(i in 1:10000){

prior2=list(mu=mu,Sigma=Sigma_prior, dispersion=dispersion2)
glmb_out1=glmb(n=1,y~x-1,family=gaussian(),prior=prior2)

b_old=glmb_out1$coefficients[1,]
beta_out[i,1:2]=glmb_out1$coefficients[1,1:2]

## sample for dispersion given beta
prior3=list(shape=shape, rate=rate,beta=b_old)

disp_out1<-rglmb_dispersion(n=1,y,x,prior_list=prior3,
offset= rep(0, length(y)),family=gaussian())

dispersion2=disp_out1$dispersion

disp_out[i,1]=disp_out1$dispersion
}

## Look at dispersion

## Maximum likelihood
dispersion

# New sampler (when prior=maximum likelihood estimate seemed to match Two-Block Gibbs)

mean(sim1$dispersion)
var(sim1$dispersion)

## Two-block Gibbs

mean(disp_out)
var(disp_out)

## Look at data precision

# maximum likelihood

1/dispersion

# new sampler

mean(1/sim1$dispersion)
var(1/sim1$dispersion)

##  Two-block Gibbs

mean(1/disp_out)
var(1/disp_out)


## Look at regression parameters

## maximum likelihood

lm.D9$coefficients

## Posterior mode and mean at maximum likelihood precision

t(glmb.D9$coef.mode)
t(glmb.D9$coef.means)

##  mean for new simulation [could be a bit too close to the assumed posterior mode]
## Estimates seem to match those from sampler with fixed dispersion
## instead of those with variable dispersion from two-block Givvs

colMeans(sim1$coefficients)

## mean for two-block Gibbs
colMeans(beta_out)

## Look at covariance
## Variance with standardized prior for intercepts seems
## slighly higher than in the model from two-block Gibbs
## This might be because iterations for two-block Gibbs were reduced 
## and variance should be bounded from above by true density

cov(glmb.D9$coefficients)  
cov(sim1$coefficients)  
cov(beta_out)

## Look at implied precision -- Results closer to two-block than
## to those from fixed dispersion - could difference be due to chance?

solve(cov(glmb.D9$coefficients) ) 
solve(cov(sim1$coefficients)  )
solve(cov(beta_out))

cor(beta_out[,1],disp_out[,1])  # -0.4168 
cor(beta_out[,2],disp_out[,1])  # 0.2615

sd1=sqrt(diag(cov(glmb.D9$coefficients)))
sd2=sqrt(diag(cov(sim1$coefficients)))
sd3=sqrt(diag(cov(beta_out)))

## Standard deviation is a tiny bit too low...
## Suggests precision for data too high (or dispersion too low)

rbind(colMeans(sim1$coefficients),sd2)
rbind(colMeans(beta_out),sd3)

## t-test suggest difference not due to chance
## Try iterating through to find better estimate for posterior mode
## similar to Block-Gibbs and might yield more consistent results
## Not sure why algorithm might require betastar to be at posterior mode
## to yield correct results

t.test(sim1$coefficients[,1],beta_out[,1])
t.test(sim1$coefficients[,2],beta_out[,2])

#####################################################################3

mean(sim1$test_out)
mean(sim1$iters)  
1/mean(sim1$iters)  ## 4.9% acceptance rate for current example - prior not too far

