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



 ## Temporarily make prior dispersion very large 
#disp_temp=100*dispersion
prior_list=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion,
                 shape=shape,rate=rate,Precision=solve(Sigma_prior))

#Note: max_disp=0.9 seened to require just 8.552 candidates per acceptance
# if max_disp =1.5, this seems much higher and algorithm seems to possibly hang

set.seed(333)

 ptm <- proc.time()
 sim2=rindependent_norm_gamma_reg(n=1000,y,x,prior_list=prior_list,
offset=NULL,weights=1,max_disp=0.9)
 proc.time()-ptm

 
 
 
 # 100 iterations took: 
#  90.53 [RSS_Old2/n_obs]
# 33.20 [rate2/shape2]
# 19.51 [rate2/(shape2-1)]  # now 12.64 seconds # now 12.41 seconds # now 10.63 seconds
# 10.10 [dispstar=0.6452614]
 
 min(sim2$dispersion)
 max(sim2$dispersion)
 summary(sim2) 

 
hist(sim2$dispersion,50)
quantile(sim2$dispersion,probs=c(0.01,0.5,0.99))
mean(sim2$iters)  # Strength of Prior did not impact much on this 
1/mean(sim2$iters)

mean(sim2$dispersion)
quantile(sim2$dispersion,c(0.1,0.5,0.9))

lm_test2=lm(sim2$weight_out~sim2$dispersion)
lmc=lm_test2$coefficients

lmc[1]+lmc[2]*0.9
lmc[1]+lmc[2]*0.45

dispstar=0.61
lm_log2=lmc[2]*dispstar
lm_log1=lmc[1]+lm_log2-lm_log2*log(dispstar)

pred3=lm_log1+lm_log2*log(sim2$dispersion)

plot(sim2$weight_out~sim2$dispersion)

disp_new <- seq(0, 2, 0.1)
pred4=lm_log1+lm_log2*log(disp_new)
lines(disp_new,pred4,lty=2)

plot(sim2$weight_out-pred3)

lm_log2
disp_UB=0.9
lmc[1]+lmc[2]*disp_UB
lmc[1]+lmc[2]*disp_UB -(lm_log1+lm_log2*log(disp_UB))

disp_UB=2
lmc[1]+lmc[2]*disp_UB
lmc[1]+lmc[2]*disp_UB -(lm_log1+lm_log2*log(disp_UB))


lines(sim2$dispersion,lm_test2$fitted.values,lty=2)


summary(sim2$dispersion)

lm_test=lm(log(sim2$weight_out)~log(1/sim2$dispersion))
lm_test$coefficients

lines(log(1/sim2$dispersion),lm_test$fitted.values,lty=2)

fit_test=2.15-1.5*log(1/sim2$dispersion)

lines(log(1/sim2$dispersion),fit_test,lty=2)

cov(sim2$coefficients)  

mean(sim2$iters)  # Strength of Prior did not impact much on this 


disp_temp=rate/shape

glmb.D9_v3=glmb(n=1000,weight ~ group,family=gaussian(),dNormal_Gamma(mu,Sigma_prior/disp_temp,
            shape=shape,rate=rate))
summary(glmb.D9_v3)

mean(glmb.D9_v3$dispersion)

cov(sim2$coefficients)
cov(glmb.D9_v3$coefficients)


## Now, replicate using two-block Gibbs sampling


## Sample for beta given dispersion

## Initialize dispersion

dispersion2=dispersion

# Loop through two-block-Gibbs sampler

low=0.2845312
upp=1.126446
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
  
  a1=0
  while(a1==0){
    disp_out1<-rGamma_reg(n=1,y,x,prior_list=prior3,offset= rep(0, length(y)),
  weights=1,family=gaussian())
  if(disp_out1$dispersion>low & disp_out1$dispersion<upp) a1=1  
   
  }  
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
a1=0
while(a1==0){
  disp_out1<-rGamma_reg(n=1,y,x,prior_list=prior3,offset= rep(0, length(y)),
                        weights=1,family=gaussian())
  if(disp_out1$dispersion>low & disp_out1$dispersion<upp) a1=1  
  
}  
dispersion2=disp_out1$dispersion

disp_out[i,1]=disp_out1$dispersion
}



## Look at dispersion

## Maximum likelihood
dispersion

# New sampler (when prior=maximum likelihood estimate seemed to match Two-Block Gibbs)
## Dispersion is a bit too low - consistent with high precision needed to be penalized more

mean(sim2$dispersion)
var(sim2$dispersion)
mean(glmb.D9_v3$dispersion)
## Two-block Gibbs

mean(disp_out)
var(disp_out)

## Look at data precision [Grid currently seems to produce estimates that are too high]
## Likely need the penalty term for the higher RSS term in rejection method

# maximum likelihood

1/dispersion

# new sampler

mean(1/sim2$dispersion)
var(1/sim2$dispersion)

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

colMeans(sim2$coefficients)

## mean for two-block Gibbs
colMeans(beta_out)

## Look at covariance
## Variance with standardized prior for intercepts seems
## slighly higher than in the model from two-block Gibbs
## This might be because iterations for two-block Gibbs were reduced 
## and variance should be bounded from above by true density

cov(glmb.D9_v2$coefficients)  
cov(glmb.D9_v3$coefficients)  
cov(sim2$coefficients)  
#cov(sim3$coefficients)  
## With 10,000 iterations, this now seems much closer to the output from the function....
cov(beta_out)

## Look at implied precision -- Results closer to two-block than
## to those from fixed dispersion - could difference be due to chance?

solve(cov(glmb.D9_v2$coefficients) ) 
solve(cov(sim2$coefficients)  )
solve(cov(beta_out))


sd1=sqrt(diag(cov(glmb.D9$coefficients)))
sd3=sqrt(diag(cov(sim2$coefficients)))
sd4=sqrt(diag(cov(beta_out)))

sd3
sd4


rbind(colMeans(sim2$coefficients),sd3)
rbind(colMeans(beta_out),sd4)

## t-tests are a bit inconclusive [Really need a multivariate test here]
## These tests don't show much difference anymore.

t.test(sim2$coefficients[,1],beta_out[,1])
t.test(sim2$coefficients[,2],beta_out[,2])

