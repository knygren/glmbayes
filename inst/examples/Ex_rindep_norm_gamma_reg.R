
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



prior_list=list(mu=mu,Sigma=Sigma_prior,dispersion=dispersion,
                 shape=shape,rate=rate,Precision=solve(Sigma_prior))


set.seed(360)

 ptm <- proc.time()
 sim2=rindependent_norm_gamma_reg(n=1000,y,x,prior_list=prior_list,
offset=NULL,weights=1,max_disp_perc=0.99)
 proc.time()-ptm

 
 
 
 
 