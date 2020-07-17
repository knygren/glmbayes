data(carinsca)
carinsca$Merit <- ordered(carinsca$Merit)
carinsca$Class <- factor(carinsca$Class)
options(contrasts=c("contr.treatment","contr.treatment"))
Claims=carinsca$Claims
Insured=carinsca$Insured
Merit=carinsca$Merit
Class=carinsca$Class
Cost=carinsca$Cost

scale<-0.1275 # SAS estimate
dispersion<-1/scale
out <- glm(Cost/Claims~Merit+Class,family=Gamma(link="log"),weights=Claims,x=TRUE)
summary(out)
disp=gamma.dispersion(out)
ps=Prior_Setup(Cost/Claims~Merit+Class)
mu=ps$mu
V=ps$Sigma
mu[1,1]=log(mean(Cost/Claims))
Prior_Check(Cost/Claims~Merit+Class,family=Gamma(link="log"),pfamily=dNormal(mu=mu,Sigma=V),weights=Claims/disp)

out3 <- glmb(Cost/Claims~Merit+Class,family=Gamma(link="log"),pfamily=dNormal(mu=mu,Sigma=V,dispersion=disp),weights=Claims)
summary(out,dispersion=gamma.dispersion(out))
summary(out3)

###########################################   Below part should not be part of final example ####################

#out2 <- glmb(Cost/Claims~Merit+Class,family=Gamma(link="log"),pfamily=dNormal(mu=mu,Sigma=V),weights=Claims/disp)
#
#summary(out2)










ps=Prior_Setup(Claims/Insured~Merit+Class)
mu=ps$mu
V=ps$Sigma
mu[1,1]=log(mean(Claims/Insured))

Prior_Check(Claims/Insured~Merit+Class,family="poisson",pfamily=dNormal(mu=mu,Sigma=V))

## The # of Insured is very large

#out1=glm(Claims/Insured~Merit+Class,family="poisson",weights=Insured,x=TRUE,y=TRUE)
#summary(out1)
#out1b=glm(Claims/Insured~Merit+Class,family="quasipoisson",weights=Insured,x=TRUE,y=TRUE)
#summary(out1b)

out2 <- glmb(Claims/Insured~Merit+Class,family="poisson",dNormal(mu=mu,Sigma=V),weights=Insured)
summary(out2)
out3 <- glmb(Claims/Insured~Merit+Class,family="quasipoisson",dNormal(mu=mu,Sigma=V),weights=Insured)
summary(out3)

mean(out3$dispersion)

# These ratios are nearly equal to the mean of out3$dispersion because the data overwhelms the prior in this model
# In the classical model, they would be identical
out2$Dbar/out3$Dbar  
out2$Dthetabar/out3$Dthetabar

out2$Dbar-out2$Dthetabar  # Effective number of parameters in original model
out3$Dbar-out3$Dthetabar  # Effective number of parameters in quasi-poisson model (in this case higher)

out2$Dbar+(out2$Dbar-out2$Dthetabar)  # DIC in original model
out3$Dbar+(out3$Dbar-out3$Dthetabar)  # DIC in quasi-poisson model without scaling (in this case lower)
(out3$Dbar+(out3$Dbar-out3$Dthetabar))*mean(out3$dispersion)

# Because the first term in the DIC calculation dominates for this model, the DIC
# For the first model is larger by nearly the estimate for the dispersion

out2$DIC/out3$DIC

out3$DIC*mean(out3$dispersion)


#extractAIC(out2)
#extractAIC(out3)

#######################################################################################################






mean(out3$dispersion)

## This should be formula for Expected Residual Deviance
mean(rowSums(residuals(out3)*residuals(out3)))
mean(rowSums(residuals(out2)*residuals(out2)))

x=out3$x
y=out3$y
nobs=length(y)
nvar=ncol(x)


help(rglmb)
xbeta=x%*%out2$coef.mode
vtemp=out3$dispersion


disp_ML=48
n_prior=1
shape=n_prior/2
rate= disp_ML*shape

n_sim=100
beta_temp=matrix(0,nrow=nobs)
b_out=matrix(0,nrow=n_sim,ncol=nobs)
beta_out=matrix(0,nrow=n_sim,ncol=nobs)

nu_out=matrix(0,nrow=n_sim,ncol=nvar)
v_out=rep(0,n_sim)


########################  Gibbs Sampler Burn-in

n_sim1=100

for(k in 1:n_sim1){
  
  ###  Update random effects
  
for(j in 1:nobs){
outtemp=rglmb(n = 1, y=y[j], x=as.matrix(1,nrow=1,ncol=1), 
      family=poisson, pfamily=dNormal(mu=xbeta[j],Sigma=vtemp), offset = NULL, weights = Insured[j])
beta_temp[j]=outtemp$coefficients


}
  
b_out[k,1:nobs]=b_temp
  
b_temp=t(beta_temp-xbeta)
beta_out[k,1:nobs]=t(beta_temp)

### Update Fixed Effects and Dispersion

outtemp2=rglmb(n = 1, y=beta_temp, x=x, 
               family=gaussian(), pfamily=dNormal_Gamma(mu=mu,Sigma=V/disp_ML,shape=shape,rate=rate), 
               offset = NULL, weights = 1)

xbeta=x%*%outtemp2$coefficients[1,1:nvar]
vtemp=outtemp2$dispersion

nu_out[k,1:nvar]=t(as.matrix(outtemp2$coefficients[1,1:nvar],nrow=1,ncol=nvar))

v_out[k]=vtemp

}

########  End of Burn-In iterations


########################  Gibbs Sampler Remaining Simulation


for(k in 1:n_sim){
  
  ###  Update random effects
  
  for(j in 1:nobs){
    outtemp=rglmb(n = 1, y=y[j], x=as.matrix(1,nrow=1,ncol=1), 
                  family=poisson, pfamily=dNormal(mu=xbeta[j],Sigma=vtemp), offset = NULL, weights = Insured[j])
    beta_temp[j]=outtemp$coefficients
    
    
  }
  
  b_out[k,1:nobs]=b_temp
  
  b_temp=t(beta_temp-xbeta)
  beta_out[k,1:nobs]=t(beta_temp)
  
  ### Update Fixed Effects and Dispersion
  
  outtemp2=rglmb(n = 1, y=beta_temp, x=x, 
                 family=gaussian(), pfamily=dNormal_Gamma(mu=mu,Sigma=V/disp_ML,shape=shape,rate=rate), 
                 offset = NULL, weights = 1)
  
  xbeta=x%*%outtemp2$coefficients[1,1:nvar]
  vtemp=outtemp2$dispersion
  
  nu_out[k,1:nvar]=t(as.matrix(outtemp2$coefficients[1,1:nvar],nrow=1,ncol=nvar))
  
  v_out[k]=vtemp
  
}

########  End of Burn-In iterations



nu_out




#colMeans(exp(beta_out))
#y

#outtemp2=rglmb(n = 1, y=beta_temp, x=x, 
#      family=gaussian(), pfamily=dNormal(mu=mu,Sigma=V), offset = NULL, weights = Insured)


outtemp2=rglmb(n = 1, y=beta_temp, x=x, 
               family=gaussian(), pfamily=dNormal_Gamma(mu=mu,Sigma=V/disp_ML,shape=shape,rate=rate), 
               offset = NULL, weights = Insured)

xbeta=x%*%outtemp2$coefficients[1,1:nvar]



outtemp2$dispersion

outtemp2=rglmb(n = 1, y=beta_temp, x=x, 
               family=gaussian(), pfamily=dIndependent_Normal_Gamma(mu=mu,Sigma=V,shape=shape,rate=rate), offset = NULL, weights = Insured)



outtemp2$coefficients




out2$coef.mode


plot(b_out[,1])

colMeans(b_out)
colMeans(beta_out)
colMeans(residuals(out2))


#t(beta_temp)
colMeans(residuals(out2))




rglmb(n = 1, y=y[2], x=as.matrix(1,nrow=1,ncol=1), 
      family="poisson", pfamily=dNormal(mu=xbeta[2],Sigma=vtemp), offset = NULL, weights = Insured[1])





Res_m=colMeans(residuals(out2))
RSS_m=sum(Res_m*Res_m)
sqrt(RSS_m)
2*sqrt(RSS_m)


### Demonstration of quasi-poisson Dispersion calculation
### Could add calculation so this is calculated on back end...
### Really should be just a diagnostic to see if the dispersion=1 seems incorrect

### https://stats.stackexchange.com/questions/62006/definition-of-dispersion-parameter-for-quasipoisson-family

m<-length(y)
k=length(out1$coefficients)
fitted(out1)
ytemp=Claims/Insured
fittemp=fitted(out1)
wt=Insured
1/(m-k)*sum((ytemp-fittemp)^2*wt/fitted(out1))
(1/(m-k))*sum((ytemp-fittemp)^2*wt/fitted(out1))
(1/(m-k))*sum((ytemp-fitted(out1))^2/fitted(out1))


###########################################


#carinsca$Merit <- ordered(carinsca$Merit)
#carinsca$Class <- factor(carinsca$Class)
#options(contrasts=c("contr.treatment","contr.treatment"))


scale<-0.1275 # SAS estimate

dispersion<-1/scale

out <- glm(Cost/Claims~Merit+Class,family=Gamma(link="log"),weights=Claims,x=TRUE)
summary(out)

disp=gamma.dispersion(out)

ps=Prior_Setup(Cost/Claims~Merit+Class)
mu=ps$mu
V=ps$Sigma

mu[1,1]=log(mean(Cost/Claims))

Prior_Check(Cost/Claims~Merit+Class,family=Gamma(link="log"),pfamily=dNormal(mu=mu,Sigma=V),weights=Claims/disp)



out2 <- glmb(Cost/Claims~Merit+Class,family=Gamma(link="log"),pfamily=dNormal(mu=mu,Sigma=V),weights=Claims/disp)

summary(out,dispersion=gamma.dispersion(out))
summary(out2)


V_Out=summary(out)$cov.scaled  
P_Data=solve(V_Out)

Like_std=summary(out)$coefficients[,2]

y1<-out$y
x1<-out$x
b1<-out$coefficients
wt1<-out$prior.weights
scale<-0.1275 # SAS estimate

## Should the above model be re-run with dispersion to get correct estimate 
## for the Precision matrix?

dispersion<-1/scale
wt2<-wt1/dispersion
alpha1<-rep(0,length(y1))

mu<-matrix(0,8)
P<-0.1*solve(diag(Like_std*Like_std))
V0<-solve(P)

mu[1,1]=-1.1

prior_list=list(mu=mu,Sigma=solve(P),dispersion=dispersion)
out2<-rglmb(n = 1000, out$y, out$x, prior=prior_list, wt = out$prior.weights, 
family = Gamma(link="log"), offset2 = rep(0, 20), start = out$coefficients, Gridtype = 3) 




summary(out2)
mean(out2$iters)


################################################################################

## ~ 2.381 candidates per iid sample [Consistent with theory]

###  Dispersion Model  ######################

a0<-0.01
m0<-0.01
b<-colMeans(out2$coefficients)

out3<-rglmbdisp(n=10000,y=out$y,x=out$x,b=b,alpha= rep(0, length(out$y)),wt=out$prior.weights,shape=m0,rate=a0,family=Gamma(link=log))

### Need to create a class for rglmbdisp so that print and summary produces
### more meaningful outputs 

summary(out3)


