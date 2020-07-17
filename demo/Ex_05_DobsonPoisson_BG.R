set.seed(360)
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
mysd<-1
mu<-matrix(0,5)
mu[1]=log(mean(counts))
V0<-((mysd)^2)*diag(5)
glm.D93<-glm(counts ~ outcome + treatment, family = poisson())
glmb.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu,Sigma=V0))
summary(glmb.D93)
glmb.D93_v2<-glmb(counts ~ outcome + treatment, family = quasipoisson(),pfamily=dNormal(mu=mu,Sigma=V0))
summary(glmb.D93_v2)

############################  Try - Two - Block Gibbs sampler  ################################

## Get estimate from glmb.D93 model as starting point

y=glmb.D93$y
x=glmb.D93$x
n_obs=length(y)
n_var=ncol(x)
offset2=rep(0,n_obs)
wt=glmb.D93$prior.weights
bbar=glmb.D93_v2$coef.means
dispbar=glmb.D93_v2$dispersion

xbeta=offset2+x%*%bbar
b1temp=rep(0,n_obs)

n_sim1=100
n_sim2=1000

b1out=matrix(0,nrow=n_sim2,ncol=n_obs)
b2out=matrix(0,nrow=n_sim2,ncol=n_var)


ptm <- proc.time()
for(i in 1:n_sim1){
  
  # Generate estimates for random effects
  
  for(j in 1:n_obs){
    b1temp[j]=rglmb(1,y[j],as.matrix(1),family=poisson(),pfamily=dNormal(xbeta[j],dispbar),offset=0,weights=wt[j])$coefficients  
  }
  
  # Generate estimates for fixed effects
  
  b2temp=rglmb(1,b1temp,x,family=gaussian(),pfamily=dNormal(mu,V0),offset=offset2,weights=1/dispbar)$coefficients 
  xbeta=offset+x%*%t(b2temp)
}
print(proc.time()-ptm)


ptm <- proc.time()
for(i in 1:n_sim2){
  
  # Generate estimates for random effects
  
  for(j in 1:n_obs){
    b1temp[j]=rglmb(1,y[j],as.matrix(1),family=poisson(),pfamily=dNormal(xbeta[j],dispbar),offset=0,weights=wt[j])$coefficients  
  }
  
  # Generate estimates for fixed effects
  
  b2temp=rglmb(1,b1temp,x,family=gaussian(),pfamily=dNormal(mu,V0),offset=offset2,weights=1/dispbar)$coefficients 
  xbeta=offset+x%*%t(b2temp)
  
  b1out[i,1:n_obs]=b1temp[1:n_obs]
  b2out[i,1:n_var]=b2temp[1:n_var]
  
}
print(proc.time()-ptm)

rbind(
  glmb.D93$fit$coefficients,
  glmb.D93$coef.means,
glmb.D93_v2$coef.means,
colMeans(b2out))

rbind(log(y),
colMeans(predict(glmb.D93))
,colMeans(predict(glmb.D93_v2)),
colMeans(b1out) )
