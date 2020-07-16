data(carinsca)
carinsca$Merit <- ordered(carinsca$Merit)
carinsca$Class <- factor(carinsca$Class)
options(contrasts=c("contr.treatment","contr.treatment"))
Claims=carinsca$Claims
Insured=carinsca$Insured
Merit=carinsca$Merit
Class=carinsca$Class
Cost=carinsca$Cost

ps=Prior_Setup(Claims/Insured~Merit+Class)
mu=ps$mu
V=ps$Sigma
mu[1,1]=log(mean(Claims/Insured))

Prior_Check(Claims/Insured~Merit+Class,family="poisson",pfamily=dNormal(mu=mu,Sigma=V))

## The # of Insured is very large

out1=glm(Claims/Insured~Merit+Class,family="poisson",weights=Insured,x=TRUE,y=TRUE)
summary(out1)
out2 <- glmb(Claims/Insured~Merit+Class,family="poisson",dNormal(mu=mu,Sigma=V),weights=Insured)
summary(out2)

#mu=out1$coefficients
#V=vcov(out1)

out1b=glm(Claims~offset(Insured)+Merit+Class,family="poisson",weights=Insured,x=TRUE,y=TRUE)
out1b=glm(Claims~offset(Insured)+Merit+Class,family="poisson",x=TRUE,y=TRUE)

log(Claims/Insured)
b=out1$coefficients
x=out1$x
y=out1$y
alpha=rep(0,length(y))  

f2_temp(b,y,x,mu,solve(V),alpha=rep(0,length(y)),wt=Insured)  
lambda<-t(exp(alpha+x%*%b))
wt=Insured
P=solve(V)
-sum(dpois(y, lambda,log=TRUE)*wt)+0.5*t((b-mu))%*%P%*%(b-mu)
-sum(dpois2(y, lambda,log=TRUE)*wt)+0.5*t((b-mu))%*%P%*%(b-mu)


sum(dpois2(y, lambda,log=TRUE))

y2=round(100*y)
lambda2=100*lambda

dpois(y2, lambda2,log=TRUE)
(-lambda2+y2*log(lambda2)-log(gamma(y2+1)))

(-lambda+y*log(lambda)-log(gamma(y+1)))*wt

y*log(lambda)
is.integer(y)

dpois(y, lambda,log=FALSE)

dpois2<-function(x,lambda,log=TRUE){
  
  if(is.integer(x)) return(dpois(x,lambda,log=TRUE))
  else{warning("Non-Integer Values to Poisson Density - Switching to Gamma Function to Evaluate Factorial")
    return(-lambda+x*log(lambda)-log(gamma(x+1)))
    
  }
}






f2_temp<-function(b,y,x,mu,P,alpha=0,wt=1){
  lambda<-t(exp(alpha+x%*%b))
  -sum(dpois(y, lambda,log=TRUE)*wt)+0.5*t((b-mu))%*%P%*%(b-mu)
}


log(Claims)

solve(V)

summary(out1,cor=F)

out1$prior.weights




data(carinsca)
carinsca$Merit <- ordered(carinsca$Merit)
carinsca$Class <- factor(carinsca$Class)
options(contrasts=c("contr.treatment","contr.treatment"))


scale<-0.1275 # SAS estimate

dispersion<-1/scale

out <- glm(Cost/Claims~Merit+Class,family=Gamma(link="log"),weights=Claims,x=TRUE)
summary(out)

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


