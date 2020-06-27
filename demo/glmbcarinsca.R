
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


