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


out1b=glm(Claims/Insured~Merit+Class,family="quasipoisson",weights=Insured,x=TRUE,y=TRUE)
summary(out1b)



out2 <- glmb(Claims/Insured~Merit+Class,family="poisson",dNormal(mu=mu,Sigma=V),weights=Insured)
summary(out2)

## This seesm to only adjust the standard errors and not the estimates themselves
## Need to update code to hanlde 

#dispersion=48.58325
#out3 <- glmb(Claims/Insured~Merit+Class,family="poisson",dNormal(mu=mu,Sigma=V),weights=Insured/dispersion)


#summary(out3)


out3 <- glmb(Claims/Insured~Merit+Class,family="quasipoisson",dNormal(mu=mu,Sigma=V),weights=Insured)


out3
summary(out2)
summary(out3)




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


