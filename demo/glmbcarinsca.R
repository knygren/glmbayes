
data(carinsca)
carinsca$Merit <- ordered(carinsca$Merit)
carinsca$Class <- factor(carinsca$Class)
options(contrasts=c("contr.treatment","contr.treatment"))
attach(carinsca)
out <- glm(Claims/Insured~Merit+Class,family="poisson")
summary(out,cor=FALSE)

wt1<-out$prior.weights
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

Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V0)),out$coefficients,Like_std)

mu
out$coefficients

mu[1,1]=-1.1

Prior_Error_Checks=Prior_Likelihood_Check(mu,sqrt(diag(V0)),out$coefficients,Like_std)
Prior_Error_Checks

out2<-rglmb(n = 1000, out$y, out$x, mu=mu, P=P, wt = out$prior.weights, dispersion = dispersion,
            family = Gamma(link="log"), offset2 = rep(0, 20), start = out$coefficients, Gridtype = 3) 

summary(out2)

###  Dispersion Model  ######################

a0<-0.01
m0<-0.01
b<-colMeans(out2$coefficients)

out3<-rglmbdisp(n=10000,y=out$y,x=out$x,b=b,alpha= rep(0, length(out$y)),wt=out$prior.weights,shape=m0,rate=a0,family=Gamma(link=log))

### Need to create a class for rglmbdisp so that print and summary produces
### more meaningful outputs 

summary(out3)

### Development Code may have attempts at calculating the DIC (not sure aout rglmbdisp)


######  Joint Estimation - Two -Block Gibbs Sampler ###


y1<-out$y
x1<-out$x
b1<-out$coefficients
wt1<-out$prior.weights

### This part likely not used 

#alpha1<-rep(0,length(y1)) 
#mu<-matrix(0,8)
#P<-10*diag(8)

scale<-0.1275 # SAS estimate
dispersion<-1/scale

a0<-0.01
m0<-0.01

outbetas<-matrix(0,nrow=11000,ncol=8)
outdisp<-matrix(0,nrow=11000)

disp<-dispersion
start.time<-Sys.time()
for(i in 1:100)
{
  start.time1<-Sys.time()
  
  outtemp1<-rglmb(n = 1, out$y, out$x, mu=mu, P=P, wt = out$prior.weights, dispersion = disp, family = Gamma(link="log"), 
        offset2 = rep(0, 20), start = out$coefficients, Gridtype = 2) 
  end.time1<-Sys.time()
  time.taken1<-end.time1-start.time1
#  print("Time for rglmb")
#  print(time.taken1)
  
  outbetas[i,1:8]<-outtemp1$coefficients[1,1:8]
  b<-outbetas[i,1:8]  
  start.time2<-Sys.time()  
  outdisp[i,1]<-rglmbdisp(n=1,y=out$y,x=out$x,b=b,alpha= rep(0, length(out$y)),wt=out$prior.weights,shape=m0,rate=a0,family=Gamma(link=log))
  end.time2<-Sys.time()
  time.taken2<-end.time2-start.time2
#  print("Time for rglmbdisp")
#  print(time.taken2)
  
#  print(i) 
  disp<-outdisp[i,1]  
}  
end.time<-Sys.time()
time.taken<-end.time-start.time
print("Time for loop")
print(time.taken)

mean(outdisp[40:100,1])
colMeans(outbetas[40:100,])
