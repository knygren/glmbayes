library(glmbayes)

### Seeds Example

germinated<-c(10,17,32,10,3 ,3,23,5,46,8, 22,23,53,10,23,15,26,55,8,0,32)
seedcount <-c(39,39,51,30,12,7,62,6,79,28,41,81,74,13,45,30,51,72,16,4,51)
seedtype  <-c(0,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,0,0,1,1,1)
extract   <-c(0,0,1,0,1,1,0,1,1,0,1,0,1,1,0,1,0,1,0,0,1)



out<-glm(cbind(germinated, seedcount-germinated) ~ seedtype+extract+seedtype*extract,family=binomial(logit), x=TRUE)

out

y<-out$y

#epsilon<-0.1

#y<-(1-epsilon)*y+epsilon*sum(germinated)/sum(seedcount)


x<-out$x
b1<-out$coefficients
wt1<-out$prior.weights
dispersion<-1
wt2<-wt1/dispersion
alpha1<-rep(0,length(y))
mu<-matrix(0,4)
P_0<-1*diag(4)


# Initialize Simulation
n<-1000

betaout<-matrix(0,nrow=n,ncol=length(y))
alphaout<-matrix(0,nrow=n,ncol=length(mu))
xtemp<-diag(1) 
P<-10



###################################  v3 ###############    

P<-14.02
sqrt(1/P)

# Smaller P improves convergence


P<-10 # Seems to impact on beta_constant but not trace constant, lambda star, or epsilon1 

# Smaller P allows tstar to be smaller so that adjusted trace_constant and adjusted gammstar
# can be kept closer to unadjusted values


pwt<-0.1

m0<-pwt/(1-pwt)


# Trace constant gets larger as m0 decreases
# epsilon1 gets smaller as m0 decreases

lambdastar<-1/(m0+1)   # Lambdastar will equal this (1-prior weight)

#lambdastar


#P_0<-diag(4)
P_0<-as.matrix(m0*P*t(x)%*%x)
P_0



num1<-length(P_0[,1])
denom1<-(m0+1)-((m0+1)/(m0+2))

trace_constant<-num1/denom1

trace_constant

gammastar_Lower<-trace_constant/(1-lambdastar)

gammastar_Lower



epsilon1<-sqrt(((m0+1)-((m0+1)/(m0+2)))^4/(1+m0)^4)
epsilon1

epsilon1*exp(-gammastar_Lower)


1-epsilon1*exp(-gammastar_Lower)


################################### Logit #################

 
qc1<-rglmb_rand(n=1000,y=y,x=x,mu=mu,P_0=P_0,P=P,wt=wt2,dispersion=dispersion,
              nu=NULL,V=NULL,family=binomial(logit),offset2=alpha1,start=mu,Gridtype=3,
              epsilon_converge=0.01)





qc1$randPostMode


trace_constant<-qc1$simconstants$trace_const
lambda_star<-qc1$simconstants$lambda_star
gammastar<-qc1$simconstants$gammastar

mu_constant<-qc1$simconstants$mu_constant
               
tstar<-qc1$simconstants$tstar

epsilonstar<-qc1$simconstants$epsilonstar
rstar<-qc1$simconstants$rstar
nstar<-qc1$simconstants$nstar

U_out<-(1+trace_constant*(1+tstar)+2*lambda_star*gammastar*(1+tstar))
U_out

alpha_out<-(1+gammastar*(1+tstar))/(1+trace_constant*(1+tstar)+lambda_star*gammastar*(1+tstar))

alpha_out



summary(qc1)
summary(qc1$randcoefficients)

library(coda)


mcmcout<-mcmc(qc1$coefficients)
plot(mcmcout)

effectiveSize(mcmcout)
autocorr.plot(mcmcout)


mcmcout2<-mcmc(qc1$randcoefficients)
plot(mcmcout2)


effectiveSize(mcmcout2)
autocorr.plot(mcmcout2)





end.timeE<-Sys.time()
time.takenE<-end.timeE-start.timeE
print("Time taken for C Version")
print(time.takenE)

m1<-matrix(0,nrow=4)
s1<-matrix(0,nrow=4)

for(i in 1:4){
m1[i,1]<-mean(qc1$alphaout[1001:10000,i])
s1[i,1]<-sd(qc1$alphaout[1001:10000,i])
  
}

frame1<-data.frame(m1,s1)

frame1

table(mean=mean1,sd=stdev1)



#############################  Probit #############################

P<-14.02
sqrt(1/P)

# Smaller P improves convergence


#P<-8 # Seems to impact on beta_constant but not trace constant, lambda star, or epsilon1 

# Smaller P allows tstar to be smaller so that adjusted trace_constant and adjusted gammstar
# can be kept closer to unadjusted values


pwt<-0.1

m0<-pwt/(1-pwt)


# Trace constant gets larger as m0 decreases
# epsilon1 gets smaller as m0 decreases

lambdastar<-1/(m0+1)   # Lambdastar will equal this (1-prior weight)

#lambdastar


#P_0<-diag(4)
P_0<-as.matrix(m0*P*t(x)%*%x)
P_0



num1<-length(P_0[,1])
denom1<-(m0+1)-((m0+1)/(m0+2))

trace_constant<-num1/denom1

trace_constant

gammastar_Lower<-trace_constant/(1-lambdastar)

gammastar_Lower



epsilon1<-sqrt(((m0+1)-((m0+1)/(m0+2)))^4/(1+m0)^4)
epsilon1

epsilon1*exp(-gammastar_Lower)


1-epsilon1*exp(-gammastar_Lower)




## Need to investigate why a large number of additional iterations end up being required
## beta_Constant ends up much larger second time around


start.timeE<-Sys.time()

qc1<-rglmb_rand(n=1000,y=y,x=x,mu=mu,P_0=P_0,P=P,wt=wt2,dispersion=dispersion,
                nu=NULL,V=NULL,family=binomial(probit),offset2=alpha1,start=mu,Gridtype=3,
                epsilon_converge=0.01)

end.timeE<-Sys.time()


summary(qc1)
summary(qc1$randcoefficients)

library(coda)


mcmcout<-mcmc(qc1$coefficients)
plot(mcmcout)

effectiveSize(mcmcout)
autocorr.plot(mcmcout)


mcmcout2<-mcmc(qc1$randcoefficients)
plot(mcmcout2)


effectiveSize(mcmcout2)
autocorr.plot(mcmcout2)


time.takenE<-end.timeE-start.timeE
print("Time taken for C Version")
print(time.takenE)

for(i in 1:4){
  m1[i,1]<-mean(qc1$alphaout[1001:10000,i])
  s1[i,1]<-sd(qc1$alphaout[1001:10000,i])
  
}

frame2<-data.frame(m1,s1)
frame2




##############################   Cloglog ######################################

P<-14.02
sqrt(1/P)

# Smaller P improves convergence


P<-5 # Seems to impact on beta_constant but not trace constant, lambda star, or epsilon1 

# Smaller P allows tstar to be smaller so that adjusted trace_constant and adjusted gammstar
# can be kept closer to unadjusted values


pwt<-0.1

m0<-pwt/(1-pwt)


# Trace constant gets larger as m0 decreases
# epsilon1 gets smaller as m0 decreases

lambdastar<-1/(m0+1)   # Lambdastar will equal this (1-prior weight)

#lambdastar


#P_0<-diag(4)
P_0<-as.matrix(m0*P*t(x)%*%x)
P_0



num1<-length(P_0[,1])
denom1<-(m0+1)-((m0+1)/(m0+2))

trace_constant<-num1/denom1

trace_constant

gammastar_Lower<-trace_constant/(1-lambdastar)

gammastar_Lower



epsilon1<-sqrt(((m0+1)-((m0+1)/(m0+2)))^4/(1+m0)^4)
epsilon1

epsilon1*exp(-gammastar_Lower)


1-epsilon1*exp(-gammastar_Lower)


start.timeE<-Sys.time()

qc1<-rglmb_rand(n=1000,y=y,x=x,mu=mu,P_0=P_0,P=P,wt=wt2,dispersion=dispersion,
                nu=NULL,V=NULL,family=binomial(cloglog),offset2=alpha1,start=mu,Gridtype=3,
                epsilon_converge=0.01)

end.timeE<-Sys.time()





time.takenE<-end.timeE-start.timeE
print("Time taken for C Version")
print(time.takenE)

for(i in 1:4){
  m1[i,1]<-mean(qc1$alphaout[1001:10000,i])
  s1[i,1]<-sd(qc1$alphaout[1001:10000,i])
  
}

frame3<-data.frame(m1,s1)

frame3


-5.19748e-010
