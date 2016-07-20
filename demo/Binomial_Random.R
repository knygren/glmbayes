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


P<-0.1 # Seems to impact on beta_constant but not trace constant, lambda star, or epsilon1 

# Smaller P allows tstar to be smaller so that adjusted trace_constant and adjusted gammstar
# can be kept closer to unadjusted values


pwt<-0.6

m0<-pwt/(1-pwt)


# Trace constant gets larger as m0 decreases
# epsilon1 gets smaller as m0 decreases

lambdastar<-1/(m0+1)   # Lambdastar will equal this (1-prior weight)

#lambdastar


#P_0<-diag(4)
P_0<-as.matrix(m0*P*t(x)%*%x)
P_0

#det(P*t(x)%*%x)


#det(P_0+P*t(x)%*%x)

#det(P_0+P*t(x)%*%x)/det(P*t(x)%*%x)

#(1+m0)^4
#((m0+1)-((m0+1)/(m0+2)))^4




#P_0<-0.01*P_0

# Note: Larger P_0 improves convergence
#P_0<-100*P_0


# lg_A1_out_New Fails for this 

#P_0<-4*P_0

# lg_A1_out Fails for this 

#P_0<-10*P_0

# Lg_alpha_out fails for this

#P_0<-0.0001*P_0

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



ropt<-function(rstar,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant){
d1<-(1-epsilonstar)^rstar
d2<-(U_out^rstar)/(alpha_out^(1-rstar))
d1^nstar+(d2^nstar)*(1+gammastar*(1+tstar)+mu_constant*(1+tstar))
}


old_opt<-ropt(rstar,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant)

low<-old_opt


# Initialize sc_temp

sc_temp1<-0.001


sc_temp<-(sc_temp1/(1+sc_temp1))*rstar





golden.section.search = function(){

golden.ratio = 2/(sqrt(5) + 1)
upper.bound<-rstar
lower.bound<-0
tolerance<-0.00001


### Use the golden ratio to set the initial test points
r1 = upper.bound - golden.ratio*(upper.bound-lower.bound)
r2 = lower.bound + golden.ratio*(upper.bound - lower.bound)


ropt(rstar,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant)
f1<-ropt(r1,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant)
f2<-ropt(r2,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant)

iteration = 0

while (abs(upper.bound - lower.bound) > tolerance){

  iteration = iteration + 1
  cat('', '\n')
  cat('Iteration #', iteration, '\n')
  cat('f1 =', f1, '\n')
  cat('f2 =', f2, '\n')
  
  if (f2 > f1){
  # then the minimum is to the left of r2
  # let r2 be the new upper bound
  # let r1 be the new upper test point
    cat('f2 > f1', '\n')
    ### Set the new upper bound
  ### Set the new upper bound
  upper.bound <- r2  
  cat('New Upper Bound =', upper.bound, '\n')
  cat('New Lower Bound =', lower.bound, '\n')
  
  ### Set the new upper test point
  ### Use the special result of the golden ratio
  r2 = r1
  cat('New Upper Test Point = ', r2, '\n')
  f2 = f1

  ### Set the new lower test point
  r1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
  cat('New Lower Test Point = ', r1, '\n')
  f1 <- ropt(r1,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant)
    
}

else
  
{
  cat('f2 < f1', '\n')
  # the minimum is to the right of x1
  # let x1 be the new lower bound
  # let x2 be the new lower test point
  
  ### Set the new lower bound
  lower.bound = r1
  cat('New Upper Bound =', upper.bound, '\n')
  cat('New Lower Bound =', lower.bound, '\n')
  
  ### Set the new lower test point
  r1 = r2
  cat('New Lower Test Point = ', r1, '\n')
  
  f1 = f2
  
  ### Set the new upper test point
  r2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
  cat('New Upper Test Point = ', r2, '\n')
  f2 = ropt(r2,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant)
}

}


  
  
}









sc_temp


rdown_bound_flag<-0
  
rup<-rstar

rprior<-rup

rstar_new<-rstar
  
while(rdown_bound_flag==0){

rdown<-rstar-sc_temp

down<-ropt(rdown,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant)
diff_down<-down-old_opt


if(diff_down<0)
{
  if(down<low){
  low<-down
  rup<-rprior
  rstar_new<-rdown
  }
  sc_temp1<-sc_temp1*10

  sc_temp<-(sc_temp1/(1+sc_temp1))*rstar
  
  }

if(diff_down>0)
{
  rdown_bound<-rdown
  rdown_bound_flag<-1
}


rprior<-rdown


}


rstar
rup
rstar_new
rdown

rdown_bound

gold_ratio<-(1+sqrt(5))/2

gold_ratio

rstart<-rdown_bound+(1/(1+gold_ratio))*(rup-rdown_bound)


rstart


(rup-rstart)/(rstart-rdown_bound)


sc_temp


old_opt
low


rdown_bound_flag
sc_temp


# Iteration 2


rdown<-rstar-sc_temp

down<-ropt(rdown,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant)
diff_down<-down-old_opt

diff_down


if(diff_down<0)
{
  low<-down
  
  sc_temp<-sc_temp*10
}

if(diff_down>0)
{
  rdown_bound<-rdown
  rdown_bound_flag<-1
}


rdown_bound_flag
sc_temp


# Iteration 3


rdown<-rstar-sc_temp

down<-ropt(rdown,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant)
diff_down<-down-old_opt

diff_down


if(diff_down<0)
{
  low<-down
  
  sc_temp<-sc_temp*10
}

if(diff_down>0)
{
  rdown_bound<-rdown
  rdown_bound_flag<-1
}


rdown_bound_flag
sc_temp



# Iteration 4


rdown<-rstar-sc_temp

down<-ropt(rdown,epsilonstar,U_out,alpha_out,nstar,gammastar,tstar,mu_constant)
diff_down<-down-old_opt

diff_down


if(diff_down<0)
{
  low<-down
  
  sc_temp<-sc_temp*10
}

if(diff_down>0)
{
  rdown_bound<-rdown
  rdown_bound_flag<-1
}


rdown_bound_flag
sc_temp



#################################################################33





min(diff1,diff2)



#Change due to rstar

(d1^nstar)*nstar*log(1-epsilonstar)
+(1+gammastar*(1+tstar)+mu_constant*(1+tstar))((1/(alpha_out^(1-rstar)))*log(U_out))

summary(qc1)


library(coda)


mcmcout<-mcmc(qc1$coefficients)
plot(mcmcout)



mcmcout2<-mcmc(qc1$randcoefficients)
plot(mcmcout2)





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

start.timeE<-Sys.time()

qc1<-rglmb_rand(n=11000,y=y,x=x,mu=mu,P_0=P_0,P=P,wt=wt2,dispersion=dispersion,
                nu=NULL,V=NULL,family=binomial(probit),offset2=alpha1,start=mu,Gridtype=3)

end.timeE<-Sys.time()
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

start.timeE<-Sys.time()

qc1<-rglmb_rand(n=11000,y=y,x=x,mu=mu,P_0=P_0,P=P,wt=wt2,dispersion=dispersion,
                nu=NULL,V=NULL,family=binomial(cloglog),offset2=alpha1,start=mu,Gridtype=3)

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
