library(glmbayes)

### Seeds Example

germinated<-c(10,17,32,10,3 ,3,23,5,46,8, 22,23,53,10,23,15,26,55,8,0,32)
seedcount <-c(39,39,51,30,12,7,62,6,79,28,41,81,74,13,45,30,51,72,16,4,51)
seedtype  <-c(0,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,0,0,1,1,1)
extract   <-c(0,0,1,0,1,1,0,1,1,0,1,0,1,1,0,1,0,1,0,0,1)

out<-glm(cbind(germinated, seedcount-germinated) ~ seedtype+extract+seedtype*extract,family=binomial(logit), x=TRUE)

out

y<-out$y
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
P<-1



###################################  v3 ###############    

P<-14.02
sqrt(1/P)


P_0<-diag(4)
P_0<-0.01*P_0

################################### Logit #################

 
qc1<-rglmb_rand(n=11000,y=y,x=x,mu=mu,P_0=P_0,P=P,wt=wt2,dispersion=dispersion,
              nu=NULL,V=NULL,family=binomial(logit),offset2=alpha1,start=mu,Gridtype=3)



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



