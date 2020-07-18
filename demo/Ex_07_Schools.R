# Gelman et. al (2003) Schools data
school <- list("A","B","C","D","E","F","G","H")
estimate<- c(28.39, 7.94, -2.75 , 6.82, -0.64, 0.63, 18.01, 12.16)
sd <- c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6)

sd=sd/10

disp_ML=var(estimate)
sqrt(disp_ML)
n_prior=0.5
shape=n_prior/2
rate= disp_ML*shape
shape
rate
J<-length(school)

y=estimate
sigma.y_2<-sd^2

theta<-rep(0,J)
mu.theta<-mean(estimate)
sigma.theta<-runif(1,0,100)
sigma.theta_2=sigma.theta^2
#sigma.mu=1/1.0e-06
sigma.mu=((max(estimate)-min(estimate))/2)^2
sigma.mu=var(estimate)
sqrt(sigma.mu)

1/sigma.mu
wt=1/sigma.y_2
mu.mu=mean(estimate)


sim_1=100
sim_2=1000
theta_out<-matrix(0,nrow=sim_2,ncol=J)
mu_out<-rep(0,sim_2)
sigma_theta_out<-rep(0,sim_2)

iters_out1=rep(0,sim_1)
iters_out2=rep(0,sim_2)

theta=c(11.22,7.766,6.476,7.513,5.519,6.287,10.33,8.681)


for(k in 1:sim_1){

  out2<-rlmb(1,y=theta,x=as.matrix(rep(1,J),nrow=J,ncol=1),
             pfamily=dIndependent_Normal_Gamma(mu.mu,sigma.mu,shape=shape,rate=rate))
  
  mu.theta=out2$coefficients
  sigma_theta_2=out2$dispersion
  
  for(j in 1:J){
  
  theta[j]<-rlmb(1,y[j],as.matrix(1),
                       pfamily=dNormal(mu.theta,Sigma=sigma.theta_2,dispersion=sigma.y_2[j]))$coefficients
  

}

  theta_out[k,1:J]=theta[1:J]
  mu_out[k]=mu.theta
  sigma_theta_out[k]=sqrt(sigma_theta_2)
  iters_out2[k]=out2$iters
  

iters_out1[k]=out2$iters
}



mean(iters_out1)  # Much better acceptance rate once prior mean more sensible

theta ## Way too large for some values
theta_out[1:100,]
mu_out[1:100]
mean(sigma_theta_out[1:100])



for(k in 1:sim_2){

  out2<-rlmb(1,y=theta,x=as.matrix(rep(1,J),nrow=J,ncol=1),
             pfamily=dIndependent_Normal_Gamma(mu.mu,sigma.mu,shape=shape,rate=rate))
  
  mu.theta=out2$coefficients
  sigma_theta_2=out2$dispersion
  
  
    for(j in 1:J){
    
    theta[j]<-rlmb(1,y[j],as.matrix(1),
                    pfamily=dNormal(mu.theta,Sigma=sigma.theta_2,dispersion=sigma.y_2[j]))$coefficients
  
  }

  theta_out[k,1:J]=theta[1:J]
  mu_out[k]=mu.theta
  sigma_theta_out[k]=sqrt(sigma_theta_2)
  iters_out2[k]=out2$iters
  
}


mean(iters_out2)


colMeans(theta_out)
mean(mu_out)
mean(sigma_theta_out)

sqrt(diag(var(theta_out)))


################################################ Testing:

### Step 1: Random effects code

mu_test=7.941  # From Winbugs output
sigma.theta_2_test=6.164^2
  theta_test<-matrix(0,nrow=1000,ncol=8)

for(j in 1:8){
theta_test[1:1000,j]<-rlmb(1000,y[j],as.matrix(1),
               pfamily=dNormal(mu_test,Sigma=sigma.theta_2_test,dispersion=sigma.y_2[j]))$coefficients
}

### Seems quite close to Winbugs output  
colMeans(theta_test)
sqrt(diag(var(theta_test)))

### Step 2: Fixed effects and variance code

theta=c(11.22,7.766,6.476,7.513,5.519,6.287,10.33,8.681)

out2<-rlmb(1000,y=theta,x=as.matrix(rep(1,J),nrow=J,ncol=1),
           pfamily=dIndependent_Normal_Gamma(as.matrix(0),sigma.mu,shape=shape,rate=rate))

colMeans(out2$coefficients)
mean(out2$dispersion)
mean(sqrt(out2$dispersion))

summary(sqrt(out2$dispersion))

p1temp=1/sigma.theta_2_test
p2temp=1/sigma.y_2

theta=c(11.22,7.766,6.476,7.513,5.519,6.287,10.33,8.681)
mutemp=7.941
wt1=p1temp/(p1temp+p2temp)

mutemp*wt1+theta*(1-wt1)

