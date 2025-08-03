# Gelman et. al (2003) Schools data
school <- list("A","B","C","D","E","F","G","H")
estimate<- c(28.39, 7.94, -2.75 , 6.82, -0.64, 0.63, 18.01, 12.16)
sd <- c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6)

shape
rate
J<-length(school)


theta<-rep(0,J)
mu.theta<-mean(estimate)
sigma.theta<-runif(1,0,100)
sigma.theta_2=sigma.theta^2
#sigma.mu=1/1.0e-06
sigma.mu=((max(estimate)-min(estimate))/2)^2
sqrt(sigma.mu)
1/sigma.mu
wt=1/sigma.y_2




## Modify the data

y=estimate
sd=sd/10  # Modifid data precision (Maybe temporary)


################################# Independent_Normal_Gamma Prior  #######################

## Prior Values

# Normal Component

mu.mu=mean(estimate) 
sigma.mu=var(estimate)

## Gamma Component

n_prior=0.5
disp_ML=var(estimate)
sqrt(disp_ML)
shape=n_prior/2
rate= disp_ML*shape

# Data Precision Component 

# set number of draws to generate during burn-in and simulation

sim_1=100
sim_2=1000

# Set up objects holding output 
theta_out<-matrix(0,nrow=sim_2,ncol=J)
mu_out<-rep(0,sim_2)
sigma_theta_out<-rep(0,sim_2)
disp_out<-rep(0,sim_2)
iters_out1=rep(0,sim_1)
iters_out2=rep(0,sim_2)

# Initial value for the Two-Block Gibbs sampler 

theta=c(27.93,7.7978,-2.545,6.814,-0.5745,0.8067,17.94,12.07)


for(k in 1:sim_1){

  out2<-rlmb(1,y=theta,x=as.matrix(rep(1,J),nrow=J,ncol=1),
             pfamily=dIndependent_Normal_Gamma(mu.mu,sigma.mu,shape=shape,rate=rate))
  
  mu.theta=out2$coefficients
  sigma_theta_2=out2$dispersion
  
  for(j in 1:J){
  
  theta[j]<-rlmb(1,y[j],as.matrix(1),
                       pfamily=dNormal(mu.theta,Sigma=sigma.theta_2,dispersion=sigma.y_2[j]))$coefficients
  

}

iters_out1[k]=out2$iters
}



mean(iters_out1)  # Much better acceptance rate once prior mean more sensible




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
  disp_out[k]=out2$dispersion
  sigma_theta_out[k]=sqrt(sigma_theta_2)
  iters_out2[k]=out2$iters
  
}


mean(iters_out2)   ## 58.167 draws per acceptance - not too bad - should be much faster if *.cpp code was used
colMeans(theta_out)
mean(mu_out)
mean(disp_out)
mean(sigma_theta_out)
sqrt(diag(var(theta_out)))

######################################## Normal_Gamma Prior  ##########################


## Prior Values

# Normal Component

mu.mu=mean(estimate) 
sigma.mu=var(estimate)

## Gamma Component

n_prior=0.5
disp_ML=var(estimate)
sqrt(disp_ML)
shape=n_prior/2
rate= disp_ML*shape

# Data Precision Component 

# set number of draws to generate during burn-in and simulation

sim_1=100
sim_2=1000

# Set up objects holding output 
theta_out2<-matrix(0,nrow=sim_2,ncol=J)
mu_out2<-rep(0,sim_2)
sigma_theta_out2<-rep(0,sim_2)
disp_out2<-rep(0,sim_2)

# Initial value for the Two-Block Gibbs sampler 

theta=c(27.93,7.7978,-2.545,6.814,-0.5745,0.8067,17.94,12.07)


for(k in 1:sim_1){
  
  out2<-rlmb(1,y=theta,x=as.matrix(rep(1,J),nrow=J,ncol=1),
             pfamily=dNormal_Gamma(mu.mu,sigma.mu/disp_ML,shape=shape,rate=rate))
  
  mu.theta=out2$coefficients
  sigma_theta_2=out2$dispersion
  
  for(j in 1:J){
    
    theta[j]<-rlmb(1,y[j],as.matrix(1),
                   pfamily=dNormal(mu.theta,Sigma=sigma.theta_2,dispersion=sigma.y_2[j]))$coefficients
    
  }
  

}



for(k in 1:sim_2){
  
  out2<-rlmb(1,y=theta,x=as.matrix(rep(1,J),nrow=J,ncol=1),
             pfamily=dNormal_Gamma(mu.mu,sigma.mu/disp_ML,shape=shape,rate=rate))
  
  mu.theta=out2$coefficients
  sigma_theta_2=out2$dispersion
  
  
  for(j in 1:J){
    
    theta[j]<-rlmb(1,y[j],as.matrix(1),
                   pfamily=dNormal(mu.theta,Sigma=sigma.theta_2,dispersion=sigma.y_2[j]))$coefficients
    
  }
  
  theta_out2[k,1:J]=theta[1:J]
  mu_out2[k]=mu.theta
  disp_out2[k]=out2$dispersion
  sigma_theta_out2[k]=sqrt(sigma_theta_2)

}


colMeans(theta_out2)
mean(mu_out2)
mean(disp_out2)
mean(sigma_theta_out2)
sqrt(diag(var(theta_out2)))


######################################## Normal Prior  ##########################
#### The distance between samples and true density here is very simple and can
#### be very easily bounded
#### Same holds when multiple explanatory variables in hierarchy
#### Or when lower level had multiple components (say if it is a regression itself)

## Prior Values

# Normal Component

mu.mu=mean(estimate) 
sigma.mu=var(estimate)

## Gamma Component

n_prior=0.5
#disp_ML=var(estimate)
#shape=n_prior/2
#rate= disp_ML*shape
disp_Post=129.6592

# Data Precision Component 

# set number of draws to generate during burn-in and simulation

sim_1=100
sim_2=1000

# Set up objects holding output 
theta_out3<-matrix(0,nrow=sim_2,ncol=J)
mu_out3<-rep(0,sim_2)
sigma_theta_out3<-rep(0,sim_2)
disp_out3<-rep(0,sim_2)

# Initial value for the Two-Block Gibbs sampler 

theta=c(27.93,7.7978,-2.545,6.814,-0.5745,0.8067,17.94,12.07)


for(k in 1:sim_1){
  
  out2<-rlmb(1,y=theta,x=as.matrix(rep(1,J),nrow=J,ncol=1),
             pfamily=dNormal(mu.mu,sigma.mu,dispersion=disp_Post))
  
  mu.theta=out2$coefficients
#  sigma_theta_2=out2$dispersion
  sigma_theta_2=disp_Post
  
  for(j in 1:J){
    
    theta[j]<-rlmb(1,y[j],as.matrix(1),
                   pfamily=dNormal(mu.theta,Sigma=sigma.theta_2,dispersion=sigma.y_2[j]))$coefficients
    
  }
  
  
}



for(k in 1:sim_2){
  
  out2<-rlmb(1,y=theta,x=as.matrix(rep(1,J),nrow=J,ncol=1),
             pfamily=dNormal(mu.mu,sigma.mu,dispersion=disp_Post))
  
  mu.theta=out2$coefficients
#  sigma_theta_2=out2$dispersion
  sigma_theta_2=disp_Post
  
  for(j in 1:J){
    
    theta[j]<-rlmb(1,y[j],as.matrix(1),
                   pfamily=dNormal(mu.theta,Sigma=sigma.theta_2,dispersion=sigma.y_2[j]))$coefficients
    
  }
  
  theta_out3[k,1:J]=theta[1:J]
  mu_out3[k]=mu.theta
  disp_out3[k]=sigma_theta_2
  sigma_theta_out3[k]=sqrt(sigma_theta_2)
  
}


colMeans(theta_out3)
mean(mu_out3)
mean(disp_out3)
mean(sigma_theta_out3)
sqrt(diag(var(theta_out3)))
