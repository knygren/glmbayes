rand_Norm_apprx<-function(n=1,y,x,P_0,P,mustar,bstar,alpha1,family,link){

  # Get dimensions 
  l1<-length(bstar)
  l2<-length(mustar)
  
  # Set up matrices of output 
  
  Pstar<-vector(mode = "numeric", length = l1)
  
  outsim2<-matrix(0,nrow=n,ncol=l1)

  # Get Family Functions
  
  myfamfunc<-glmbfamfunc(list(family=family,link=link))  # Make this flexible (linked to family and link)
  f7_temp<-myfamfunc$f7
  f3_temp<-myfamfunc$f3
  f2_temp<-myfamfunc$f2
  f1_temp<-myfamfunc$f1
  
  # Calculate constants for simulation
  
  xmustar<-t(t(alpha1))+x%*%mustar

  
  for(i in 1:l1){
    
    Pstar[i]<-f7_temp(bstar[i],y[i],1,xmustar[i],P)   
  }  
  

  Ppost<-Pstar+P
  
  # First error appears here
  
  qcystar1<-Ppost*bstar
  qcystar2<-P[1]*xmustar

  qcystar3<-Ppost*bstar-P[1]*xmustar
  
#    l3<-length(P)
  
  ystar<-(Ppost*bstar-P[1]*xmustar)/Pstar
  
  Pstar2<-matrix(0,nrow=l1,ncol=l1)
  
  for(i in 1:l1){
    Pstar2[i,i]<-Pstar[i]
  }
  
  

  varcovout<-solve(P_0+P[1]*t(x)%*%x -t(x)%*%(P[1]*diag(l1))%*%solve(Pstar2+P[1]*diag(l1))%*%(P[1]*diag(l1)%*%x))  # Variance -covariance from normal approximation

  outsim<-rmnorm(n = n, mean = mustar, varcov=varcovout, sqrt=NULL)
  

  # Run simulation for beta
  
  for(i in 1:n){
    
    xmutemp<-t(t(alpha1))+x%*%outsim[i,1:l2]
    
#   bstar2_mu[,j]<-(P[1]*xmutemp+Pstar*ystar)/Ppost
    bstartemp<-(P[1]*xmutemp+Pstar*ystar)/Ppost
    outsim2[i,]<-rnorm(n=l1,mean=bstartemp,sd=sqrt(1/Ppost))
    
  }
  
  
  return(list(mu_sim=outsim,beta_sim=outsim2,Pstar=Pstar,ystar=ystar,Ppost=Ppost))
  
#  return(list(n=n,y=y,x=x,P_0=P_0,P=P,mustar=mustar,bstar,alpha1,family,link))
}





KL_Norm_apprx<-function(y,x,P,Pstar,ystar,Ppost,mu_sim,family,link,eps_in){

#KL_Norm_apprx<-function(y,x,P,Pstar,ystar,Ppost,mu_sim,family,link){
  
  # Get dimensions
  
  n<-length(mu_sim[,1])
  l1<-length(x[,1])
  
#  eps_in<-matrix(1,nrow=l1,ncol=1) # figure out how to move outside of function
  
  
  
  # Setup output matrices
  
  bstar2_mu<-matrix(0,nrow=l1,ncol=n) # Matrix of conditional posterior modes (and means) of normal approximation for posterior density
  bstar_mu<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  bstar_mu_eps<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  bstar_mu_half<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  
  Pstar_mu<-matrix(0,nrow=l1,ncol=n)
  Pstar2_mu<-matrix(0,nrow=l1,ncol=n)
  Pstar2_mu_eps<-matrix(0,nrow=l1,ncol=n)
  Pstar_mu_half<-matrix(0,nrow=l1,ncol=n)
  
  
  
  P_eps<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  P_eps_half<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  
  t1<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  t2<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  t3<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  t4<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  
  t5<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  t6<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  
  t7<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  
  # Get family functions
  
  myfamfunc<-glmbfamfunc(list(family=family,link=link))  # Make this flexible (linked to family and link)
  f7_temp<-myfamfunc$f7
  f3_temp<-myfamfunc$f3
  f2_temp<-myfamfunc$f2
  f1_temp<-myfamfunc$f1
  
  # Loop through iterations
  
  for(j in 1:n){
    
    # Calculate conditional means for b under normal approximation
    
    xmutemp<-t(t(alpha1))+x%*%mu_sim[j:j,]
    
    
    bstar2_mu[,j]<-(P[1]*xmutemp+Pstar*ystar)/Ppost
    
    # use bstar2_mu[,j] as starting point for finding bstar_mu[,j]
    
    bstar_mu[,j]<-bstar2_mu[,j]
    
    bstar_mu_eps[,j]<-bstar2_mu[,j]
    bstar_mu_half[,j]<-bstar2_mu[,j]
    
    # Loop Through dimensions

    for(i in 1:l1){
      
      # Leverage newtown's method (currently five iterations) to find bstar_mu
      
      for(k in 1:5){  
        a<-f3_temp(bstar_mu[i,j],y[i],1,xmutemp[i],P)
        
        b<-P+f7_temp(bstar_mu[i,j],y[i],1,xmutemp[i],P)   # need to modify to use Precision determined by family and link
        
        bstar_mu[i,j]<-bstar_mu[i,j]-(a/b)
      }
      
      for(k in 1:5){  
        a<-(1-(eps_in[i]))*f3_temp(bstar_mu_eps[i,j],y[i],1,xmutemp[i],P)+((eps_in[i]))*(P%*%(bstar_mu_eps[i,j]-xmutemp[i])+Pstar[i]%*%(bstar_mu_eps[i,j]-ystar[i]))
        b<-P+(1-(eps_in[i]))*f7_temp(bstar_mu_eps[i,j],y[i],1,xmutemp[i],P)+(eps_in[i])*Pstar[i]   # need to modify to use Precision determined by family and link
        
        
        bstar_mu_eps[i,j]<-bstar_mu_eps[i,j]-(a/b)
      }
      
      
      for(k in 1:5){  
        a<-(1-(eps_in[i]/2))*f3_temp(bstar_mu_half[i,j],y[i],1,xmutemp[i],P)+((eps_in[i]/2))*(P%*%(bstar_mu_half[i,j]-xmutemp[i])+Pstar[i]%*%(bstar_mu_half[i,j]-ystar[i]))
        b<-P+(1-(eps_in[i]/2))*f7_temp(bstar_mu_half[i,j],y[i],1,xmutemp[i],P)+(eps_in[i]/2)*Pstar[i]   # need to modify to use Precision determined by family and link
        
        
        bstar_mu_half[i,j]<-bstar_mu_half[i,j]-(a/b)
      }
      
      
      Pstar_mu[i,j]<-f7_temp(bstar_mu[i,j],y[i],1,xmutemp[i],P) 
      Pstar2_mu[i,j]<-f7_temp(bstar2_mu[i,j],y[i],1,xmutemp[i],P) 
      Pstar2_mu_eps[i,j]<-f7_temp(bstar_mu_eps[i,j],y[i],1,xmutemp[i],P) 
      Pstar_mu_half[i,j]<-f7_temp(bstar_mu_half[i,j],y[i],1,xmutemp[i],P) 
      
      
        }
    
    
        
      }

  
  for(j in 1:n){
    for(i in 1:l1){
      
      P_eps[i,j]<-P+(1-eps_in[i])*Pstar2_mu_eps[i,j]+eps_in[i]*Pstar[i]
      P_eps_half[i,j]<-P+(1-(eps_in[i]/2))*Pstar_mu_half[i,j]+(eps_in[i]/2)*Pstar[i]
      
      t1[i,j]<-0.5*(1-eps_in[i])*(1/P_eps[i,j])*(Pstar_mu[i,j]-Pstar2_mu_eps[i,j])
      t2[i,j]<-0.5*(eps_in[i])*(1/P_eps[i,j])*(Pstar_mu[i,j]-Pstar[i])
      
      t3[i,j]<-0.5*(1-eps_in[i])*(1/P_eps_half[i,j])*(Pstar_mu[i,j]-Pstar2_mu_eps[i,j])
      t4[i,j]<-0.5*(eps_in[i])*(1/P_eps_half[i,j])*(Pstar_mu[i,j]-Pstar[i])
      
      t5[i,j]<-0.5*(bstar_mu_eps[i,j]-bstar_mu[i,j])*(P+Pstar_mu[i,j])*(bstar_mu_eps[i,j]-bstar_mu[i,j])
      t6[i,j]<-0.5*(bstar_mu_half[i,j]-bstar_mu[i,j])*(P+Pstar_mu[i,j])*(bstar_mu_half[i,j]-bstar_mu[i,j])
      
      t7[i,j]<-0.5*(bstar_mu_half[i,j]-bstar_mu_eps[i,j])*((1-eps_in[i])*(P+Pstar2_mu_eps[i,j])+(eps_in[i])*(P+Pstar[i]))*(bstar_mu_half[i,j]-bstar_mu_eps[i,j])
      
    }
  }

  
  KL_out1<-t(t1)+t(t2)-t(t3)-t(t4)+t(t5)-t(t6)+t(t7)
  KL_out2<-matrix(0,nrow=n,ncol=1)
  KL_out2b<-matrix(0,nrow=1,ncol=l1)
  
  for(i in 1:n){
    KL_out2[i]<-sum(KL_out1[i,1:l1])
  }
  
  for(j in 1:l1){
    KL_out2b[1,j]<-mean(KL_out1[,j])
  }
  
  
  KL<-mean(KL_out2)
  
  TV<-sqrt(0.5*KL)
  
  
    
  
#  return(list(bstar2_mu=bstar2_mu))    

  return(list(bstar2_mu=bstar2_mu,bstar_mu=bstar_mu,bstar_mu_eps=bstar_mu_eps,bstar_mu_half=bstar_mu_half,t2=t2,t4=t4,t5=t5,t6=t6,t7=t7,KL=KL,TV=TV,KL_out2b=KL_out2b))    
  
  
  }





KL_Norm_apprx_old<-function(y,x,P,Pstar,ystar,Ppost,mu_sim,family,link){
    
  # Get dimensions
  
  n<-length(mu_sim[,1])
  l1<-length(x[,1])
  
  eps_in<-matrix(1,nrow=l1,ncol=1) # figure out how to move outside of function
  
  
  # Setup output matrices
  
  bstar2_mu<-matrix(0,nrow=l1,ncol=n) # Matrix of conditional posterior modes (and means) of normal approximation for posterior density
  bstar_mu<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  bstar_mu_eps<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  bstar_mu_half<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  
  Pstar_mu<-matrix(0,nrow=l1,ncol=n)
  Pstar2_mu<-matrix(0,nrow=l1,ncol=n)
  Pstar2_mu_eps<-matrix(0,nrow=l1,ncol=n)
  Pstar_mu_half<-matrix(0,nrow=l1,ncol=n)
  
  
  
  P_eps<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  P_eps_half<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  
  t1<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  t2<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  t3<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  t4<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  
  t5<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  t6<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  
  t7<-matrix(0,nrow=l1,ncol=n)  # Matrix of conditional posteriod modes of posterior density
  
  # Get family functions
  
  
  
  myfamfunc<-glmbfamfunc(list(family=family,link=link))  # Make this flexible (linked to family and link)
  f7_temp<-myfamfunc$f7
  f3_temp<-myfamfunc$f3
  f2_temp<-myfamfunc$f2
  f1_temp<-myfamfunc$f1
  
  # Loop through iterations
  
  for(j in 1:n){
    
    # Calculate conditional means for b under normal approximation
    
    xmutemp<-t(t(alpha1))+x%*%mu_sim[j:j,]
    
    
    bstar2_mu[,j]<-(P*xmutemp+Pstar*ystar)/Ppost
    
    # use bstar2_mu[,j] as starting point for finding bstar_mu[,j]
    
    bstar_mu[,j]<-bstar2_mu[,j]
    
    bstar_mu_eps[,j]<-bstar2_mu[,j]
    bstar_mu_half[,j]<-bstar2_mu[,j]
    
    # Loop Through dimensions
    
    for(i in 1:l1){
      
      # Leverage newtown's method (currently five iterations) to find bstar_mu
      
      for(k in 1:5){  
        a<-f3_temp(bstar_mu[i,j],y[i],1,xmutemp[i],P)
        
        b<-P+f7_temp(bstar_mu[i,j],y[i],1,xmutemp[i],P)   # need to modify to use Precision determined by family and link
        
        bstar_mu[i,j]<-bstar_mu[i,j]-(a/b)
      }
      
      for(k in 1:5){  
        a<-(1-(eps_in[i]))*f3_temp(bstar_mu_eps[i,j],y[i],1,xmutemp[i],P)+((eps_in[i]))*(P%*%(bstar_mu_eps[i,j]-xmutemp[i])+Pstar[i]%*%(bstar_mu_eps[i,j]-ystar[i]))
        b<-P+(1-(eps_in[i]))*f7_temp(bstar_mu_eps[i,j],y[i],1,xmutemp[i],P)+(eps_in[i])*Pstar[i]   # need to modify to use Precision determined by family and link
        
        
        bstar_mu_eps[i,j]<-bstar_mu_eps[i,j]-(a/b)
      }
      
      
      for(k in 1:5){  
        a<-(1-(eps_in[i]/2))*f3_temp(bstar_mu_half[i,j],y[i],1,xmutemp[i],P)+((eps_in[i]/2))*(P%*%(bstar_mu_half[i,j]-xmutemp[i])+Pstar[i]%*%(bstar_mu_half[i,j]-ystar[i]))
        b<-P+(1-(eps_in[i]/2))*f7_temp(bstar_mu_half[i,j],y[i],1,xmutemp[i],P)+(eps_in[i]/2)*Pstar[i]   # need to modify to use Precision determined by family and link
        
        
        bstar_mu_half[i,j]<-bstar_mu_half[i,j]-(a/b)
      }
      
      
      Pstar_mu[i,j]<-f7_temp(bstar_mu[i,j],y[i],1,xmutemp[i],P) 
      Pstar2_mu[i,j]<-f7_temp(bstar2_mu[i,j],y[i],1,xmutemp[i],P) 
      Pstar2_mu_eps[i,j]<-f7_temp(bstar_mu_eps[i,j],y[i],1,xmutemp[i],P) 
      Pstar_mu_half[i,j]<-f7_temp(bstar_mu_half[i,j],y[i],1,xmutemp[i],P) 
      
      
    }
    
  }

  
  #################################################
    
  for(j in 1:n){
    for(i in 1:l1){
      
      P_eps[i,j]<-P+(1-eps_in[i])*Pstar2_mu_eps[i,j]+eps_in[i]*Pstar[i]
      P_eps_half[i,j]<-P+(1-(eps_in[i]/2))*Pstar_mu_half[i,j]+(eps_in[i]/2)*Pstar[i]
      
      t1[i,j]<-0.5*(1-eps_in[i])*(1/P_eps[i,j])*(Pstar_mu[i,j]-Pstar2_mu_eps[i,j])
      t2[i,j]<-0.5*(eps_in[i])*(1/P_eps[i,j])*(Pstar_mu[i,j]-Pstar[i])
      
      t3[i,j]<-0.5*(1-eps_in[i])*(1/P_eps_half[i,j])*(Pstar_mu[i,j]-Pstar2_mu_eps[i,j])
      t4[i,j]<-0.5*(eps_in[i])*(1/P_eps_half[i,j])*(Pstar_mu[i,j]-Pstar[i])
      
      t5[i,j]<-0.5*(bstar_mu_eps[i,j]-bstar_mu[i,j])*(P+Pstar_mu[i,j])*(bstar_mu_eps[i,j]-bstar_mu[i,j])
      t6[i,j]<-0.5*(bstar_mu_half[i,j]-bstar_mu[i,j])*(P+Pstar_mu[i,j])*(bstar_mu_half[i,j]-bstar_mu[i,j])
      
      t7[i,j]<-0.5*(bstar_mu_half[i,j]-bstar_mu_eps[i,j])*((1-eps_in[i])*(P+Pstar2_mu_eps[i,j])+(eps_in[i])*(P+Pstar[i]))*(bstar_mu_half[i,j]-bstar_mu_eps[i,j])
      
    }
  }
  
  
  
  KL_out1<-t(t1)+t(t2)-t(t3)-t(t4)+t(t5)-t(t6)+t(t7)
  KL_out2<-matrix(0,nrow=n,ncol=1)
  KL_out2b<-matrix(0,nrow=1,ncol=l1)
  
  for(i in 1:n){
    KL_out2[i]<-sum(KL_out1[i,1:l1])
  }
  
  for(j in 1:l1){
    KL_out2b[1,j]<-mean(KL_out1[,j])
  }
  
  
  KL<-mean(KL_out2)
  
  TV<-sqrt(0.5*KL)
  
  return(list(bstar2_mu=bstar2_mu,bstar_mu=bstar_mu,bstar_mu_eps=bstar_mu_eps,bstar_mu_half=bstar_mu_half,t2=t2,t4=t4,t5=t5,t6=t6,t7=t7,KL=KL,TV=TV,KL_out2b=KL_out2b))    
}


