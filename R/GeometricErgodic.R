set_nstar<-function(par_in,epsilon_converge,trace_constant,lambda_star,epsilon1){
  
  gammastar<-par_in
  
  alpha_out<-(1+gammastar)/(1+trace_constant+lambda_star*gammastar)
  U_out<-(1+trace_constant+2*lambda_star*gammastar)
  A1_out<-(1-exp(log(epsilon1)-gammastar))
    
  r_star1<-log(alpha_out)/(log(U_out)+log(alpha_out)-log(A1_out))
  
  A3<-A1_out^r_star1
  
  nstar<-((log(epsilon_converge))-log(2+gammastar))/log(A3)
  
  if(nstar>100000000){nstar<-100000000}
  if(nstar<0){nstar<-100000000}
  
  return(nstar)
  
  
}
