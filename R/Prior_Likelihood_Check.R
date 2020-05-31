Prior_Likelihood_Check<-function(Prior_mean, Prior_std, Like_est, Like_std){
  
  Prior_diff=(Like_est-Prior_mean)/Prior_std
  Likelihood_diff=(Prior_mean-Like_est)/Like_std
  Abs_Prior_error=sqrt(Prior_diff^2)
  Abs_Likelihood_error=sqrt(Likelihood_diff^2)
  
  Error_Checks=cbind(Prior_diff,Likelihood_diff,Abs_Prior_error,Abs_Likelihood_error)
  
  colnames(Error_Checks)=c("Prior_std_diff","Likelihood_std_diff","Abs_Prior_std_error","Abs_Likelihood_std_error")
  
  Max_Prior_error=max(Abs_Prior_error)
  Max_Likelihood_error=max(Abs_Likelihood_error)
  
  if (Max_Prior_error>3){ 
    warning("At least one of the maximum likelihood estimates is more than 3 standard deviations 
            away from the prior point estimate (using the prior standard deviation). This suggest that the 
            prior point estimate may be poorly choosen and/or that the prior could be set too strong 
            (the prior standard deviation is too small).")}
  
  if (Max_Likelihood_error>3){ 
    warning("At least one of the prior point estimates is more than 3 standard deviations away from 
            the maximum likeihood estimate (using the likelihood standard deviation). This suggest that the 
            prior point estimate may be poorly choosen and inconsistent with the data.")}
  
    return(Error_Checks)
  }

