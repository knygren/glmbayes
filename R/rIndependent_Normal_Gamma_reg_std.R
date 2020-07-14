#' The Standard Bayesian Indendepent Normal-Gamma Regression Distribution
#'
#' Generates iid samples from the posterior density for the standard
#' independent normal-gamma regression
#' @param mu Prior mean
#' @param P Prior Precision
#' @param alpha offset
#' @param wt weights
#' @param link link function
#' @param progbar progress bar
#' @param f2 Gradient function
#' @param Envelope Envelope object
#' @param gamma_list list with parameters used by to generate candidates from a restricted gamma
#' @param UB_list list with constants used during accept-reject calculations.
#' @inheritParams rindependent_norm_gamma_reg
#' @return Currently mainly the draws for the dispersion and the regression coefficients
#' will be updated to return outputs consistent with other function
#' @example inst/examples/Ex_rindep_norm_gamma_reg.R
#' @export 


rindep_norm_gamma_reg_std_R <- function(n, y, x, mu, P, alpha, wt, f2, Envelope, 
                                        gamma_list, 
                                        UB_list,
                                        family, link, progbar = 1L) {
  
  ## Pull constants from Env2, gamma_list, and UB_list
  
  Env2=Envelope
  shape3=gamma_list$shape3
  rate2=gamma_list$rate2
  disp_upper=gamma_list$disp_upper
  disp_lower=gamma_list$disp_lower
  gs=nrow(Env2$cbars) # Size of grid
  New_LL=c(1:gs)
  
  RSS_ML=UB_list$RSS_ML
  max_New_LL_UB=UB_list$max_New_LL_UB
  max_LL_log_disp=UB_list$max_LL_log_disp
  lm_log1=UB_list$lm_log1
  lm_log2=UB_list$lm_log2 
  log_P_diff=UB_list$log_P_diff
  
  
  
  # Initialized Objects that will be populated during the loop
  
  iters_out<-rep(1,n)
  disp_out<-matrix(0,nrow=n,ncol=1)
  beta_out<-matrix(0,nrow=n,ncol=ncol(x))
  weight_out<-0*c(1:n)
  
  for(i in 1:n){  
    a1=0  
    
    while(a1==0){
      
      # Generate Candidate for dispersion from Restricted Gamma
      dispersion=r_invgamma(1,shape=shape3,rate=rate2,
                            disp_upper=disp_upper,disp_lower=disp_lower)
      #p=1/dispersion   # This may not be needed
      
      
      ## Setup for Candidate Generation from Restricted Multivariate Normal
      wt2=wt/rep(dispersion,length(y))
      New_thetabars=Inv_f3_gaussian(t(Env2$cbars), y, as.matrix(x),as.matrix(mu,ncol=1), as.matrix(P), 
                                    as.vector(alpha), as.vector(wt2))
      LL_New=-f2_gaussian_vector(t(New_thetabars), y, as.matrix(x), as.matrix(mu,ncol=1),
                                 as.matrix(P), as.vector(alpha), as.vector(wt2))
      
      
      # Generate Candidate for Restricted Multivariate Normal
      
      sim=.rindep_norm_gamma_reg_std_cpp(n=1,y=y,x=x,mu=mu,P=P,alpha=alpha,wt=wt2,
                                         f2=f2,Envelope=Env2,family="gaussian",link="identity",as.integer(0))
      
      # Get Key output from *.cpp function
      J_out=sim$J+1
      log_U2=sim$log_U2
      
      # Store Candidates in permananent array - should move towards end of function
      
      disp_out[i,1]=dispersion
      beta_out[i,1:ncol(x)]=sim$out[1,1:ncol(x)]
      
      # Calculate LL function for candidate
      
      LL_Test=-f2_gaussian_vector(as.matrix(beta_out[i,1:ncol(x)],ncol=1), y, as.matrix(x), 
                                  as.matrix(mu,ncol=1),as.matrix(P), as.vector(alpha), as.vector(wt2))
      
      ###Setup to calculate acceptance-rejection [This may not be needed]
      
      #ntheta_star=as.matrix(thetastars1[J_out,],ncol=1) --> May not be used here
      
      ## Series of Upper Bounds
      
      # Block 1: UB1 [Same form as in fixed dispersion case]
      betadiff=as.matrix(sim$out[1,1:ncol(x)],ncol=1)-as.matrix(New_thetabars[J_out,1:ncol(x)],ncol=1)
      UB1=LL_New[J_out]-Env2$cbars[J_out,1:ncol(x)]%*%(betadiff)
      
      ## Block 2: UB2 [Bounded by shifting RSS_ML to Gamma]
      
      ntheta=as.matrix(New_thetabars[J_out,1:ncol(x)],ncol=1)
      yxbeta=(y-alpha-x%*%ntheta)*sqrt(wt)
      RSS_ntheta=t(yxbeta)%*%(yxbeta)
      UB2  =  0.5*(1/dispersion)*(RSS_ntheta-RSS_ML)
      
      
      ## Block 3: UB3A [Shifts factors for components in grid as a function of New_thetabars and New_LL]
      ##          If the variance for the dispersion is low, this term should be small
      
      ##    UB1 contains a term that is quadratic in New_thetabars[J_out,1:ncol(x)]
      ##    When there is more than one Likelihood subgradient density, 
      ##    This term varies across the likelihood subgradient densities
      ##    and therefore needs to be handled across them
      
      for(j in 1:gs){
        
        theta_bars_temp=as.matrix(New_thetabars[j,1:ncol(x)],ncol=1)
        cbars_temp=as.matrix(Env2$cbars[j,1:ncol(x)],ncol=1)
        New_LL[j]=-0.5*t(theta_bars_temp)%*%P%*%theta_bars_temp+t(cbars_temp)%*% theta_bars_temp
        #logP(i,1)=logP(i,0)-NegLL(i)+0.5*arma::as_scalar(cbarrow.t() * cbarrow)+arma::as_scalar(G3row.t() * cbarrow);
      }
      
      LL_temp=log_P_diff+New_LL
      max_New_LL=max(LL_temp)
      UB3A =max_New_LL -LL_temp[J_out]
      
      ## Block UB3B - Adjusts for different max factors for different dispersions
      
      New_LL_log_disp=lm_log1+lm_log2*log(dispersion)
      UB3B=(max_New_LL_UB-max_LL_log_disp)-(max_New_LL-New_LL_log_disp)
      
      test1= LL_Test-UB1
      test2= LL_Test-(UB1+UB2)
      test3A= LL_Test-(UB1+UB2+UB3A)
      test3B= LL_Test-(UB1+UB2+UB3A+UB3B)
      
      
      test3B= LL_Test-(UB1+UB2+UB3A+UB3B)
      
      if(test3B>0){ 
        warning("Candidate with Positive Log-Acceptance. Bounding function may need adjustment")  
      }
      
      #    print("tests 1, 2, 3A, and 3B")
      #    print(test1)
      #    print(test2)
      #    print(test3A)
      #    print(test3B)
      
      test=test3B-log_U2
      
      #    print("test")
      #    print(test)
      
      #    stop("test 3B and test above")
      
      weight_out[i]=max_New_LL
      
      #test=1  # edit out later
      #a1=1
      if(test>=0) a1=1
      else{iters_out[i]=iters_out[i]+1}        
      
      
    }
    
    
  }
  
  return(list(beta_out=beta_out,disp_out=disp_out,iters_out=iters_out, weight_out=weight_out))
  
}

