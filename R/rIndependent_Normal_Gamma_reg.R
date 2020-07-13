#' The Bayesian Indendepent Normal-Gamma Regression Distribution
#'
#' Generates iid samples from the posterior density for the 
#' independent normal-gamma regression
#' @param prior_list a list with the prior parameters (mu, Sigma, shape and rate) for the 
#' prior distribution.
#' @param offset an offset parameter
#' @param weights a weighting variable
#' @param max_disp parameter currently used to control upper bound in accept-reject procedure 
#' @inheritParams glmb
#' @return Currently mainly the draws for the dispersion and the regression coefficients
#' will be updated to return outputs consistent with other function
#' @family simfuncs 
#' @example inst/examples/Ex_rindep_norm_gamma_reg.R
#' @export 


## Single Likelihood subgradient model below seems like it may be working but slow
## To move towards a Grid based approach, this function goes part way 
## by still using a single point but using a likelihood function that includes 
## a shifted part from the prio


##########################  Grid Approach ################################
##
##   1) Use betastar and pstar from single tangent as starting point
##   2) Find Envelope for model where p is fixed at pstar 
##   3) Find all the points of tangencies and compute RSS
##   4) Find the smallest of all the RSS estimates (RSS_LB) and set
##      pstar_UB = 1/(RSS_LB)
##   5) Adjust the candidates for p to use RSS_LB instead of RSS_post to
##      determine the shape
##   6) When implementing accept-reject sampling using the grid, adjust the accept-reject
##      by adding the test1 term back in with value test1=-0.5p(RSS_post-RSS_LB)
##    Note 1: This should ensure that the accept-reject approach still is valie
##    Note 2: Need to compare single point to grid to ensure validity


rindependent_norm_gamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                                         max_disp=0.9){
  
  call<-match.call()
  
  offset2=offset
  wt=weights
  
  ### Initial implementation of Likelihood subgradient Sampling 
  ### Currently uses as single point for conditional tangencis
  ### (at conditional posterior modes)
  ### Verify this yields correct results and then try to implement grid approach
  
  ## Use the prior list to set the prior elements if it is not missing
  ## Error checking to verify that the correct elements are present
  ## Shold be implemented
  
  if(missing(prior_list)) stop("Prior Specification Missing")
  if(!missing(prior_list)){
    if(!is.null(prior_list$mu)) mu=prior_list$mu
    if(!is.null(prior_list$Sigma)) Sigma=prior_list$Sigma
    if(!is.null(prior_list$dispersion)) dispersion=prior_list$dispersion
    else dispersion=NULL
    if(!is.null(prior_list$shape)) shape=prior_list$shape
    else shape=NULL
    if(!is.null(prior_list$rate)) rate=prior_list$rate
    else rate=NULL
  }
  
  lm_out=lm(y ~ x-1) # run classical regression to get maximum likelhood estimate
  RSS=sum(residuals(lm_out)^2)
  
  RSS_ML=sum(residuals(lm_out)^2)
  
  n_obs=length(y)
  
  
  dispersion2=dispersion
  
  for(j in 1:10){
    prior=list(mu=mu,Sigma=Sigma, dispersion=dispersion2)
    glmb_out1=glmb(y~x-1,family=gaussian(),dNormal(mu=mu,Sigma=Sigma,dispersion=dispersion))
    b_old=glmb_out1$coef.mode
    
    xbetastar=x%*%b_old
    
    ## Residual SUM of SQUARES at current conditional posterior mode estimate
    
    RSS2_post=t(y-xbetastar)%*%(y-xbetastar)  

    shape2= shape + n_obs/2
    rate2 =rate + RSS2_post/2
    
        
    # Inverse of implied expectation for precision used here
    # Seems consistent with first order Taylor series approximation for mean of betas 
    
    dispersion2=rate2/shape2
    
  }
  
  # Temporarily make dispersion2 very large
  

  prior_list2=list(mu=mu,Sigma=Sigma, dispersion=dispersion2,RSS2_post=RSS2_post)
  glmb_out1=glmb(y~x-1,family=gaussian(),dNormal(mu=mu,Sigma=Sigma,dispersion=dispersion2))
  
  #  glmb_out1=glmb(n=1,y~x-1,family=gaussian(),prior=prior_list2)
  
  ## Temporarily use maximum likelihood estimates as betastar
  
  betastar=glmb_out1$coef.mode
  
  
  xbetastar=x%*%betastar
  RSS2_post=t(y-xbetastar)%*%(y-xbetastar)  ## Residual SUM of SQUARES at post model
  
  ## Updated version (when using the likelihood subgradient approach
  ## The RSS at the tangent point gets shifted to the gamma distribution
  
  shape2= shape + n_obs/2
  rate2 =rate + RSS2_post/2
  
  #  print("dispersion from optimization")
  #  print(dispersion2)
  
  #gm<-gamma_median1(shape2,rate2)
  #print("dispersion Median estimate")
  #print(gm$v3)
  
  #dispstar=RSS2_post/n_obs
  #dispstar=rate2/shape2
  dispstar=rate2/(shape2-1)
  #dispstar=0.6452614  # Hard code to result of calculated further below
  
  
  ## Get family functions for gaussian()  
  
  famfunc<-glmbfamfunc(gaussian())
  
  f1<-famfunc$f1
  f2<-famfunc$f2
  f3<-famfunc$f3
  f5<-famfunc$f5
  f6<-famfunc$f6
  
  start <- mu
  
  if(is.null(offset2))  offset2=rep(as.numeric(0.0),length(y))
  P=solve(Sigma)
  
  ###### Adjust weight for dispersion
  
  dispersion2=dispstar
  
  if(is.null(wt))  wt=rep(1,length(y))
  if(length(wt)==1)  wt=rep(wt,length(y))
  
  wt2=wt/rep(dispersion2,length(y))
  
  ######################### Shift mean vector to offset so that adjusted model has 0 mean
  
  alpha=x%*%as.vector(mu)+offset2
  mu2=0*as.vector(mu)
  P2=P
  x2=x
  
  parin=start-mu
  #parin=start

  #parin=glmb_out1$glm$coefficients
  
  # This step only used to get posterior precision it seems
  # Since normal case, can likely be computed instead
  
  opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x2),mu=as.vector(mu2),
                P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
                method="BFGS",hessian=TRUE
  )

  bstar=opt_out$par  ## Posterior mode for adjusted model
  
  ## Temporarily use bstar as posterior mode
  
  
  #  bstar
  #  bstar+as.vector(mu)  # mode for actual model
  A1=opt_out$hessian # Approximate Precision at mode
  
  
  Standard_Mod=glmb_Standardize_Model(y=as.vector(y), x=as.matrix(x2),
                                      P=as.matrix(P2),bstar=as.matrix(bstar,ncol=1), A1=as.matrix(A1))
  
  bstar2=Standard_Mod$bstar2  
  A=Standard_Mod$A
  x2=Standard_Mod$x2
  mu2=Standard_Mod$mu2
  P2=Standard_Mod$P2
  L2Inv=Standard_Mod$L2Inv
  L3Inv=Standard_Mod$L3Inv
  
  
  ## Note, use Gridtype =4 here temporarily (Single Likelihood subgradient)
  
  Gridtype=as.integer(3)
  
  ## Pull the initial Envelope based on optimized values above
  
  Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
                     as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),
                     as.vector(wt2),
                     family="gaussian",link="identity",
                     Gridtype=Gridtype, n=as.integer(n),
                     sortgrid=TRUE)
  
  
  #### End of Standardization - Steps below might change  
  print("RSS - Post and ML")
  print(RSS2_post)
  print(RSS_ML)
  
  ## Update Gamma parameters
  
  shape2= shape + n_obs/2
  rate2 =rate + (RSS_ML/2)
  rate3 =rate + (RSS2_post/2)

#  lm_log2=5.555805   # Also used above
  
#  shape3=shape2-lm_log2
  
  print("Old Candidate shape and rate are")
  print(shape2)
  print(rate2)
  
  
  print("shape and rate for approximate posterior are")
  print(shape2)
  print(rate3)
  
  
  print("99st percentile for Inverse Gamma candidates density is")
  print(1/qgamma(c(0.01),shape2,rate2))
  max_disp1=1/qgamma(c(0.01),shape2,rate2)
  

  print("99th percentile for Inverse Gamma prior density is")
  print(1/qgamma(c(0.01),shape,rate))
  max_disp2=1/qgamma(c(0.01),shape,rate)
  
  print("99th percentile for approximate Inverse  Gamma posterior density is")
  
  upp=1/qgamma(c(0.01),shape2,rate3)
  low=1/qgamma(c(0.99),shape2,rate3)
  
  print(upp)
  
  
  print("1st percentile for approximate Inverse  Gamma posterior density is")
  #print(1/qgamma(c(0.99),shape2,rate3))
  print(low)
  
  a1=2
  b1=(upp-low)
  c1=-log(upp/low)
  
  dispstar_new= (-b1+ sqrt(b1^2-4*a1*c1))/(2*a1)

  print("dispstar_old")
  print(dispstar)
  
  dispstar=dispstar_new
  print("dispstar_new")
  print(dispstar_new)
  
    
  max_disp3=1/qgamma(c(0.01),shape2,rate3)
  
 
  print("max_disp3")
  print(max_disp3)
  
  print("quantile function lower bound for 1/dispersion under approximate posterior")
  pgamma(1/max_disp3,shape2,rate3)
  
  ## Copy Envelope for use in simulation (as needed update components)
  
  Env_temp=Env2
  
  ## get cbars
  
  cbars=Env2$cbars
  
  ## Calculate theta_star [Inv_f3_gaussian suggest positive sign for thetastar]
  ##
  theta_star=t(as.matrix(solve(P2)%*%as.matrix(t(cbars),ncol=1),ncol=1)) 
  
  gs=nrow(theta_star)
  
  
  logP1=Env2$logP
  New_LL=c(1:gs)
  New_logP2=c(1:gs)
  
#  print("logP - Old Envelope")
#  print(Env2$logP)
  
  for(i in 1:gs){
    
    cbars_temp=as.matrix(cbars[i,1:ncol(x)],ncol=1)
    
    #New_logP2[i]=logP1[i]+New_LL[i]+0.5*t(cbars_temp)%*%cbars_temp
    New_logP2[i]=logP1[i]+0.5*t(cbars_temp)%*%cbars_temp
    
    #logP(i,1)=logP(i,0)-NegLL(i)+0.5*arma::as_scalar(cbarrow.t() * cbarrow)+arma::as_scalar(G3row.t() * cbarrow);
  }
  
  
  
  New_thetabars=Inv_f3_gaussian(t(Env2$cbars), y, as.matrix(x2),as.matrix(mu2,ncol=1), as.matrix(P2), as.vector(alpha), as.vector(wt2))
  
#  print("New_thetabars that should match")
#  print(New_thetabars)

  #max_disp2=max_disp/2
  #max_disp_f=max(max_disp1)
  
  # max_disp3 is the approximation to posterior- makes more sense
  # than the prior or the candidata distribution as starting point
  # for now, link this to the 99th percentile for dispersion under approximate posterior
  
  
  max_disp_f=max_disp3
  wt2=wt/rep(max_disp_f,length(y))
    
  thetastars1=Inv_f3_gaussian(t(Env2$cbars), y, as.matrix(x2),as.matrix(mu2,ncol=1), as.matrix(P2), 
                              as.vector(alpha), as.vector(wt2))
  

  maxlogP=max(New_logP2)
  
  PLSD_new=exp(New_logP2-maxlogP)
  
  sumP=sum(PLSD_new)
  
  PLSD_new=PLSD_new/sumP
  
  Env_temp$PLSD=PLSD_new

  log_P_diff=log(Env_temp$PLSD)-log(Env2$PLSD)
  
  #################  this Block sets overall upper bound for LL_temp   #################
  
## Upper boundvalue for LL_temp using thetastars (based on upper bound dispersion)  
  
  for(j in 1:gs){
    
    theta_bars_temp=as.matrix(thetastars1[j,1:ncol(x)],ncol=1)
    cbars_temp=as.matrix(cbars[j,1:ncol(x)],ncol=1)
    
    # New formula - This should likely replace: -NegLL(i)+arma::as_scalar(G3row.t() * cbarrow)
    
    New_LL[j]=-0.5*t(theta_bars_temp)%*%P2%*%theta_bars_temp+t(cbars_temp)%*% theta_bars_temp
    
  }
  
  
  LL_temp=log_P_diff+New_LL
  
  
  max_New_LL_UB=max(LL_temp)
  
  print("max_New_LL_UB")  
  print(max_New_LL_UB)  
  
  #######################################  Use this block to find linear function
  
  
  
  dep_out<-matrix(0,nrow=100,ncol=1)
  disp_temp_out<-matrix(0,nrow=100,ncol=1)
    
  for(k in 1:100){
  
  # Generate some values based on unrestricted candidate density
  
  p=rgamma(1,shape=shape2,rate=rate2)  
  dispersion=1/p
  
  disp_temp_out[k]=dispersion
  
  wt2=wt/rep(dispersion,length(y))
  
  New_thetabars=Inv_f3_gaussian(t(Env2$cbars), y, as.matrix(x2),as.matrix(mu2,ncol=1), 
                                as.matrix(P2), as.vector(alpha), as.vector(wt2))
  
  LL_New=-f2_gaussian_vector(t(Env_temp$thetabars), y, as.matrix(x2), as.matrix(mu2,ncol=1),
                             as.matrix(P2), as.vector(alpha), as.vector(wt2))
  
  
  
  for(j in 1:gs){
    
    theta_bars_temp=as.matrix(New_thetabars[j,1:ncol(x)],ncol=1)
    cbars_temp=as.matrix(cbars[j,1:ncol(x)],ncol=1)
    
    # New formula - This should likely replace: -NegLL(i)+arma::as_scalar(G3row.t() * cbarrow)
    
    New_LL[j]=-0.5*t(theta_bars_temp)%*%P2%*%theta_bars_temp+t(cbars_temp)%*% theta_bars_temp
    
    #logP(i,1)=logP(i,0)-NegLL(i)+0.5*arma::as_scalar(cbarrow.t() * cbarrow)+arma::as_scalar(G3row.t() * cbarrow);
  }
  
  
  LL_temp=log_P_diff+New_LL
  max_New_LL=max(LL_temp)
  
  dep_out[k]=max_New_LL
  
  }
  
  lm_test2=lm(dep_out~disp_temp_out)
  lmc=lm_test2$coefficients

  print("lmc output - check if matches outside of function")  
  print(lmc)
  
  print("dispstar")
  print(dispstar)
  #dispstar=0.61
  lm_log2=lmc[2]*dispstar
  lm_log1=lmc[1]+lm_log2-lm_log2*log(dispstar)
  

  print("lm_log")
  
  print(lm_log1)
  print(lm_log2)
  
    
  shape3=shape2-lm_log2

  print("New Candidate shape and rate are")
  print(shape3)
  print(rate2)
  

  max_LL_log_disp=lm_log1+lm_log2*log(max_disp_f) ## From above
  
  
  
  
  sim_temp=rindep_norm_gamma_reg_std_R(n=n,y=y,x=x2,mu=mu2,P=P2,alpha=alpha,wt=wt,
  f2=f2,Envelope=Env2,
  gamma_list=list(shape3=shape3,rate2=rate2,disp_upper=upp,disp_lower=low),
  RSS_ML=RSS_ML,max_New_LL_UB=max_New_LL_UB,
  max_LL_log_disp=max_LL_log_disp,lm_log1=lm_log1,lm_log2=lm_log2, 
  log_P_diff=log_P_diff,
  family="gaussian",link="identity",as.integer(0))

#########################################  Post Processing for simulation function handling loop
  
  beta_out=sim_temp$beta_out
  disp_out=sim_temp$disp_out
  iters_out=sim_temp$iters_out
  weight_out=sim_temp$weight_out
  

  out=L2Inv%*%L3Inv%*%t(beta_out)

  for(i in 1:n){
    out[,i]=out[,i]+mu
  }
  
  
  famfunc=glmbfamfunc(gaussian())  
  f1=famfunc$f1
  
  outlist=list(
    coefficients=t(out), 
    coef.mode=betastar,  ## For now, use the conditional mode (not universal)
    dispersion=disp_out,
    ## For now, name items in list like this-eventually make format/names
    ## consistent with true prior (current names needed by summary function)
    Prior=list(mean=mu,Sigma=Sigma,shape=shape,rate=rate,Precision=solve(Sigma)), 
    family=gaussian(),
    prior.weights=wt,
    y=y,
    x=x,
    call=call,
    famfunc=famfunc,
    iters=iters_out,
    Envelope=NULL,
    loglike=NULL,
    weight_out=weight_out
    #,test_out=test_out
  )
  
  colnames(outlist$coefficients)<-colnames(x)
  outlist$offset2<-offset2
  class(outlist)<-c(outlist$class,"rglmb")
  
  return(outlist)  
  
  ########################################################################################3
  
  #################################################################################
  
  ########################  End of test ############################
  
  
  iters_out<-c(1:n)
  disp_out<-matrix(0,nrow=n,ncol=1)
  beta_out<-matrix(0,nrow=n,ncol=ncol(x))
  weight_out<-0*c(1:n)
  
  ## Loop through and accept/reject based on test from internal function
  
  for(i in 1:n){  
    a1=0  
    
    iters_out[i]=1
    
    while(a1==0){
      
      #p=rgamma(1,shape=shape2,rate=rate2)  
      dispersion=r_invgamma(1,shape=shape3,rate=rate2,disp_upper=upp,disp_lower=low)
      p=1/dispersion
      wt2=wt/rep(dispersion,length(y))
      New_thetabars=Inv_f3_gaussian(t(Env2$cbars), y, as.matrix(x2),as.matrix(mu2,ncol=1), as.matrix(P2), 
      as.vector(alpha), as.vector(wt2))
      LL_New=-f2_gaussian_vector(t(New_thetabars), y, as.matrix(x2), as.matrix(mu2,ncol=1),
                                 as.matrix(P2), as.vector(alpha), as.vector(wt2))
      
      

      sim=.rindep_norm_gamma_reg_std_cpp(n=1,y=y,x=x2,mu=mu2,P=P2,alpha=alpha,wt=wt2,
      f2=f2,Envelope=Env2,family="gaussian",link="identity",as.integer(0))
      
      J_out=sim$J+1
      log_U2=sim$log_U2
      disp_out[i,1]=dispersion
      beta_out[i,1:ncol(x)]=sim$out[1,1:ncol(x)]
      
      LL_Test=-f2_gaussian_vector(as.matrix(beta_out[i,1:ncol(x)],ncol=1), y, as.matrix(x2), as.matrix(mu2,ncol=1),
                                  as.matrix(P2), as.vector(alpha), as.vector(wt2))
      
      ntheta=as.matrix(New_thetabars[J_out,1:ncol(x)],ncol=1)
      ntheta_star=as.matrix(thetastars1[J_out,],ncol=1)
      betadiff=as.matrix(sim$out[1,1:ncol(x)],ncol=1)-as.matrix(New_thetabars[J_out,1:ncol(x)],ncol=1)

      ## Series of Upper Bounds
  
      # Block 1: UB1
      
      UB1=LL_New[J_out]-Env2$cbars[J_out,1:ncol(x)]%*%(betadiff)

      ## UB1 contains a term that involves -0.5*(1/dispersion) * RSS(ntheta)
      ## This term is re-entered here and bounded using term also shifted to gamma candidate
      
      ## Block 2: UB2
      
      yxbeta=(y-alpha-x2%*%ntheta)*sqrt(wt)
      RSS_ntheta=t(yxbeta)%*%(yxbeta)
      UB2  =  0.5*(1/dispersion)*(RSS_ntheta-RSS_ML)
      
      ## Block 3: UB3A
      
      ##    UB1 contains a term that is quadratic in New_thetabars[J_out,1:ncol(x)]
      ##    When there is more than one Likelihood subgradient density, 
      ##    This term varies across the likelihood subgradient densities
      ##    and therefore needs to be handled across them

      for(j in 1:gs){
        
        theta_bars_temp=as.matrix(New_thetabars[j,1:ncol(x)],ncol=1)
        cbars_temp=as.matrix(cbars[j,1:ncol(x)],ncol=1)
        New_LL[j]=-0.5*t(theta_bars_temp)%*%P2%*%theta_bars_temp+t(cbars_temp)%*% theta_bars_temp
        #logP(i,1)=logP(i,0)-NegLL(i)+0.5*arma::as_scalar(cbarrow.t() * cbarrow)+arma::as_scalar(G3row.t() * cbarrow);
      }
    
      LL_temp=log_P_diff+New_LL
      max_New_LL=max(LL_temp)
      UB3A =max_New_LL -LL_temp[J_out]
      
      ## max_new_LL depends on the dispersion and in turn needs to be bounded
      
      ## max_new_LL appears to be linear in the dispersion
      ## Let's consider a function New_LL_log_disp (linear in log_dispersion)
      ## that is tangent to max_new_LL near center of distribution
      ## and bounds max_new_LL from below.
      ## such a log term can be used to modify density that generates
      ## dispersion candidates (by changing shape parameter)
      ## and can then be included here. then we instead try to bound max_new_LL - New_LL_log_disp
      ## This can not be achieved either, but we can use a reasonable upper bound for dispersion
      ## so that difference is bounded below that upper bound
      ## As a test, we hard code this function here (will then generalize)
      
      #lm_log1=6.637078
      #lm_log2=5.555805   # Also used above
      
      
      #New_LL_log_disp=0
      New_LL_log_disp=lm_log1+lm_log2*log(dispersion)
      max_LL_log_disp=lm_log1+lm_log2*log(max_disp_f) ## From above
      
      
      UB3B=max_New_LL_UB-max_New_LL
      
      #print("UB3B - Old Method")
      #print(UB3B)
      UB3B=(max_New_LL_UB-max_LL_log_disp)-(max_New_LL-New_LL_log_disp)
      
      if(UB3B<0){
        
        # This does fail for some values close to lower bound
        # could be small inacurracies - adjust sampling
        # to also use lower bound
        #print("UB3B was negative with dispersion:")
        #print("Dispersion for Unexpected UB3B")
        #print(dispersion)
      }

      test1= LL_Test-UB1
      test2= LL_Test-(UB1+UB2)
      test3A= LL_Test-(UB1+UB2+UB3A)
      test3B= LL_Test-(UB1+UB2+UB3A+UB3B)

      if(test3B>0){ 
        warning("Candidate with Positive Log-Acceptance. Bounding function may need adjustment")  
        #cat("\nAssociated value for UB3B was:",UB3B)
        #cat("\nAssociated value for candidate dispersion was:",dispersion)
        #cat("\nOverall Log-Accepance rate:",test3B)
      }
      
      
      #test4= LL_Test-(UB1+UB2+UB3A+UB4)
      
      ## For now use test3A here (will likely need to further bound across dispersion values)
      ## but this could be close
      
      test=test3A-log_U2
      test=test3B-log_U2
      
      weight_out[i]=max_New_LL
          
      #a1=1
      if(test>=0) a1=1
      else{iters_out[i]=iters_out[i]+1}        
      
    }        
  }
  
  
  out=L2Inv%*%L3Inv%*%t(beta_out)
  
  for(i in 1:n){
    out[,i]=out[,i]+mu
  }
  
  
  famfunc=glmbfamfunc(gaussian())  
  f1=famfunc$f1
  
  outlist=list(
    coefficients=t(out), 
    coef.mode=betastar,  ## For now, use the conditional mode (not universal)
    dispersion=disp_out,
    ## For now, name items in list like this-eventually make format/names
    ## consistent with true prior (current names needed by summary function)
    Prior=list(mean=mu,Sigma=Sigma,shape=shape,rate=rate,Precision=solve(Sigma)), 
    family=gaussian(),
    prior.weights=wt,
    y=y,
    x=x,
    call=call,
    famfunc=famfunc,
    iters=iters_out,
    Envelope=NULL,
    loglike=NULL,
    weight_out=weight_out
    #,test_out=test_out
  )
  
  colnames(outlist$coefficients)<-colnames(x)
  outlist$offset2<-offset2
  class(outlist)<-c(outlist$class,"rglmb")
  
  return(outlist)  
  
}



rindep_norm_gamma_reg_std_R <- function(n, y, x, mu, P, alpha, wt, f2, Envelope, 
  gamma_list, RSS_ML,max_New_LL_UB,max_LL_log_disp,
  lm_log1,lm_log2, log_P_diff,
  family, link, progbar = 1L) {

  ## Pull constants from Env2 and gamma_list
  
  Env2=Envelope
  shape3=gamma_list$shape3
  rate2=gamma_list$rate2
  disp_upper=gamma_list$disp_upper
  disp_lower=gamma_list$disp_lower
  gs=nrow(Env2$cbars) # Size of grid
  New_LL=c(1:gs)
  
  # Initialized Objects that will be populated during the loop
  
  iters_out<-rep(1,length(y))
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
  


## This function is used by the above (not sure why Neg_logLik is not working)
## Could be because it is not exported - replace with Neg_logLik
#' @rdname Neg_logLik
#' @export 

Neg_logLik2<-function(b, y, x, alpha, wt,family){
  
  ## Add required checks on other inputs at the top
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  okfamilies <- c("gaussian","poisson","binomial","quasipoisson","quasibinomial","Gamma")
  if(family$family %in% okfamilies){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-c("log")		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-c("logit","probit","cloglog")		
    if(family$family=="Gamma") oklinks<-c("log")		
    if(family$link %in% oklinks){
      
      ## This may be the R version of these files so may not be using the efficiency of C++
      ## This may be safer
      
      famfunc<-glmbfamfunc(family)
      f1<-famfunc$f1
      f2<-famfunc$f2
      f3<-famfunc$f3
      #      f5<-famfunc$f5
      #      f6<-famfunc$f6
    }
    else{
      stop(gettextf("link \"%s\" not available for selected family; available links are %s", 
                    family$link , paste(sQuote(oklinks), collapse = ", ")), 
           domain = NA)
      
    }	
    
  }		
  else {
    stop(gettextf("family \"%s\" not available in glmb; available families are %s", 
                  family$family , paste(sQuote(okfamilies), collapse = ", ")), 
         domain = NA)
  }
  
  return(f1(b, y, x, alpha, wt)) 
  
}



################################## Utility functions used by the above  #################

p_inv_gamma<-function(dispersion,shape,rate){
  1-pgamma(1/dispersion,shape=shape,rate=rate)
}

q_inv_gamma<-function(p,shape,rate,disp_upper,disp_lower){
  p_upp=p_inv_gamma(disp_upper,shape=shape,rate=rate)
  p_low=p_inv_gamma(disp_lower,shape=shape,rate=rate)
  p1=p_low+p*(p_upp-p_low)
  p2=1-p1
  1/qgamma(p2,shape,rate)
}

r_invgamma<-function(n,shape,rate,disp_upper,disp_lower){
  p=runif(n)
  q_inv_gamma(p=p,shape=shape,rate=rate,disp_upper=disp_upper,disp_lower)
}



# From gamma distribution Wikipedia page

gamma_median1<-function(k,rate){
  
  v1=k-(1/3) +(8/(405*k)) +184/(25515*k*k)+2248/(3444525*k*k*k)-19006408/(15345358875*k*k*k*k) 
  
  v2=1/v1
  
  v3=v2*rate
  
  return(list(v1=v1,v2=v2,v3=v3))
}
