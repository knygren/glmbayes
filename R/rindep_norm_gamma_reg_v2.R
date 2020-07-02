#' The Bayesian Indendepent Normal-Gamma Regression Distribution
#'
#' Generates iid samples from the posterior density for the 
#' independent normal-gamma regression
#' @param prior_list a list with the prior parameters (mu, Sigma, shape and rate) for the 
#' prior distribution.
#' @param offset an offset parameter
#' @param weights a weighting variable
#' @inheritParams glmb
#' @return Currently mainly the draws for the dispersion and the regression coefficients
#' will be updated to return outputs consistent with other function
#' @family simfuncs 
#' @example inst/examples/Ex_confint.glmb.R
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

rindependent_norm_gamma_reg_v2<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian()){

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

  #rate2 =rate + RSS/2
  
  dispersion2=dispersion

  for(j in 1:10){
  prior=list(mu=mu,Sigma=Sigma, dispersion=dispersion2)
  glmb_out1=glmb(n=1,y~x-1,family=gaussian(),prior=prior)
  b_old=glmb_out1$coef.mode

  xbetastar=x%*%b_old
  
  ## Residual SUM of SQUARES at current conditional posterior mode estimate
  
  RSS2_post=t(y-xbetastar)%*%(y-xbetastar)  
  
  dispersion2=RSS2_post/n_obs
  
    }

  prior_list2=list(mu=mu,Sigma=Sigma, dispersion=dispersion2,RSS2_post=RSS2_post)
  glmb_out1=glmb(n=1,y~x-1,family=gaussian(),prior=prior_list2)
  
  ## Use this as betastar
  
  betastar=glmb_out1$coef.mode
  
  
  xbetastar=x%*%betastar
  RSS2_post=t(y-xbetastar)%*%(y-xbetastar)  ## Residual SUM of SQUARES at post model

  ## Updated version (when using the likelihood subgradient approach
  ## The RSS at the tangent point gets shifted to the gamma distribution
  
  shape2= shape + n_obs/2
  rate2 =rate + RSS2_post/2

#  print("dispersion from optimization")
#  print(dispersion2)
  
  dispstar=RSS2_post/n_obs

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
  #n=1000
  
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
  
  opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x2),mu=as.vector(mu2),
                P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
                method="BFGS",hessian=TRUE
  )
  
  bstar=opt_out$par  ## Posterior mode for adjusted model
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
  
  Gridtype=as.integer(4)
  
  ## Note: At the posterior mode only, thetabaras and cbars are just the opposite sign of each other
  
  Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
                     as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),
                     as.vector(wt2),
                     family="gaussian",link="identity",
                     Gridtype=Gridtype, n=as.integer(n),
                     sortgrid=TRUE)
 
  
  
  
  
#  print("Original_thetabars")
#  print(Env2$thetabars)
  
  
  ll_Check=f2_gaussian_vector(t(Env2$thetabars), y, as.matrix(x2), as.matrix(mu2,ncol=1),
                              as.matrix(P2), as.vector(alpha), as.vector(wt2))
  
  
  
  #print("NegLL - Original")
  #print(Env2$NegLL)

#  print("RSS2_post")
#  print(RSS2_post)
  
#  print("0.5*RSS2_post/dispstar")
#  print(0.5*RSS2_post/dispstar)
  
  #  Just the Log-Likelihood
  
  test_const1=-Env2$NegLL   
  
  ## Removing normalizing constant and shifting log p to gamma
  
  test_const2=-Env2$NegLL+(n_obs/2)*log(dispstar)-(n_obs/2)*log(2*pi)
  
  ## Not sure why Env2$RSS is present - list seemingly only returns RSS_Out [from lm it seems perhaps]
  
  #stop("This was the envelope")
  
  RSS=Env2$RSS*dispstar[1,1]
  RSS_Min=min(RSS)  ## This should likely be moved to the gamma distribution rate parameter
  RSS_Diff=RSS-RSS_Min   
  

  test_const3=test_const2+0.5*(1/dispstar[1,1])*RSS_Min
  
  # remaining differencees are now from the prior components and do no depend on dispersion
  
  test_const4=test_const3+0.5*(1/dispstar[1,1])*RSS_Diff  
  
  

  ## Set the updated shape and rate parameters - may or may not coincide with RSS2_post
  
  shape2= shape + n_obs/2
  rate2 =rate + RSS_Min/2

  ## These can be initialized here or even at the top of the function
  
  disp_out<-matrix(0,nrow=n,ncol=1)
  beta_out<-matrix(0,nrow=n,ncol=ncol(x))
  test_out<-matrix(0,nrow=n,ncol=3)
  iters_out<-matrix(0,nrow=n,ncol=1)
  
  # Calculate cbars_int - can actually be quite large (~50% of the likelihood cbars but with opposite sign)
  # In part because epsilon is currently not being brought up close to largest possible matrix in code
  #
  
  Env2$cbars_int=Env2$cbars-Env2$cbars_slope
  
#  print(Env2$cbars_int)
#  print(Env2$cbars_slope)  # includes precision multiplier
#  print(Env2$cbars)        # includes precision multiplier
  
#  stop("The above are the different cbars")  
  
  Env2$NegLL_int=Env2$NegLL-Env2$NegLL_slope

  # This should match the slope
  
  NegLL_temp_part2=0.5*RSS*(1/dispstar[1,1])+rep((n_obs/2)*log(2*pi),length(RSS))-rep((n_obs/2)*log(1/dispstar[1,1]),length(RSS) )
  
#   print(Env2$NegLL_int)
#   print(Env2$NegLL_slope)  # includes precision multiplier
#   print(NegLL_temp_part2)  # recomputed slope [should match previous row]
#   print(Env2$NegLL)        # includes precision multiplier
  
  NegLL_int=Env2$NegLL_int
  
#  stop("The above are the different NegLL")  
  

  ######  Call Set_Grid function with   
  
  GIndex=Env2$GridIndex
  #cbars_temp=Env2$cbars
  Lint=Env2$Lint1
  
  ## This is the confusing part - need to understand how to update this 
  ## May not matter when using a single Likelihood subgradient density
  ## Likely does not make much sense to call this function for intercept only
  ## or for just the "slope" of cbars
  ## Could possbly make sense to update these for different dispersion values as slopes change
  ## For single likelihood subgradient density, logP=0 regardless so not relevant
  
  outgrid=.Set_Grid_cpp(GIndex, Env2$cbars, Lint)
#  outgrid_int=.Set_Grid_cpp(GIndex, Env2$cbars_int, Lint)

  
  logP=outgrid$logP

#  print(logP)

#    stop("The above are the different logP values")  

  ### Call the setlogP function to compute the probability with which
  ### Each component in the grid should be visited
  
  outlogP=.setlogP_cpp(outgrid$logP, Env2$NegLL,  Env2$cbars,Env2$thetabars)
#  outlogP_int=.setlogP_cpp(logP_int, NegLL_int, Env2$cbars_int, Env2$thetabars)

  ## Probably should adjust this
    
  Env2$LLconst_int=outlogP$LLconst
  

  logP2=outlogP$logP
  #logP2_int=outlogP_int$logP
  
  LLconst=outlogP$LLconst
  logP2_2 = logP2[,2]
  
  maxlogP=max(logP2_2)
  
  PLSD=exp(logP2_2-maxlogP)
  PLSD=PLSD/sum(PLSD)
  

  Env_temp=Env2

#  print("Difference from old Neg_LL")
  
## Note: If the mean of the likelihood subgradient densities are to be kept constant, then one likely needs the 
##       following process
  
##      1) Find the initial means (equals - cbars and not thetabars!)
##      2) For each simulated dispersion find the cbar corresponding to the original mean and
##      3) then back into the implied thetabar

########################  End of test ############################

  ## Because of Acceptance procedure below, this should use RSS_ML not RSS_Post from above
  
  shape2= shape + n_obs/2
  rate2 =rate + RSS_ML/2
  
  ## Loop through and accept/reject based on test from internal function

  for(i in 1:n){  
    a1=0  
    
    iters_out[i]=1
    
    while(a1==0){
      p=rgamma(1,shape=shape2,rate=rate2)  
      dispersion=1/p
      
      ## Update wt2

      wt2=wt/rep(dispersion,length(y))

      # calculate new_thetabars for new dispersion - should be similar
      # but slightly different from original thetabars
      
      New_thetabars=Inv_f3_gaussian(t(Env2$cbars), y, as.matrix(x2),as.matrix(mu2,ncol=1), as.matrix(P2), as.vector(alpha), as.vector(wt2))
      Env_temp$thetabars=New_thetabars
      
      ## This may have been the problem (was using the old thetabars) 
      
      NegLL_New=f2_gaussian_vector(t(Env_temp$thetabars), y, as.matrix(x2), as.matrix(mu2,ncol=1),
                                   as.matrix(P2), as.vector(alpha), as.vector(wt2))
      



      ## In single likelihood subgradient case, we likely don't need to call the Set_Grid 
      ## and Set_logP function but can turn directly to the simulation
      
      ## The mean for the restricted normal is determined to cbars. When the precision (dispersion)
      ## changes, then question is if this should change...
      ## The points of tangencies should not change, but the slope likely should, 
      ## Hence we should use Env_temp here
      ## Changing cbars each time should increase variance
      ## This should likely be changed to call an updated envelope
      
      
      sim=.rindep_norm_gamma_reg_std_cpp(n=1,y=y,x=x2,mu=mu2,P=P2,alpha=alpha,wt=wt2,
                         f2=f2,Envelope=Env2,family="gaussian",link="identity",as.integer(0))
      
      
  
      
      ## Add 1 here becasuse arrays in *.cpp start with 0 instead of 1

      J_out=sim$J+1
      log_U2=sim$log_U2
      
      disp_out[i,1]=dispersion
      beta_out[i,1:ncol(x)]=sim$out[1,1:ncol(x)]
    
      LL_Test=-f2_gaussian_vector(as.matrix(beta_out[i,1:ncol(x)],ncol=1), y, as.matrix(x2), as.matrix(mu2,ncol=1),
                                  as.matrix(P2), as.vector(alpha), as.vector(wt2))
      

      # RSS got replaced somewhere along the way above so use RSS_ML 
      # UB1 uses only terms that can be shifted to the Gamma distribution and the constant with pi

      f3=famfunc$f3  ## These should match original - can run check
      
      cbars_new=f3(as.matrix(New_thetabars[J_out,],ncol=1), y, as.matrix(x2), as.matrix(mu2,ncol=1),
                                  as.matrix(P2), as.vector(alpha), as.vector(wt2))
      

      betadiff=as.matrix(sim$out[1,1:ncol(x)],ncol=1)-as.matrix(New_thetabars[J_out,1:ncol(x)],ncol=1)
      
      ## UB1 --> Should be test if we are using prior for beta but shift the two components to gamm
      ## UB2 --> Should be test if we were sampling for the conditional density using the grid positioned at betastar, dispstar
      ## UB3 --> Should be test if we are sampling from the density with different variances but keeping mean of
      ##         normal densities fixed 
      ##         unclear if/how the probabilities of visiting different parts of the grid should be adjusted based on what the dipserson
      ##         is. Next step is to figure this part out
      
      UB1=-0.5*(1/dispersion)*RSS_ML-(n_obs/2)*log(dispersion)-(n_obs/2)*log(2*pi)  
      UB2=-NegLL_New[J_out]-t(cbars_new)%*%(betadiff)
      UB3=UB1-t(cbars_new)%*%(betadiff)
      
      Diff1=LL_Test-UB1
      Diff2=LL_Test-UB2
      Diff3=LL_Test-UB3
      
      test=Diff3-log_U2
      

      # a1=1
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
test_out=test_out)

    colnames(outlist$coefficients)<-colnames(x)
    outlist$offset2<-offset2
    class(outlist)<-c(outlist$class,"rglmb")
    
return(outlist)  

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



