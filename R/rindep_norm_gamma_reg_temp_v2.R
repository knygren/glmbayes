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
## To use Grid based approach for this model, try the Below
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

rindependent_norm_gamma_reg_temp_v2<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian()){

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
  dispersion2=dispersion
  
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
  
  
  Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
                     as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),
                     as.vector(wt2),
                     family="gaussian",link="identity",
                     Gridtype=as.integer(3), n=as.integer(n),
                     sortgrid=TRUE)
 
  ## Return list needed by standard simulation function
  ## Standard simulation function may need to be modified to generate samples
  ## for the gaussian

  RSS=Env2$RSS*dispstar[1,1]
  RSS_Min=min(RSS)  ## This should likely be moved to the gamma distribution rate parameter

  ## This should likely adjust PLSD, be part of rejection, or both
  
  RSS_Diff=RSS-RSS_Min   
    
  ## Set the updated shape and rate parameters - may or may not coincide with RSS2_post
  
  shape2= shape + n_obs/2
  rate2 =rate + RSS_Min/2

  ## These can be initialized here or even at the top of the function
  
  disp_out<-matrix(0,nrow=n,ncol=1)
  beta_out<-matrix(0,nrow=n,ncol=ncol(x))
  test_out<-matrix(0,nrow=n,ncol=3)
  iters_out<-matrix(0,nrow=n,ncol=1)
  
  # Copy Envelope into temporarily envelope (that is modifed based on the simulated dispersion)
  # Added Env2$cbars_int and Env2$cbars_slope
  
  
  Env2$cbars_int=Env2$cbars-Env2$cbars_slope
  
  Env_temp=Env2

  ## This Should correspond to the intercept terms for the LL (from prior)
  
  NegLL_int=Env2$NegLL-Env2$NegLL_slope

##  print("NegLL_temp_part1")
##  print(NegLL_temp_part1)  

  ## This part comes from the likelihood function - may want to eventually
  ## only calculate the part that depends on beta here.....
  
  ##################################################################

    ### Test .set_Grid_Cpp and .setlogP_cpp functions to be removed once confirmed 
    ## Everything is working

  NegLL_temp_part2=0.5*RSS*(1/dispstar[1,1])+rep((n_obs/2)*log(2*pi),length(RSS))-rep((n_obs/2)*log(1/dispstar[1,1]),length(RSS) )
  
  #Set_Grid_cpp(GIndex, cbars, Lint)
  
  GIndex=Env2$GridIndex
  cbars_temp=Env2$cbars
  Lint=Env2$Lint1
  
  
  outgrid=.Set_Grid_cpp(GIndex, cbars_temp, Lint)
  outgrid_int=.Set_Grid_cpp(GIndex, Env2$cbars_int, Lint)
  
  logP_temp=outgrid$logP
  logP_int=outgrid_int$logP
  
  NegLL_temp=NegLL_int+NegLL_temp_part2
  G3_temp=Env2$thetabars
  
  #outlogP=.setlogP_cpp(logP, NegLL, cbars, G3)
  
  # compute this for the intercept as well to get LLconst_int
  
  outlogP=.setlogP_cpp(logP_temp, NegLL_temp, cbars_temp, G3_temp)
  outlogP_int=.setlogP_cpp(logP_int, NegLL_int, Env2$cbars_int, G3_temp)
  
  logP_temp2=outlogP$logP
  #logP2_int=outlogP_int$logP
  
  
  LLconst_temp=outlogP$LLconst
  LLconst_int=outlogP_int$LLconst

  Env2$LLconst_int=outlogP_int$LLconst
  
  
#  print("LL Const- full")
#  print(LLconst_temp)
# print("LL Const- intercept")
#  print(LLconst_int)
  
  
      
  logP2 = logP_temp2[,2]
  
  maxlogP=max(logP2)
  
  PLSD=exp(logP2-maxlogP)
  
  PLSD=PLSD/sum(PLSD)
  
#  print("PLSD_Original")
#  print(PLSD)
  
  ########################  End of test ############################
  
  
## Loop through and accept/reject based on test from internal function

  for(i in 1:n){  
    a1=0  
    
    iters_out[i]=1
    
    while(a1==0){
      p=rgamma(1,shape=shape2,rate=rate2)  
      dispersion=1/p
      
      ## Update wt2

      wt2=wt/rep(dispersion,length(y))

      
      ## Update the Grid
      
      # Step 1: Update Cbars and NegLL

      disp_ratio=(dispstar/dispersion)
      temp=disp_ratio[1,1]*Env2$cbars_slope
      Env_temp$cbars=Env2$cbars_int+temp
      
      NegLL_temp_part2=0.5*RSS*(1/dispersion)+rep((n_obs/2)*log(2*pi),length(RSS))-rep((n_obs/2)*log(1/dispersion),length(RSS) )
      NegLL_temp=NegLL_int+NegLL_temp_part2

      Env_temp$NegLL=NegLL_temp
        
      ## Step 2:  Call Set_Grid using revised cbars but original GIndex and Lint (intervals)
      ##          and update some parts of Envelope using the results
      
      GIndex=Env_temp$GridIndex
      Lint=Env_temp$Lint1
      
      outgrid=.Set_Grid_cpp(GIndex, Env_temp$cbars, Env_temp$Lint1)
      
      Env_temp$loglt=outgrid$lglt
      Env_temp$logrt=outgrid$lgrt
      Env_temp$logU=outgrid$logU
      
      ## Step 3:  Call set_logP
      
      logP_temp=outgrid$logP

      #outlogP=.setlogP_cpp(logP, NegLL, cbars, G3)
      outlogP=.setlogP_cpp(logP_temp, Env_temp$NegLL, Env_temp$cbars, Env_temp$thetabars)
      
      logP_temp2=outlogP$logP
      Env_temp$LLconst=outlogP$LLconst

      
            
      ## Compute final constants
      ## Be careful with if brining back to *.cpp --> Index is one less
      
      ## Check if correction can be made here - difference form Normal case (be careful)!
      ## Unclear if type of objects are the same (length should match)

      logP2 = logP_temp2[,2]
      #logP2 = logP_temp2[,2]-0.5*p*RSS_Diff
      
      maxlogP=max(logP2)
      
      PLSD=exp(logP2-maxlogP)

      ## Be careful with if brining back to *.cpp --> Index is one less

      Env_temp$logP=logP_temp2[,1]
      Env_temp$PLSD=PLSD/sum(PLSD)
  
##      print("precision for candidate")
##      print(p)
##      print("PLSD Inside loop")
##      print(Env_temp$PLSD)

##      logP2_alt = logP_temp2[,2]-0.5*p*RSS_Diff
##      maxlogP_alt=max(logP2_alt)
##      PLSD_alt=exp(logP2_alt-maxlogP_alt)
##      PLSD_alt=PLSD_alt/sum(PLSD_alt)
      
##      print("PLSD Alternative Inside loop")
##      print(PLSD_alt)
      
      ## Last constant is likely not changed (or used in sampling)  
      
      ## Verify envelope with updated cbars and dispersion can be sampled 
      ### (should not be the final sampling as envelope must be updated each iteration)
      
      ## This now likely uses the correct envelope but should only generate one candidate
      ## and not continue on until acceptance
      
##      print("thetabars used by envelope")
##      print(Env_temp$thetabars)
      
      
      ### Do away with adjusted evelope use original
      
      #sim=.rindep_norm_gamma_reg_std_cpp(n=1,y=y,x=x2,mu=mu2,P=P2,alpha=alpha,wt=wt2,
      #                                   f2=f2,Envelope=Env_temp,family="gaussian",link="identity",
      #                                  as.integer(0))
      sim=.rindep_norm_gamma_reg_std_cpp(n=1,y=y,x=x2,mu=mu2,P=P2,alpha=alpha,wt=wt2,
                         f2=f2,Envelope=Env2,family="gaussian",link="identity",as.integer(0))
      
      
      ## Add 1 here becasuse arrays in *.cpp start with 0 instead of 1
      J_out=sim$J+1
      log_U2=sim$log_U2
      
      #print("Selected component")
      #print(J_out)
      
      disp_out[i,1]=dispersion
      beta_out[i,1:ncol(x)]=sim$out[1,1:ncol(x)]
      #iters_out[i]=sim$draws[1]
      
      ## Simulation from updated grid goes here (including rejection calculation)
      ## Update test to include penalty for high dispersion
      ## we need to pass element of grid for this to work
      
      # This might penalize too much - mean gets shifted to the prior mean
      
      test1=-0.5*p*RSS_Diff[J_out]

      #print("Implied RSS penalty")
      #print(test1)
      
      ## Just scale the test by the dispersion ratio and add test1
      ## This would likely be correct if there was no prior component part of the test
      
      #test=disp_ratio*sim$test +test1-log_U2

#      print("Test components")
#      print(sim$test)
#      print(sim$test_int)
#      print(sim$test_data)
      
      #test=disp_ratio*sim$test -log_U2
      test=sim$test_int+disp_ratio*sim$test_data -log_U2
   
      ## We should subtract squared term here that is added to the normal
      
      
         
      # See what adding test 1 would do
      #test=sim$test_int+disp_ratio*sim$test_data +test1 -log_U2
      
      #test=sim$test+test1
      
      if(test>=0) a1=1
      else{iters_out[i]=iters_out[i]+1}        
      
    }        
  }

 ## Undo standardization

  ## implied mean corresponding to mu2
  
#  thetabars_real=L2Inv%*%L3Inv%*%t(Env_temp$thetabars)
  
#  for(i in 1:9){
#    thetabars_real[,i]=thetabars_real[,i]+mu
#  }
  
#  print("thetabars_real")
#  print(thetabars_real)
     
  out=L2Inv%*%L3Inv%*%t(beta_out)
  
  for(i in 1:n){
    out[,i]=out[,i]+mu
  }
  
#    return(list(coefficients=out,dispersion=disp_out,iters=iters_out))
  
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



