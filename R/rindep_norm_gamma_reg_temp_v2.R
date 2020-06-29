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
  
  ## Set the updated shape and rate parameters - may or may not coincide with RSS2_post
  
  shape2= shape + n_obs/2
  rate2 =rate + RSS_Min/2

  ## These can be initialized here or even at the top of the function
  
  disp_out<-matrix(0,nrow=n,ncol=1)
  beta_out<-matrix(0,nrow=n,ncol=ncol(x))
  test_out<-matrix(0,nrow=n,ncol=3)
  iters_out<-matrix(0,nrow=n,ncol=1)
  
  # Copy Envelope into temporarily envelope (that is modifed based on the simulated dispersion)
  
  Env2$cbars_int=Env2$cbars-Env2$cbars_slope
  
  Env_temp=Env2

  ## This Should correspond to the intercept terms for the LL (from prior)
  
  NegLL_temp_part1=Env2$NegLL-Env2$NegLL_slope

  print("NegLL_temp_part1")
  print(NegLL_temp_part1)  

  ## This part comes from the likelihood function - may want to eventually
  ## only calculate the part that depends on beta here.....
  
  NegLL_temp_part2=0.5*RSS*(1/dispstar[1,1])+rep((n_obs/2)*log(2*pi),length(RSS))-rep((n_obs/2)*log(1/dispstar[1,1]),length(RSS) )
 
  print("NegLL_temp_part2")
  print(NegLL_temp_part2)  
  
  print("NegLL_temp_total")
  print(NegLL_temp_part1+NegLL_temp_part2)  
  
  print("NegLL from Optimized Envelope call")
  print(Env2$NegLL)
  
## For now, just loop through n and accept all candidates
  
  print("Looping and printing calculated NegLL_temp")
  
  for(i in 1:n){  
    a1=0  
    while(a1==0){
      p=rgamma(1,shape=shape2,rate=rate2)  
      dispersion=1/p
      
      ## Calculations of updated Grid to here
      
      # Update cbars and wt2 also need to update LL and all related constants

      disp_ratio=(dispstar/dispersion)
      temp=disp_ratio[1,1]*Env2$cbars_slope
      
      Env_temp$cbars=Env2$cbars_int+temp
      wt2=wt/rep(dispersion,length(y))
      
      NegLL_temp_part2=0.5*RSS*(1/dispersion)+rep((n_obs/2)*log(2*pi),length(RSS))-rep((n_obs/2)*log(1/dispersion),length(RSS) )
      
      NegLL_temp=NegLL_temp_part1+NegLL_temp_part2
      print(NegLL_temp)
      
      ## Verify envelope with updated cbars and dispersion can be sampled 
      ### (should not be the final sampling as envelope must be updated each iteration)
      
      ## As this uses "incorrect" envelope, simulation could be slow....
      
      sim=rnnorm_reg_std(n=1,y=y,x=x2,mu=mu2,P=P2,alpha=alpha,wt=wt2,
                         f2=f2,Envelope=Env_temp,family="gaussian",link="identity",as.integer(0))
      
      
      disp_out[i,1]=dispersion
      beta_out[i,1:ncol(x)]=sim$out[1,1:ncol(x)]
      iters_out[i]=sim$draws[1]
      ## Simulation from updated grid goes here (including rejection calculation)
      
  
      a1=1
            }
  }
  
  return(list(coefficients=beta_out,dispersion=disp_out,iters=iters_out))
  
  return(list(n=as.integer(n),y=as.vector(y),x2=as.matrix(x2),
              mu2=as.matrix(mu2,ncol=1),P2=as.matrix(P2),
              alpha=as.vector(alpha),wt2=as.vector(wt2),f2=f2,
              Env=Env2))
  
  
#  return(Env2)
  
  
    
#  print("dispersion that maximizes log-likelihood")
#  print(disp_temp)
  
  ## Set updated gamma parameters for the candidate generation
  ## shape matches posterior while rate will be too low (adjusted using accept-rejecte)
  
  
  # set up matrices to hold output
  
  disp_out<-matrix(0,nrow=n,ncol=1)
  beta_out<-matrix(0,nrow=n,ncol=ncol(x))
  test_out<-matrix(0,nrow=n,ncol=3)
  iters_out<-matrix(0,nrow=n,ncol=1)
  
  # Internal function used below to determine acceptance rate
  
  testfun<-function(beta,betastar,p,y,x,RSS){
    
    beta=as.matrix(beta,ncol=1) ## Simulated candidate
    betastar=as.matrix(betastar,ncol=1)   ## Conditional posterior mode
    xbeta=x%*%beta
    xbetastar=x%*%betastar
    RSS2_test=t(y-xbeta)%*%(y-xbeta)  ## Residual SUM of SQUARES for test candidate
    RSS2_post=t(y-xbetastar)%*%(y-xbetastar)  ## Residual SUM of SQUARES at post model
    
    # Log-Acceptance rate
    
    # 1)if we used the prior to generate candidate for beta, 
    # we would have test1=-0.5*p*(RSS2_test-RSS_ML) and 
    ## RSS_ML would shift to gamma distribution
    # 
    # 2) when we replace prior with likelihood sub-gradient density, we replace  
    # RSS2_test with (RSS_post+2*t(y-xbetastar)%*%x%*%(beta-betastar))
    # and replace RSS_ML with RSS_Post (which gets shifted to the gamma distribution)  
    
    ## This would be test if sampled beta from prior
    
#    test1=-p*0.5*(RSS2_test-RSS_ML)
    test1=0
    
    ## The below should be a concave function that obtains its max when beta=betastar
    ## When beta=betastar, RSS2_test=RSS_Post and the slope term=0
    ## So whole terms should be 0

    test2=p*0.5*(RSS2_post-2*t(y-xbetastar)%*%x%*%(beta-betastar)-RSS2_test)

    ## This should likely be -p*0.5*(RSS2_post - RSS_LB)
    ## when using grid where RSS_LB is lower bound among tangent points
      
#    test1=-p*0.5*(RSS2_post-RSS) 
    


     test=test1+test2
     

    return(list(test=test[1,1],test1=test1,test2=test2))
  }
  
  
  ## Look through iterations to generate candidates until acceptance
  
  for(i in 1:n){  
    
    # Initialize flag for acceptace to 0 and count of iters to 1
    
    a1=0  
    iters_out[i,1]=1  
    
    while(a1==0){
      
      p=rgamma(1,shape=shape2,rate=rate2)  
      dispersion=1/p
      
      ## Set prior list in order to get the conditional posterior mode from glmb function
      ## inefficient but works for now.
      
      prior_list3=list(mu=mu,Sigma=Sigma, dispersion=dispersion)
      glmb_out1=glmb(n=1,y~x-1,family=gaussian(),prior=prior_list3)
      
      betatest=as.matrix(mvrnorm(n = 1, mu=betastar, Sigma=Sigma, tol = 1e-6, empirical = FALSE),ncol=1)

      testtemp=testfun(betatest,betastar,p,as.matrix(y,ncol=1),x,RSS)
      test=exp(testtemp$test)
      
      disp_out[i,1]=dispersion
      
      beta_out[i,1:ncol(x)]=betatest 
      test_out[i,1]=exp(testtemp$test)
      test_out[i,2]=exp(testtemp$test1)
      test_out[i,3]=exp(testtemp$test2)
      
      # set a1 to 1 if draw was accepted to end while loop
      
      if(runif(1)<test) a1=1
      ## increment iters count if candidate not accepted 
      
      iters_out[i,1]=iters_out[i,1]+1  
      
    }
    
  }
  
  # Return list with elements similar to the rnorm_gamma function
  # For now, leave PostMode as NULL - Needs iterative procedure or
  # customized optimization to determine accurately (betastar depends on dispersion)
  # Also leave envelope, and loglike null for now 
  # add test_out to list for now to assess accept-reject implementation
  
  famfunc=glmbfamfunc(gaussian())  
  f1=famfunc$f1
  
    outlist=list(
coefficients=beta_out, 
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



