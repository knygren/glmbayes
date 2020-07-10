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

rindependent_norm_gamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian()){
  
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
    
    dispersion2=RSS2_post/n_obs
    
  }
  
  prior_list2=list(mu=mu,Sigma=Sigma, dispersion=dispersion2,RSS2_post=RSS2_post)
  glmb_out1=glmb(y~x-1,family=gaussian(),dNormal(mu=mu,Sigma=Sigma,dispersion=dispersion2))
  
  #  glmb_out1=glmb(n=1,y~x-1,family=gaussian(),prior=prior_list2)
  
  ## Use this as betastar
  
  betastar=glmb_out1$coef.mode
  
  
  xbetastar=x%*%betastar
  RSS2_post=t(y-xbetastar)%*%(y-xbetastar)  ## Residual SUM of SQUARES at post model
  
  ## Updated version (when using the likelihood subgradient approach
  ## The RSS at the tangent point gets shifted to the gamma distribution
  
  #shape2= shape + n_obs/2
  #rate2 =rate + RSS2_post/2
  
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
  
  ## Pull the initial Envelope based on optimized values above
  
  Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
                     as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),
                     as.vector(wt2),
                     family="gaussian",link="identity",
                     Gridtype=Gridtype, n=as.integer(n),
                     sortgrid=TRUE)
  
    
  #### End of Standardization - Steps below might change  


  ## Update Gamma parameters
  
  shape2= shape + n_obs/2
  rate2 =rate + RSS_ML/2
  
  ## Copy Envelope for use in simulation (as needed update components)
  
  Env_temp=Env2

  ### Add step here to find p_prec which optimizes cbars*thetabars (cbars(p_prec)) for each dimension
  ### Then select the thetabars_opt which are implied to derive bound for cbars*thetabars
    
  ########################  End of test ############################
  
  ## Because of Acceptance procedure below, this should use RSS_ML not RSS_Post from above
  
  iters_out<-c(1:n)
  
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
      

      LL_New=-f2_gaussian_vector(t(Env_temp$thetabars), y, as.matrix(x2), as.matrix(mu2,ncol=1),
                                   as.matrix(P2), as.vector(alpha), as.vector(wt2))
      
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
      
      #f3=famfunc$f3        
      
      ## These should match original - can run check
      
      cbars_new=f3(as.matrix(New_thetabars[J_out,],ncol=1), y, as.matrix(x2), as.matrix(mu2,ncol=1),
                   as.matrix(P2), as.vector(alpha), as.vector(wt2))
      
      ntheta=as.matrix(New_thetabars[J_out,1:ncol(x)],ncol=1)
      ntheta_star=as.matrix(-P2%*%as.matrix(cbars_new,ncol=1),ncol=1)  
      
      betadiff=as.matrix(sim$out[1,1:ncol(x)],ncol=1)-as.matrix(New_thetabars[J_out,1:ncol(x)],ncol=1)
      
      test=0
      
      ## This is an Upper bound for LL in terms of LL at tangent point and the gradient 
      
      UB1=LL_New[J_out]-t(cbars_new)%*%(betadiff)

      ##    UB1 contains a term that is quadratic in New_thetabars[J_out,1:ncol(x)]
      ##    UB2 subtracts that same term and bounds the term using ntheta_star 
      ##    Bounding term involving ntheta_star should be pre-calculated
      ##    for now, do it here
      
      UB2= (0.5 * t(ntheta)%P2*%ntheta +t(cbars_new)%*% ntheta 
            -(0.5*t(ntheta_star)%P2*%ntheta_star+t(cbars_new)%*% ntheta_star ) )
      
      ## UB1 also contains a term that involves -0.5*(1/dispersion) * RSS(ntheta)
      ## This term is then re-entered here and bounded
     
      yxbeta=(y-alpha-x2%*%ntheta)*sqrt(wt)
       
      RSS_ntheta=t(yxbeta)%*%(yxbeta)
      
      # -0.5*(1/dispersion)*RSS_ntheta <= -0.5*(1/dispersion) RSS_ML
      
      UB3  =  0.5*(1/dispersion)*(RSS_ntheta-RSS_ML)
            

      test1= LL_Test-UB1
      test2= LL_Test-(UB1+UB2)
      test3= LL_Test-(UB1+UB2+UB3)

      print("tests 1, 2, and 3 - Each test should get a bit more negative")
      print(test1)
      print(test2)
      print(test3)
      
      stop("test values printed above")
      
      ## UB1 --> Should be test if we are using prior for beta but shift the two components to gamm
      ## UB2 --> Should be test if we were sampling for the conditional density using the grid positioned at betastar, dispstar
      ## UB3 --> Should be test if we are sampling from the density with different variances but keeping mean of
      ##         normal densities fixed 
      ##         unclear if/how the probabilities of visiting different parts of the grid should be adjusted based on what the dipserson
      ##         is. Next step is to figure this part out
      
      
      #UB1=-0.5*(1/dispersion)*RSS_ML-(n_obs/2)*log(dispersion)-(n_obs/2)*log(2*pi)  
      #UB2=-NegLL_New[J_out]-t(cbars_new)%*%(betadiff)
      #UB3=UB1-t(cbars_new)%*%(betadiff)
      
      #Diff1=LL_Test-UB1
      #Diff2=LL_Test-UB2
      #Diff3=LL_Test-UB3
      
      #test=Diff3-log_U2
      
      
      
     a1=1
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
    loglike=NULL
    #,test_out=test_out
    )
  
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



