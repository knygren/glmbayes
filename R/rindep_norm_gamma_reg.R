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
#' @example inst/examples/Ex_rindep_norm_gamma_reg.R
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
  n_obs=length(y)

  shape2= shape + n_obs/2
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

  prior_list2=list(mu=mu,Sigma=Sigma, dispersion=dispersion2)
  glmb_out1=glmb(n=1,y~x-1,family=gaussian(),prior=prior_list2)
  
  ## Use this as betastar
  
  betastar=glmb_out1$coef.mode
  xbetastar=x%*%betastar
  RSS2_post=t(y-xbetastar)%*%(y-xbetastar)  ## Residual SUM of SQUARES at post model

  ## Updated version (when using the likelihood subgradient approach
  ## The RSS at the tangent point gets shifted to the gamma distribution
  
  rate2 =rate + RSS2_post/2

#  print("dispersion from optimization")
#  print(dispersion2)
  
  disp_temp=RSS2_post/n_obs
  
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



