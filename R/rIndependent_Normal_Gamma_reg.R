#' The Bayesian Indendepent Normal-Gamma Regression Distribution
#'
#' Generates iid samples from the posterior density for the 
#' independent normal-gamma regression
#' @param prior_list a list with the prior parameters (mu, Sigma, shape and rate) for the 
#' prior distribution.
#' @param offset an offset parameter
#' @param weights a weighting variable
#' @param max_disp_perc parameter currently used to control upper bound in accept-reject procedure 
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
                                         max_disp_perc=0.99){
  
  call<-match.call()
  
  offset2=offset
  wt=weights
  
  if(length(wt)==1) wt=rep(wt,length(y))
  
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
  
  lm_out=lm(y ~ x-1,weights=wt,offset=offset2) # run classical regression to get maximum likelhood estimate
  RSS=sum(residuals(lm_out)^2)
  
  RSS_ML=sum(residuals(lm_out)^2)
  n_obs=length(y)
  
  
  dispersion2=dispersion
    RSS_temp<-matrix(0,nrow=1000)
  
  for(j in 1:10){
    glmb_out1=glmb(y~x-1,family=gaussian(),
    dNormal(mu=mu,Sigma=Sigma,dispersion=dispersion2),weights=wt,offset=offset2)
    
    res_temp=residuals(glmb_out1)
    
    for(k in 1:1000){
      RSS_temp[k]=sum(res_temp[k,1:length(y)]*res_temp[k,1:length(y)])
    }
    
    RSS_Post2=mean(RSS_temp)
    b_old=glmb_out1$coef.mode
    xbetastar=x%*%b_old
    RSS2_post=t(y-xbetastar)%*%(y-xbetastar)  
    shape2= shape + n_obs/2
  #  rate2 =rate + RSS2_post/2  # 38 candidates per acceptance
    rate2=rate + RSS_Post2/2    # 35.7 candidates per acceptance [though some with positive acceptance]
  
    dispersion2=rate2/(shape2-1)

  }
  
  
  betastar=glmb_out1$coef.mode
  dispstar=rate2/(shape2-1)

  
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

  ## Update Gamma parameters
  
  shape2= shape + n_obs/2
  rate2 =rate + (RSS_ML/2)
  #rate3 =rate + (RSS2_post/2)
  rate3 =rate + (RSS_Post2/2)
  

  #################  this Block Does evaluations at lower and upper bounds   #################

  cbars=Env2$cbars
  
  gs=nrow(Env2$cbars)
  logP1=Env2$logP
  New_LL=c(1:gs)
  New_logP2=c(1:gs)
  
  upp=1/qgamma(c(1-max_disp_perc),shape2,rate3)
  low=1/qgamma(c(max_disp_perc),shape2,rate3)
  wt_upp=wt/rep(upp,length(y))
  wt_low=wt/rep(low,length(y))
  
  theta_upp=Inv_f3_gaussian(t(Env2$cbars), y, as.matrix(x2),as.matrix(mu2,ncol=1), as.matrix(P2), 
                            as.vector(alpha), as.vector(wt_upp))

  theta_low=Inv_f3_gaussian(t(Env2$cbars), y, as.matrix(x2),as.matrix(mu2,ncol=1), as.matrix(P2), 
                            as.vector(alpha), as.vector(wt_low))

  
  New_LL_upp=c(1:gs)
  New_LL_low=c(1:gs)
  
  for(j in 1:gs){
    
    theta_temp_upp=as.matrix(theta_upp[j,1:ncol(x)],ncol=1)
    theta_temp_low=as.matrix(theta_low[j,1:ncol(x)],ncol=1)
    cbars_temp=as.matrix(cbars[j,1:ncol(x)],ncol=1)
    
    # New formula - This should likely replace: -NegLL(i)+arma::as_scalar(G3row.t() * cbarrow)
    
    New_logP2[j]=logP1[j]+0.5*t(cbars_temp)%*%cbars_temp
    New_LL_upp[j]=-0.5*t(theta_temp_upp)%*%P2%*%theta_temp_upp+t(cbars_temp)%*% theta_temp_upp
    New_LL_low[j]=-0.5*t(theta_temp_low)%*%P2%*%theta_temp_low+t(cbars_temp)%*% theta_temp_low
    
        
  }

  maxlogP=max(New_logP2)
  PLSD_new=exp(New_logP2-maxlogP)
  sumP=sum(PLSD_new)
  PLSD_new=PLSD_new/sumP
  log_P_diff=log(PLSD_new)-log(Env2$PLSD)
  
    
  LL_temp_upp=log_P_diff+New_LL_upp
  LL_temp_low=log_P_diff+New_LL_low
  
  max_New_LL_upp=max(LL_temp_upp)
  max_New_LL_low=max(LL_temp_low)

  slope_out=(max_New_LL_upp-max_New_LL_low)/(upp-low)
  int_out=max_New_LL_low-slope_out*low
  lmc=c(int_out,slope_out)
  
  a1=2
  b1=(upp-low)
  c1=-log(upp/low)
  dispstar= (-b1+ sqrt(b1^2-4*a1*c1))/(2*a1) # Point of tangency
   
  
  lm_log2=lmc[2]*dispstar
  #lm_log1=lmc[1]+lm_log2-lm_log2*log(dispstar)
  lm_log1=lmc[1]+lmc[2]*dispstar-lm_log2*log(dispstar)
  
  shape3=shape2-lm_log2
  
  # Compute log_P_diff
  
  max_New_LL_UB = max_New_LL_upp
  max_LL_log_disp=lm_log1+lm_log2*log(upp) ## From above
  
  ## Outputs from this block
  
  
  ## 1) upp, low
  ## 2) log_P_diff
  ## 3) lm_log1,lm_log2
  ## 4) shape3
  ## 5) max_New_LL_UB
  ## 6) max_LL_log_disp
    
  ## Not currently from block: RSS_ML, rate2 
  
  #####################################################################################

  
  gamma_list=list(shape3=shape3,rate2=rate2,disp_upper=upp,disp_lower=low)
  
  UB_list=list(RSS_ML=RSS_ML,max_New_LL_UB=max_New_LL_UB,
               max_LL_log_disp=max_LL_log_disp,lm_log1=lm_log1,lm_log2=lm_log2, 
               log_P_diff=log_P_diff)
  
##  ptm <- proc.time()
  sim_temp=.rindep_norm_gamma_reg_std_V2_cpp (n=n, y=y, x=x2, mu=mu2, P=P2, alpha=alpha, wt,
  f2=f2, Envelope=Env2, 
  gamma_list=gamma_list,
  UB_list=UB_list,
  family="gaussian",link="identity", progbar = as.integer(0))
  
##  print("time for *.cpp function")
##  print(proc.time()-ptm)

##  print("mean candidates per acceptance - *.cpp function")
##  print(mean(sim_temp$iters_out))
  
    
  beta_out=sim_temp$beta_out
  disp_out=sim_temp$disp_out
  iters_out=sim_temp$iters_out
  weight_out=sim_temp$weight_out

  out=L2Inv%*%L3Inv%*%t(beta_out)
  
  for(i in 1:n){
    out[,i]=out[,i]+mu
  }
  
  
    ### Call standard simulation function
##  ptm <- proc.time()
  
##  sim_temp=rindep_norm_gamma_reg_std_R(n=n,y=y,x=x2,mu=mu2,P=P2,alpha=alpha,wt=wt,
##  f2=f2,Envelope=Env2,
##  gamma_list=gamma_list,
##  UB_list=UB_list,
##  family="gaussian",link="identity",as.integer(0))

##  print("time for *.R function")
##  print(proc.time()-ptm)
  
##  print("mean candidates per acceptance *.R function")
##  print(mean(sim_temp$iters_out))
  
  ##  *.cpp function took 6.54 seconds and R function 83.71
  ## Both had ~38 candidates per accepted simulated value
  
#########################################  Post Processing for simulation function handling loop
  
##  beta_out=sim_temp$beta_out
##  disp_out=sim_temp$disp_out
##  iters_out=sim_temp$iters_out
##  weight_out=sim_temp$weight_out

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
    weight_out=weight_out,
    sim_bounds=list(low=low,upp=upp)
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


