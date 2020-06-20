#' Simulate Responses
#'
#' Simulate responses from the distribution corresponding to a fitted \code{glmb} object.
#' @param object An object of class \code{glmb}, typically the result of a call to the 
#' function \code{glmb}.
#' @param nsim Defunct (see below).
#' @param seed an object specifying if and how the random number generator should be 
#' initialized (seeded).
#' @param \ldots Additional arguments passed to the function. Will 
#' frequently include a matrix pred of simulated predictions from
#' the predict function, the family (e.g., binomial) and an optional
#' vector of weights specifying prior.weights for the simulated values (default is 1)
#' @return Simulated values for data corresponding to simulated model predictions that correspond either 
#' to the original data or to a \code{newdata} data frame provided to the predict function. 
#' @example inst/examples/Ex_confint.glmb.R

simulate.glmb<-function(object,nsim=1,seed=NULL,...){
  
  method_args=list(...)
  pred=method_args$pred
  family=method_args$family
  wt=method_args$wt
#  print(method_args)
  #stop()
    ## pred
  ## family
  ## wt
  
  ## For Poisson and binomial, dispersion2=1
  ## wt2=wt/dispersion
  
  ## Note:
  ## gaussian
  ## Xb<-alpha+x%*%b
  ## --> sd=sqrt(1/(wt))
  
  
  ## binomial
  ## size=round(wt)
  ## y = sim/size
  
  
  ##  gamma
  ##  dispersion2=1/wt 
  ##  shape=1/dispersion2
  ##  scale=mu*dispersion2
  
  ## Set up temporary matrix to hold these elements
  nvars=ncol(pred)
  nsims=nrow(pred)
  y_temp<-matrix(0,nrow=nrow(pred),ncol=ncol(pred))
  
#  print(nvars)
#  print(nsims)
#  print(y_temp)
  
  for(i in 1:nsims){
    
    #rgamma(n, shape, rate = 1, scale = 1/rate)
    #rnorm(n, mean = 0, sd = 1)
    #rbinom(n, size, prob)
    ## Unclear about quasibinomial
    
    if(family=="poisson") y_temp[i,1:nvars]=rpois(n=nvars,pred[i,1:nvars])              
    if(family=="quasipoisson") y_temp[i,1:nvars]=rpois(n=nvars,pred[i,1:nvars])             
    ### Verify this part - rather complicated
    if(family=="Gamma") y_temp[i,1:nvars]=rgamma(n=nvars,shape=wt,(1/wt)*pred[i,1:nvars])              

    ## Simulate from binomial and then divide by the weight!    
    if(family=="binomial") y_temp[i,1:nvars]=rbinom(n=nvars,size=round(wt),prob=pred[i,1:nvars])/round(wt)              
    if(family=="quasibinomial") y_temp[i,1:nvars]=rbinom(n=nvars,size=round(wt),prob=pred[i,1:nvars])/round(wt)                            
    
    ## Verify this part - rather complicated    
    if(family=="gaussian") y_temp[i,1:nvars]=rnorm(n=nvars,mean=pred[i,1:nvars],sd=sqrt(1/wt))              
    
    
  }
  
  return(y_temp)
}
