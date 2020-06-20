#' Simulate New Data
#'
#' Simulates New Data tied to the model predictions 
#' @param pred Predictions from the predict function or the original glmb object
#' @param family family used in model
#' @param wt weighting variable (that incorporates the dispersion)
#' @return Simulated values for data corresponding to the model predictions
#' @example inst/examples/Ex_confint.glmb.R


glmb_simulate<-function(pred,family,wt=1){
  
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
  
  for(i in 1:nsims){
    
    #rgamma(n, shape, rate = 1, scale = 1/rate)
    #rnorm(n, mean = 0, sd = 1)
    #rbinom(n, size, prob)
    ## Unclear about quasibinomial
    
    if(family=="poisson") y_temp[i,1:nvars]=rpois(n=nvars,pred[i,1:nvars])              
    if(family=="quasipoisson") y_temp[i,1:nvars]=rpois(n=nvars,pred[i,1:nvars])             
    ### Verify this part - rather complicated
    if(family=="Gamma") y_temp[i,1:nvars]=rgamma(n=nvars,shape=wt,(1/wt)*pred[i,1:nvars])              
    if(family=="binomial") y_temp[i,1:nvars]=rbinom(n=nvars,size=round(wt),prob=pred[i,1:nvars])              
    if(family=="quasibinomial") y_temp[i,1:nvars]=rbinom(n=nvars,size=round(wt),prob=pred[i,1:nvars])                            
    
    ## Verify this part - rather complicated    
    if(family=="gaussian") y_temp[i,1:nvars]=rnorm(n=nvars,mean=pred[i,1:nvars],sd=sqrt(1/wt))              
    
    
  }
  
}
