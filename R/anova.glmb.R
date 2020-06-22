#' Analysis of Deviance for Bayesian Generalized Linear Model Fits
#'
#' Compute an analysis of deviance table for one (current implementation) or more (future)
#' Bayesian generalized linear model fits,
#' @param object an object of class \code{glmb}, typically the result of a call to \link{glmb}
#' @param \ldots Other arguments passed to or from other methods. 
#' @return An object of class \code{"anova"} inheriting from class \code{"data.frame"}.
#' @details Specifying a single object (currently only implementation) gives a sequential 
#' analysis of deviance table for that fit. That is the reductions in residual deviance as each
#' term of the formula is addded in turn are given as the rows of a table, plus the residual deviances
#' themselves.  
#' @example inst/examples/Ex_residuals.glmb.R
#' @export
#' @method anova glmb



anova.glmb<-function(object,...){
  
  # Gather information for full model (ideally would grab gridtype as well but not available)
  
  n=nrow(object$coefficients)
  mu=object$Prior$mean
  V=object$Prior$Variance
  obj_family=family(object)
  
  
  n_obs=nobs(object)
  ff_all=formula(object)
  mf_all=model.frame(ff_all)  
  nvar_all=ncol(mf_all)
  terms_all=terms(ff_all)
  tl_all=attr(terms_all,"term.labels") ## terms
  nterms_all=length(tl_all)
  
  
  # Set up data frame that will hold all information and popualte with
  # available information from full model
  #"Df"         "Deviance"   "Resid. Df"  "Resid. Dev"
  
  anova_out=data.frame(pD=rep(0,(nterms_all+1)),
                       Deviance=rep(0,(nterms_all+1)),
                       'Resid. Df'=rep(nobs(object),(nterms_all+1)),
                       'Resid. Dev'=rep(0,nterms_all+1),
                       'Mod. pD'=rep(0,(nterms_all+1)),
                       DIC=rep(0,(nterms_all+1)))
  rownames(anova_out)[1]="NULL"
  rownames(anova_out)[2:(nterms_all+1)]=tl_all
  
  # Initialize nterms_left and tt2
  nterms_left=nterms_all
  tt2=terms_all
  
  while(nterms_left>0){
    
    ## Update the formula with one less term
    
    if(nterms_left>1){
      tt2=drop.terms(tt2,nterms_left,keep.response=TRUE)
      newff=formula(tt2)
    }
    else{ newff[3]=1}
    
    ## Update the data and prior with fewer variables 
    
    mf=model.frame(newff)
    mm=model.matrix(mf)
    nvar=ncol(mm)
    mu2=matrix(mu[1:nvar,1],ncol=1)
    V2=matrix(V[1:nvar,1:nvar],nrow=nvar,ncol=nvar)
    
    # Run glmb model for smaller model
    
    object2<-glmb(n=n,newff, family = obj_family,mu=mu2,Sigma=V2,Gridtype=2)
    
    # Update anova_out table
    
    anova_out[(nterms_left),3]=nobs(object)-object2$pD
    anova_out[(nterms_left),4]=colMeans(object2$deviance)
    anova_out[(nterms_left),5]=object2$pD
    anova_out[(nterms_left),6]=object2$DIC
    
    # decrement nterms_left by 1
    
    nterms_left=nterms_left-1
  }    
  
  # Wrapup by populatibg anova_out matrix
  
  anova_out[1,1]=anova_out[1,5]
  
  for(i in 1:nterms_all){
    anova_out[(i+1),1]=anova_out[(i+1),5]-anova_out[i,5]
    anova_out[(i+1),2]=anova_out[i,4]-anova_out[(i+1),4]
    
  }
  
  # Consider whether to add class info here
  
  structure(anova_out, heading = c("Analysis of Variance Table\n",
                                   paste("Response:", deparse(formula(object)[[2L]]))),
            class = c("anova", "data.frame"))# was "tabular"
  return(anova_out)
  
}



