rglmb_rand<-function(n=1,y,x,mu,P_0,P,wt=1,dispersion=NULL,
                     nu=NULL,V=NULL,family=gaussian(),offset2=NULL,start=NULL,Gridtype=2)
{
  
  if(is.numeric(n)==FALSE||is.numeric(y)==FALSE||is.numeric(x)==FALSE||
       is.numeric(mu)==FALSE||is.numeric(P)==FALSE) stop("non-numeric argument to numeric function")
  
  x <- as.matrix(x)
  mu<-as.matrix(as.vector(mu))
  P<-as.matrix(P)    
  P_0<-as.matrix(P_0)    
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) 
    rownames(y)
  else names(y)
  if(length(n)>1) n<-length(n)     
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(offset2)) 
    offset2 <- rep(0, nobs)
  nvars2<-length(mu)  
  if(!nvars==nvars2) stop("incompatible dimensions")
  if (!all(dim(P_0) == c(nvars2, nvars2))) 
    stop("incompatible dimensions")
  if(!isSymmetric(P_0))stop("matrix P_0 must be symmetric")
  tol<- 1e-06 # Link this to Magnitude of P  
  eS <- eigen(P_0, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'P_0' is not positive definite")
  
  if (is.null(start)) 
    start <- mu
  if (is.null(offset2)) 
    offset2 <- rep.int(0, nobs)
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
      
      famfunc<-glmbfamfunc(family)
      f1<-famfunc$f1
      f2<-famfunc$f2
      f3<-famfunc$f3
      f5<-famfunc$f5
      f6<-famfunc$f6
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
  if(family$family=="poisson"||family$family=="binomial")dispersion2<-1
  else dispersion2<-dispersion
  
  mutemp<-x%*%start
  xtemp<-diag(1) 
  
  betaout<-matrix(0,nrow=n,ncol=length(y))
  alphaout<-matrix(0,nrow=n,ncol=length(mu))
  betatemp<-mutemp
  
  
  check<-rglmb_rand_cpp(n=n,y=y,x=x,
                        mu=mu,P_0=P_0,P=P,offset2=offset2
                        ,wt=wt,
                        dispersion=dispersion,
                        famfunc=famfunc,f1=f1,f2=f2,f3=f3,
                        start=mu,Gridtype=Gridtype) 
  return(check)
  
}
