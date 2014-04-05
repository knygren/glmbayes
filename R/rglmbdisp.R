rglmbdisp<-function(n,y,x,b,alpha=0,wt=1,shape,rate,family=gaussian()){
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  okfamilies <- c("gaussian","Gamma")
  if(family$family %in% okfamilies){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="Gamma") oklinks<-c("log")		
    if(family$link %in% oklinks)  {}
    else{stop(gettextf("link \"%s\" not available for selected family; available links are %s", 
                                                         family$link , paste(sQuote(oklinks), collapse = ", ")), 
                                                domain = NA)
    }
  }
  
  else{
    stop(gettextf("family \"%s\" not available in glmbdisp; available families are %s", 
                  family$family , paste(sQuote(okfamilies), collapse = ", ")), 
         domain = NA)
    
  }
  
  n1<-length(y)
  
  if(family$family=="gaussian"){
  y1<-as.matrix(y)-alpha
  xb<-x%*%b
  res<-y1-xb
  SS<-res*res
  
  a1<-shape+n1/2
  b1<-rate+sum(SS)/2
  
  out<-1/rgamma(n,shape=a1,rate=b1) 
  }

if(family$family=="Gamma")
  {
  
  mu1<-t(exp(alpha+x%*%b))
  
  testfunc<-function(v,wt){  
    -sum(lgamma(wt*v)+0.5*log(wt*v)+wt*v-wt*v*log(wt*v))
  }

  
  shape2=shape + 0.5 *n1
  rate1=rate +sum(wt*((y/mu1)-log(y/mu1)-1))
  
  vstar1<-shape2/rate1
  
  vout<-function(v){
    vstar1-(v/rate1)*sum((wt*digamma(wt*v) -wt*log(wt*v) + 0.5/v) )  
  }
  
  # Initialize vstar2
  vstar<-vstar1
  for(j in 1:20){
    vstar<-vout(vstar)
  }

  testbar<-testfunc(vstar,wt)
  cbar<--sum((wt*digamma(wt*vstar) -wt*log(wt*vstar) + 0.5/vstar))
  
  
  
  rate2=  rate +sum(wt*((y/mu1)-log(y/mu1)-1))-sum((wt*digamma(wt*vstar) -wt*log(wt*vstar) + 0.5/vstar) )
  
  out<-matrix(0,n)
  test<-matrix(0,n)
  a<-matrix(0,n)
  
  for(i in 1:n)
  {
    while(a[i]==0){
      out[i]<-rgamma(1,shape=shape2,rate=rate2)
      
      test[i]<-testfunc(out[i],wt)-(testbar+cbar*(out[i]-vstar))-log(runif(1,0,1))
      if(test[i]>0) a[i]<-1
    }
  }
  
  out<-1/out
  
  }


  return(out)

}
