#' Calculates DIC and Deviance Information
#'
#' Caculates DIC and Deviance Information for a fitted model.
#' @param coefficients A matrix with coefficients from the rglmb function
#' @param y a vector of observations of length m
#' @param x a design matrix of dimension m*p
#' @param alpha an offset parameter or vector
#' @param f1 Function with signature `f1(b, y, x, alpha, wt)` returning negative log-likelihood (scalar).
#' @param f4 Function with signature `f4(b, y, x, alpha, wt, dispersion)` returning deviance (scalar).
#' @param wt a vector of weights
#' @param dispersion dispersion parameter
#' 
#' @details Calculates DIC and Deviance Information
#' @return A list with the following components
#' \item{Deviance}{A \code{n * 1} matrix with the deviance for each draw}
#' \item{Dbar}{Mean for negative 2 times negative log-likelihood}
#' \item{Dthetabar}{Negative 2 times log-likelihood evaluated at mean parameters}
#' \item{pD}{Effective number of parameters}
#' \item{DIC}{DIC statistic}
#' @example inst/examples/Ex_glmbdic.R
#' @noRd 

DIC_Info<-function(coefficients,y,x,alpha=0,f1,f4,wt=1,dispersion=1){
  
  l1<-length(coefficients[1,])
  l2<-length(coefficients[,1])
  
  D<-matrix(0,nrow=l2,ncol=1)
  D2<-matrix(0,nrow=l2,ncol=1)
  
  if(length(dispersion)==1){
    for(i in 1:l2){
      b<-as.vector(coefficients[i,])
      D[i,1]<-f4(b=b,y=y,x=x,alpha=alpha,wt=wt,dispersion=dispersion)
      
      D2[i,1]<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt)
    }
    
    Dbar<-mean(D2)
    
    b<-colMeans(coefficients)
    
    Dthetabar<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt)
    
  }
  
  if(length(dispersion)>1){
    for(i in 1:l2){
      b<-as.vector(coefficients[i,])
      D[i,1]<-f4(b=b,y=y,x=x,alpha=alpha,wt=wt,dispersion=dispersion[i])
      
      D2[i,1]<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt/dispersion[i])
    }
    
    Dbar<-mean(D2)
    
    b<-colMeans(coefficients)
    dispbar<-mean(dispersion)
    Dthetabar<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt/dispbar)
    
    
  }
  
  
  pD<-Dbar-Dthetabar
  DIC<-pD+Dbar
  list(Deviance=D,Dbar=Dbar,Dthetabar=Dthetabar,pD=pD,DIC=DIC)
}
