
glmbdic<-function(coefficients,y,x,alpha=0,f1,f4,wt=1,dispersion=1){

l1<-length(coefficients[1,])
l2<-length(coefficients[,1])

D<-matrix(0,nrow=l2,ncol=1)
D2<-matrix(0,nrow=l2,ncol=1)

if(length(dispersion)==1){
for(i in 1:l2){
b<-as.vector(coefficients[i,])
# No offset! [Edit to multiply by dispersion to get residual deviance]

D[i,1]<-dispersion*f4(b=b,y=y,x=x,alpha=alpha,wt=wt/dispersion,dispersion=dispersion)
D2[i,1]<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt)
        }

Dbar<-mean(D2)

b<-colMeans(coefficients)

Dthetabar<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt)

}

if(length(dispersion)>1){
  for(i in 1:l2){
    b<-as.vector(coefficients[i,])
    D[i,1]<-f4(b=b,y=y,x=x,alpha=alpha,wt=wt/dispersion[i],dispersion=dispersion[i])
    D2[i,1]<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt/dispersion[i])
  }

  Dbar<-mean(D2)
  
  b<-colMeans(coefficients)
  # Use the mean for the dispersion in Dthetabar calculation
  # using mean for 1/dispersion (precision) would change the DIC
  dispbar<-mean(dispersion)
  Dthetabar<-2*f1(b=b,y=y,x=x,alpha=alpha,wt=wt/dispbar)
  
  
}


pD<-Dbar-Dthetabar
DIC<-pD+Dbar
list(Deviance=D,Dbar=Dbar,Dthetabar=Dthetabar,pD=pD,DIC=DIC)
}


