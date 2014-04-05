
glmbdic<-function(coefficients,y,x,f1,f4,wt=1,dispersion=1){

l1<-length(coefficients[1,])
l2<-length(coefficients[,1])

D<-matrix(0,nrow=l2,ncol=1)
D2<-matrix(0,nrow=l2,ncol=1)


for(i in 1:l2){
b<-as.vector(coefficients[i,])
D[i,1]<-f4(b=b,y=y,x=x,alpha=0,wt=wt,dispersion=dispersion)
D2[i,1]<-2*f1(b=b,y=y,x=x,alpha=0,wt=wt)


}
Dbar<-mean(D2)
#Dbar<-mean(D)

b<-colMeans(coefficients)
Dthetabar<-2*f1(b=b,y=y,x=x,alpha=0,wt=wt)
pD<-Dbar-Dthetabar
DIC<-pD+Dbar
list(Deviance=D,Dbar=Dbar,Dthetabar=Dthetabar,pD=pD,DIC=DIC)
}

