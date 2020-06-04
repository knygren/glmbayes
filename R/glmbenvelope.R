
glmbenvelope<-function(bStar,A,f2,f3,f5,f6,y,x,mu,P,alpha,wt=1,Gridtype=3,n=1){

a1<-diag(A)
l1<-length(a1)

# Omega, Grid, and Lint defined as on p.1149 in Nygren (2006)

omega<-as.matrix((sqrt(2)-exp(-1.20491-0.7321*sqrt(0.5+a1)))/sqrt(1+a1))
G1<-as.matrix(c(1,1,1))%*%t(bStar)+  as.matrix(c(-1,0,1))%*%t(omega)
Lint<-as.matrix(c(1,1))%*%t(bStar)+  as.matrix(c(-0.5,0.5))%*%t(omega)

G2<-list()
GIndex1<-list()
length(G2)<-l1
length(GIndex1)<-l1

if(Gridtype==2){
gridindex<-optgrid(a1,n=n)

}


for(i in 1:l1){
if(Gridtype==1){
if((1+a1[i])<=(2/sqrt(pi))){G2[i]<-list(G1[2,i])
				    GIndex1[i]<-4}
if((1+a1[i])>(2/sqrt(pi))){
				 G2[i]<-list(G1[,i])
                         GIndex1[i]<-list(c(1,2,3))}
		  }

if(Gridtype==2){
if(gridindex[i]==1){G2[i]<-list(G1[2,i])
				    GIndex1[i]<-4}
if(gridindex[i]==3){
				 G2[i]<-list(G1[,i])
                         GIndex1[i]<-list(c(1,2,3))}
		  }



if(Gridtype==3){
G2[i]<-list(G1[,i])
GIndex1[i]<-list(c(1,2,3))
}
if(Gridtype==4){
G2[i]<-list(G1[2,i])
GIndex1[i]<-4}
}

G3<-expand.grid(G2)
GIndex<-expand.grid(GIndex1)


G3<-as.matrix(G3)
GIndex<-as.matrix(GIndex)
G3<-as.matrix(G3)

l2<-length(GIndex[,1])


cbars<-matrix(0,nrow=l2,ncol=l1)
NegLL<-matrix(0,nrow=l2,ncol=1)
Up<-matrix(0,nrow=l2,ncol=l1)
Down<-matrix(0,nrow=l2,ncol=l1)
logP<-matrix(0,nrow=l2,ncol=2)
logU<-matrix(0,nrow=l2,ncol=l1)
loglt<-matrix(0,nrow=l2,ncol=l1)
logrt<-matrix(0,nrow=l2,ncol=l1)


PLSD<-matrix(0,nrow=l2,ncol=1)
LLconst<-matrix(0,nrow=l2,ncol=1)




cbars<-f6(b=t(G3),y=y,x=x,mu=mu,P=P,alpha=alpha,wt=wt)
NegLL<-f5(b=t(G3),y=y,x=x,mu=mu,P=P,alpha=alpha,wt=wt)



# New Set Grid starts here


outloop<-Set_Grid(GIndex,cbars,Lint)


Down<-outloop$Down
Up<-outloop$Up
loglt<-outloop$lglt
logrt<-outloop$lgrt
logct<-outloop$logct
logU<-outloop$logU
logP<-outloop$logP



logout<-setlogP(logP,NegLL,cbars,G3)

logP<-logout$logP
LLconst<-logout$LLconst



maxlogP<-max(logP[,2])
PLSD<-exp(logP[,2]-maxlogP)
sumP<-sum(PLSD)
PLSD<-PLSD/sumP

PLSD<-as.matrix(PLSD)



Envelope<-data.frame(GIndex=GIndex,G3=G3,cbars=cbars,logU=logU,logrt=logrt,loglt=loglt,logP=logP[,1],LLconst=LLconst,PLSD=PLSD)

Envelope=Envelope[order(-Envelope$PLSD),]





if(l1==1){
outlist<-list(GridIndex=Envelope[1:l2,1:l1],thetabars=Envelope[1:l2,(l1+1):(2*l1)],
cbars=matrix(Envelope[1:l2,(2*l1+1):(3*l1)],nrow=l2,ncol=l1),logU=Envelope[1:l2,(3*l1+1):(4*l1)],
logrt=matrix(Envelope[1:l2,(4*l1+1):(5*l1)],nrow=l2,ncol=l1),
loglt=matrix(Envelope[1:l2,(5*l1+1):(6*l1)],nrow=l2,ncol=l1),LLconst=Envelope$LLconst,logP=Envelope$logP,PLSD=Envelope$PLSD,a1=a1)
}

else{
  outlist<-list(GridIndex=Envelope[1:l2,1:l1],thetabars=Envelope[1:l2,(l1+1):(2*l1)],
                cbars=as.matrix(Envelope[1:l2,(2*l1+1):(3*l1)]),logU=Envelope[1:l2,(3*l1+1):(4*l1)],
                logrt=as.matrix(Envelope[1:l2,(4*l1+1):(5*l1)]),
                loglt=as.matrix(Envelope[1:l2,(5*l1+1):(6*l1)]),LLconst=Envelope$LLconst,logP=Envelope$logP,PLSD=Envelope$PLSD,a1=a1)
  
}



return(outlist)

}



