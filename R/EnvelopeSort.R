EnvelopeSort<-function(l1,l2,GIndex,G3,cbars,logU,logrt,loglt,logP,LLconst,PLSD,a1){
  
  Envelope<-data.frame(GIndex=GIndex,G3=G3,cbars=cbars,logU=logU,logrt=logrt,loglt=loglt,logP=logP[,1],LLconst=LLconst,PLSD=PLSD)
  Envelope<-Envelope[order(-Envelope$PLSD),]
  
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
