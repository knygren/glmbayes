
optgrid<-function(a1,n){

a1rank<-rank(1/(1+a1))
l1<-length(a1)

dimcount<-matrix(0,(l1+1),l1)
scaleest<-matrix(0,(l1+1),l1)
intest<-c(1:(l1+1))
slopeest<-c(1:(l1+1))

dimcount[1,]<-diag(diag(l1))
scaleest[1,]<-1+a1
slopeest[1]<-prod(scaleest[1,])

for(i in 2:(l1+1)){
dimcount[i,]<-dimcount[i-1,]
scaleest[i,]<-scaleest[i-1,]
for(j in 1:l1){
if(a1rank[j]==i-1){ 
	dimcount[i,j]<-3
	scaleest[i,j]<-2/sqrt(pi) 
			}
		}
intest[i]<-3^(i-1)
slopeest[i]<-prod(scaleest[i,])
}
evalest<-intest+n*slopeest
minindex<-0
for(j in 1:(l1+1)){if(evalest[j]==min(evalest)){minindex<-j}}


dimcount[minindex,]

}