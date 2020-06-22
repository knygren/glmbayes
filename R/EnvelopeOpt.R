#' Optimizes Envelope function for simulation
#'
#' Optimizes the size of the grid to try to limit the combined time of
#' the envelope construction and the simulation phase.
#' @param a1 Diagonal elements of data precision matrix for a model in standard form
#' @param n  Number of draws to generate
#' @details This function attempts to find a computationally optimal 
#' gridsize by using information on the strenght of the prior and the 
#' number of iterations desired. Generally, more data (i.e., larger values 
#' for the diagonal elements of the precision matrix) will require a larger grid.
#' The same also holds when the number of desired draws is higher 
#' (as the setup costs associated with the larger grid is offset by the 
#' savings in the number of candidates per sample).
#' @return A vector containing information on how many component each 
#' dimension should be split into.
#' @seealso \code{\link{rglmb}}, \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeSort}} 
#' @example inst/examples/Ex_EnvelopeOpt.R
#' @export


EnvelopeOpt<-function(a1,n){

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