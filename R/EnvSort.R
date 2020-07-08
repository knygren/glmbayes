#' Sorts Envelope function for simulation
#'
#' Sorts Enveloping function for simulation.
#' how frequently each component of the resulting grid should be sampled
#' during simulation.
#' @param l1 dimension for model (number of independent variables in X matrix)
#' @param l2 dimension for Envelope (number of components)
#' @param GIndex matrix containing information on how each dimension should be sampled (1 means left tail of a restricted normal, 2 center, 3 right tail, and 4 the entire line)
#' @param G3 A matrix containing the points of tangencies associated with each component of the grid
#' @param cbars A matrix containing the gradients for the negative log-likelihood at each tangency
#' @param logU A matrix containing the log of the cummulative probability associated with each dimension
#' @param logrt A matrix containing the log of the probability associated with the right tail (i.e. that to the right of the lower bound)
#' @param loglt A matrix containing the log of the probability associated with the left tail (i.e., that to the left of the upper bound)
#' @param logP A matrix containing log-probabilities related to the components of the grid
#' @param LLconst A vector containing constant for each component of the grid used during the accept-reject procedure
#' @param PLSD A vector containing the probability of each component in the Grid
#' @param a1 A vector containing the diagonal of the standardized precision matrix
#' @return The function returns a list consisting of the following components (the first six of which are matrics with number of rows equal to the number of components in the Grid and columns equal to the number of parameters):
#' \item{GridIndex}{A matrix containing information on how each dimension should be sampled (1 means left tail of a restricted normal, 2 center, 3 right tail, and 4 the entire line)}
#' \item{thetabars}{A matrix containing the points of tangencies associated with each component of the grid}
#' \item{cbars}{A matrix containing the gradients for the negative log-likelihood at each tangency}
#' \item{logU}{A matrix containing the log of the cummulative probability associated with each dimension}
#' \item{logrt}{A matrix containing the log of the probability associated with the right tail (i.e. that to the right of the lower bound)}
#' \item{loglt}{A matrix containing the log of the probability associated with the left tail (i.e., that to the left of the upper bound)}
#' \item{LLconst}{A vector containing constant for each component of the grid used during the accept-reject procedure}
#' \item{logP}{A matrix containing log-probabilities related to the components of the grid}
#' \item{PLSD}{A vector containing the probability of each component in the Grid}
#' @details This function sorts the envelope in descending order based on the 
#' probability associated with each component in the Grid. Sorting helps 
#' speed up simulation once the envelope is constructed.
#' @example inst/examples/Ex_EnvelopeSort.R
#' @export
#' @keywords internal



EnvelopeSort<-function(l1,l2,GIndex,G3,cbars,logU,logrt,loglt,logP,LLconst,PLSD,a1){
  
  Envelope<-data.frame(GIndex=GIndex,G3=G3,cbars=cbars,logU=logU,logrt=logrt,loglt=loglt,logP=logP[,1],LLconst=LLconst,PLSD=PLSD)
  Envelope<-Envelope[order(-Envelope$PLSD),]
  
  if(l1==1){
    outlist<-list(GridIndex=Envelope[1:l2,1:l1],thetabars=Envelope[1:l2,(l1+1):(2*l1)],
                  cbars=matrix(Envelope[1:l2,(2*l1+1):(3*l1)],nrow=l2,ncol=l1),logU=Envelope[1:l2,(3*l1+1):(4*l1)],
                  logrt=matrix(Envelope[1:l2,(4*l1+1):(5*l1)],nrow=l2,ncol=l1),
                  loglt=matrix(Envelope[1:l2,(5*l1+1):(6*l1)],nrow=l2,ncol=l1),
                  LLconst=Envelope$LLconst,logP=Envelope$logP,PLSD=Envelope$PLSD,a1=a1)
  }
  
  else{
    outlist<-list(GridIndex=Envelope[1:l2,1:l1],thetabars=Envelope[1:l2,(l1+1):(2*l1)],
                  cbars=as.matrix(Envelope[1:l2,(2*l1+1):(3*l1)]),logU=Envelope[1:l2,(3*l1+1):(4*l1)],
                  logrt=as.matrix(Envelope[1:l2,(4*l1+1):(5*l1)]),
                  loglt=as.matrix(Envelope[1:l2,(5*l1+1):(6*l1)]),LLconst=Envelope$LLconst,logP=Envelope$logP,PLSD=Envelope$PLSD,a1=a1)
    
  }
  
  return(outlist)
  
}
