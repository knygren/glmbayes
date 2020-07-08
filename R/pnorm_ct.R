#' The Central Normal Distribution
#'
#' Distribution function and random generation for the center (between a lower and an upper bound)
#' of the normal distribution with mean equal to mu and standard deviation equal to sigma.
#' 
#' @name 
#' Normal_ct
#' @aliases
#' Normal_ct
#' pnorm_ct
#' rnorm_ct
#' @param a lower bound
#' @param b upper bound
#' @param n number of draws to generate. If \code{length(n) > 1}, 
#' the length is taken to be the number required
#' @param lgrt log of the distribution function between the 
#' lower bound and infinity
#' @param lglt log of the distribution function between negative 
#' infinity and the upper bound
#' @param mu mean parameter
#' @param sigma standard deviation
#' @param log.p Logical argument. If \code{TRUE}, the log probability is provided
#' @param Diff Logical argument. If \code{TRUE} the second parameter is the difference between the lower and upper bound
#' @return For \code{pnorm_ct}, vector of length equal to length of \code{a} and for
#' \code{rnorm_ct}, a vector with length determined by \code{n} containing draws from 
#' the center of the normal distribution.
#' @details  The distribution function pnorm_ct finds the probability of the 
#' center of a normal density (the probability of the area between a lower 
#' bound a and an upper bound b) while the random number generator rnorm_ct 
#' samples from a restricted normal density where lgrt is the log of the 
#' distribution between the lower bound and infinity and lglt is the log of 
#' the distribution function between negative infinity and the upper bound. 
#' The sum of the exponentiated values for the two (exp(lgrt)+exp(lglt)) must 
#' sum to more than 1.
#' 
#' These functions are mainly used to handle cases where the differences 
#' between the upper and lower bounds \code{b-a} are small. In such cases,
#' using \code{pnorm(b)-pnorm(a)} may result in 0 being returned even when the 
#' difference is supposed to be positive.
#' @keywords internal 
#' @example inst/examples/Ex_Normal_ct.R
#' @rdname Normal_ct
#' @order 1
#' @export

pnorm_ct<-function(a=-Inf,b=Inf,mu=0,sigma=1,log.p=TRUE,Diff=FALSE){

if(is.vector(a)==FALSE||is.vector(b)==FALSE){stop("Arguments a and b must be vectors")}
if(is.numeric(a)==FALSE||is.numeric(b)==FALSE||is.numeric(mu)==FALSE||is.numeric(sigma)==FALSE){stop("Non-numeric argument to mathematical function")}
if(is.logical(log.p)==FALSE||is.logical(Diff)==FALSE){stop("Arguments log.P and Diff must be logical")}
if(length(a)!=length(b)){stop("a and b must be equal length vectors")}
if(length(mu)>1||length(sigma)>1){stop("mu and sigma must have length 1")}
if(sigma<=0){stop("sigma must be positive")}

# Verify arguments are numeric vectors

l1<-length(a)

output<-0*c(1:length(a))

for(i in 1:length(a)){

# Diff == True means that second argument is a difference
# This is useful when the difference between a and b is small  
  
if(Diff==FALSE){
b2<-b[i]
if(length(a)>1){Diff2<-b[i]-a[i]}
if(length(a)==1){Diff2<-b[i]-a}


if(Diff2<0){stop("Argument b must be greater than argument a")}
}
if(Diff==TRUE){
b2<-a[i]+b[i]
Diff2<-b[i]
if(Diff2<=0){stop("Argument b Must be positive")}
}

g1<-mu-a
g2<-b2-mu

if(log.p==FALSE){
if(g1>=g2){output[i]<-pnorm(q=b2,mean=mu,sd=sigma)*(1-exp(pnorm(q=a,mean=mu,sd=sigma,log.p=TRUE)-pnorm(q=b2,mean=mu,sd=sigma,log.p=TRUE))) }
if(g1<g2){output[i]<-pnorm(q=a,mean=mu,sd=sigma,lower.tail=FALSE)*(1-exp(pnorm(q=b2,mean=mu,sd=sigma,log.p=TRUE,lower.tail=FALSE)-pnorm(q=a,mean=mu,sd=sigma,log.p=TRUE,lower.tail=FALSE)))}
if(g1==Inf && g2==Inf){output[i]<-1} 
if(output[i]==0){output[i]<-Diff2*dnorm(x=0.5*(b2+a),mean=mu,sd=sigma)  
if(output[i]==0){stop("Function ctpnorm evaluates to 0")}
}
}

if(log.p==TRUE){
output[i]<-NA
output[i]<-log(pnorm(q=b,mean=mu,sd=sigma)-pnorm(q=a,mean=mu,sd=sigma))
}
}
output
}




#if(g1>=g2){output[i]<-pnorm(q=b2,mean=mu,sd=sigma,log.p=TRUE)*log((1-exp(pnorm(q=a,mean=mu,sd=sigma,log.p=TRUE)-pnorm(q=b2,mean=mu,sd=sigma,log.p=TRUE)))) }
#if(g1<g2){output[i]<-pnorm(q=a,mean=mu,sd=sigma,lower.tail=FALSE,log.p=TRUE)*log((1-exp(pnorm(q=b2,mean=mu,sd=sigma,log.p=TRUE,lower.tail=FALSE)-pnorm(q=a,mean=mu,sd=sigma,log.p=TRUE,lower.tail=FALSE))))}
#if(g1==Inf && g2==Inf){output[i]<-0} }  
#if(output[i]==-Inf){
#output[i]<-log(Diff2)+dnorm(x=0.5*(b2+a),mean=mu,sd=sigma,log=TRUE)
#if(output[i]==-Inf){stop("Function ctpnorm evaluates to -Inf")}
#}
