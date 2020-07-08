#' @keywords internal
#' @rdname Normal_ct
#' @order 2
#' @export

rnorm_ct<-function(n,lgrt,lglt,mu=0,sigma=1)
{

if(is.numeric(lgrt)==FALSE||is.numeric(lglt)==FALSE||is.numeric(mu)==FALSE||is.numeric(sigma)==FALSE){stop("Non-numeric argument to mathematical function")}


pcheck<-exp(lgrt)+exp(lglt)
#if(1>=pcheck){stop("Combined probabilities of left and right tails must exceed 1")}
if(1>pcheck){stop("Combined probabilities of left and right tails must sum to at least 1")}


if(lgrt>=lglt){
		U<-runif(n,0,1)
		u1<-1-exp(lgrt)
		lgu1<-log(u1)
		if(lgu1>lglt) stop("lg(1-exp(lgrt)) must be less than lglt")
		check<-1-exp(lgu1-lglt)
#		if(is.nan(lgu1-lglt))
#			{ check<--Inf
#			  warning("lgu1-lglt evaluates to NaN")
#			}
#		if(check==0)stop("1-exp(lgu1-lglt) must exceed 0")
		lgU2<-log(U)+lglt+log(1-exp(lgu1-lglt))
#		lgU2<-log(U)+lglt+check

		lgU3<-lgU2+log(1+exp(lgu1-lgU2))

		if(is.nan(lgU3)) lgU3<--Inf

		out<-qnorm(lgU3,mean=mu,sd=sigma,lower.tail=TRUE,log.p=TRUE)
		}

if(lgrt<lglt){
		U<-runif(n,0,1)
		e1mu2<-1-exp(lglt)
		lg1mu2<-log(e1mu2)
		lgU2<-log(U)+lgrt+log(1-exp(lg1mu2-lgrt))
		lgU3<-lgU2+log(1+exp(lg1mu2-lgU2))
		out<-qnorm(lgU3,mean=mu,sd=sigma,lower.tail=FALSE,log.p=TRUE)
		}



out
}
