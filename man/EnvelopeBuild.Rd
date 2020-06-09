\name{EnvelopeBuild_c}
%\alias{glmbenvelope}
\alias{EnvelopeBuild_c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds Envelope function for simulation
}
\description{
Constructs Enveloping function for simulation
}
\usage{
EnvelopeBuild_c(bStar,A,y,x,mu,P,alpha,wt,family="binomial",link = "logit", 
Gridtype=2,n = 1, sortgrid = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bStar}{Point at which envelope should be centered (typically posterior mode)
}
  \item{A}{Diagonal precision matrix for log-likelihood function associated with model in standard form
}
%  \item{f2}{Negative log-posterior function}
%  \item{f3}{Gradient function for Negative log-posterior function  }
%  \item{f5}{Vectorized version of Negative log-posterior function}
%  \item{f6}{Vectorized version of Gradient function for Negative log-posterior %function}
  \item{y}{
a vector of observations of length m
}
  \item{x}{
a design matrix of dimension m*p
}
  \item{mu}{
a vector giving the prior means of the variables.
}
  \item{P}{
a positive-definite symmetric matrix specifying the prior precision matrix of the variables.
}
  \item{alpha}{offset vector 
}
  \item{wt}{
a vector of weights
}
  \item{family}{Family for the envelope. Can take on values binomial, quasibinomial, poisson, quasibinomial, and Gamma}
  \item{link}{Link function for the envelope. Valid links for the binomial and quasibinomial families are logit, probit, and cloglog. The valid value for the poisson, quasipoisson, and the Gamma families is the log}
  
  \item{Gridtype}{
an optional argument specifying the method used to determine number of likelihood subgradient 
densities used to construct the enveloping function.
}
  \item{n}{
number of draws to generate from posterior density. Used here to help determine the size of the grid (Gridtype 1 and 2 only)
}
  \item{sortgrid}{ 
an optional Logical argument determining whether the final envelope should be sorted in descending order based on probability for each part of the envelope}
}
\details{
This function identfies tangency points and intervals for the enveloping function, evaluates the log-likelihood function and gradients for the log-likelihood at each point, and calculates additional constants.
}
\value{
The function returns a list consisting of the following components (the first six of which are matrics with number of rows equal to the number of components in the Grid and columns equal to the number of parameters):
 \item{GridIndex}{A matrix containing information on how each dimension should be sampled (1 means left tail of a restricted normal, 2 center, 3 right tail, and 4 the entire line) }
  \item{thetabars }{A matrix containing the points of tangencies associated with each component of the grid}
  \item{cbars}{A matrix containing the gradients for the negative log-likelihood at each tangency}
  \item{logU}{A matrix containing the log of the cummulative probability associated with each dimension}
  \item{logrt}{A matrix containing the log of the probability associated with the right tail (i.e. that to the right of the lower bound)}
  \item{loglt}{A matrix containing the log of the probability associated with the left tail (i.e., that to the left of the upper bound)}
  \item{LLconst}{A vector containing constant for each component of the grid used during the accept-reject procedure}  
  \item{logP}{A matrix containing log-probabilities related to the components of the grid }  
  \item{PLSD}{A vector containing the probability of each component in the Grid}  


%% ...
}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rglmb}}, \code{\link{expand.grid}}, \code{\link{EnvelopeOpt}}
}
\examples{

data(menarche)

summary(menarche)
plot(Menarche/Total ~ Age, data=menarche)

Age2=menarche$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche$Menarche/menarche$Total
wt=menarche$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 1 times as large as point estimate 

famfunc<-glmbfamfunc(binomial(logit))

f1<-famfunc$f1
f2<-famfunc$f2
f3<-famfunc$f3
f5<-famfunc$f5
f6<-famfunc$f6

dispersion2<-as.numeric(1.0)
start <- mu
offset2=rep(as.numeric(0.0),length(y))
P=solve(V1)
n=1000

## Appears that the type for some of these arguments are important/problematic

outlist<-rnnorm_reg_cpp(n=as.integer(n),y=as.vector(y),
x=as.matrix(x),mu=as.vector(mu),P=as.matrix(P),
offset2=as.vector(offset2),wt=as.vector(wt),dispersion=as.numeric(dispersion2),
famfunc=famfunc,f1=f1,f2=f2,f3=f3,start=as.vector(start),family="binomial",
link="logit",Gridtype=as.integer(3))

### This allows use of the rglmb summary function 
### add interface for glmbsim_NGauss_cpp later

#outlist$call<-match.call()
colnames(outlist$coefficients)<-colnames(x)
class(outlist)<-c(outlist$class,"rglmb")
#summary(outlist)
Env1=outlist$Envelope


###### Adjust weight for dispersion

wt2=wt/dispersion2

######################### Shift mean vector to offset so that adjusted model has 0 mean

alpha=x\%*\%as.vector(mu)+offset2
mu2=0*as.vector(mu)
P2=P
x2=x


#####  Optimization step to find posterior mode and associated Precision

parin=start-mu

opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x),mu=as.vector(mu2),
      P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
      method="BFGS",hessian=TRUE
      )

bstar=opt_out$par  ## Posterior mode for adjusted model
bstar
bstar+as.vector(mu)  # mode for actual model
A1=opt_out$hessian # Approximate Precision at mode

## Standardize Model

Standard_Mod=glmb_Standardize_Model(y=as.vector(y), x=as.matrix(x),P=as.matrix(P),
                       bstar=as.matrix(bstar,ncol=1), A1=as.matrix(A1))

bstar2=Standard_Mod$bstar2  
A=Standard_Mod$A
x2=Standard_Mod$x2
mu2=Standard_Mod$mu2
P2=Standard_Mod$P2
L2Inv=Standard_Mod$L2Inv
L3Inv=Standard_Mod$L3Inv

Env2=EnvelopeBuild_c(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
    as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),as.vector(wt2),
    family="binomial",link="logit",Gridtype=as.integer(3), n=as.integer(n),sortgrid=TRUE)

## These now seem to match

Env1
Env2


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
