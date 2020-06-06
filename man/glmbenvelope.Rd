\name{glmbenvelope}
\alias{glmbenvelope}
\alias{glmbenvelope_c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructs Enveloping function for simulation
}
\description{
Constructs Enveloping function for simulation
}
\usage{
glmbenvelope_c(bStar,A,y,x,mu,P,alpha,wt,family="binomial",link = "logit", 
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
\code{\link{rglmb}}, \code{\link{expand.grid}}, \code{\link{optgrid}}, \code{\link{Set_Grid}}, \code{\link{setlogP}}
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
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 

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

#####

### Adjust weight for dispersion

wt2=wt/dispersion2

################################     #####################

# Shift mean vector to offset so that adjusted model has 0 mean

alpha=x\%*\%as.vector(mu)+offset2
mu2=0*as.vector(mu)
P2=P
x2=x

############################     ####################################

#####  Optimization step to find posterior mode and associated Precision

parin=start-mu

opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x),mu=as.vector(mu),
      P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
      method="BFGS",hessian=TRUE
      )

#opt_out
b2=opt_out$par  ## Posterior mode for adjusted model
b2+as.vector(mu)  # mode for actual model
A1=opt_out$hessian # Approximate Precision at mode


########################   #################################

# Find eigenvalues and standardize to model with ~ Identity precision at posterior mode

A1_eigen=eigen(A1)
D1=diag(A1_eigen$values)
L2=sqrt(D1)\%*\%t(A1_eigen$vectors)
L2Inv=A1_eigen$vectors\%*\%sqrt(solve(D1))

## Apply standardization

b3=L2\%*\%b2   
mu3=L2\%*\%as.vector(mu2)

x3=x2\%*\%L2Inv
P3=t(L2Inv)\%*\%P2\%*\%L2Inv

P3  # Standardized prior precision - "smaller" than identity matrix
L2Inv\%*\%b3+as.vector(mu)

#########################            ##########################

###  Find diagonal matrix epsilon so that P3-epsilon is still positive definite

### For this demo, simply use epsilon=0.5*P3
### In C++ code, searches until valid values (multiplying by 0.5 multiple times if needed)



P3Diag=diag(diag(P3))  #   diagonal part of P3
epsilon=0.5*P3Diag
P4=P3-epsilon



epsilon  ## Modified prior precision after P4 shifted to likelihood
P4       ## Precision for Multivariate normal term added to log-likelihood

#########################  Transform model again so that 
###  modified prior is the Standard multivariate normal

A3=diag(length(mu))-epsilon

A3_eigen=eigen(A3)
D2=diag(A3_eigen$values)
L3=sqrt(D2)\%*\%t(A3_eigen$vectors)
L3Inv=A3_eigen$vectors\%*\%sqrt(solve(D2))

b4=L3\%*\%b3 
mu4=L3\%*\%mu3 
x4=x3\%*\%L3Inv
A4=t(L3Inv)\%*\%A3\%*\%L3Inv #   Should be transformed data precision matrix
P5=t(L3Inv)\%*\%P4\%*\%L3Inv #   Should be precision matrix without epsilon
P6Temp=P5+diag(length(mu)) #  Should be precision matrix for posterior

L3Inv\%*\%L2Inv\%*\%b4+as.vector(mu)  ## Check that posterior mode still is the same

Envelope=glmbenvelope_c(as.vector(b4), as.matrix(A4),y, as.matrix(x4),
      as.matrix(mu4,ncol=1),as.matrix(P5),as.vector(alpha),as.vector(wt2),
      family="binomial",link="logit",Gridtype=as.integer(3), n=as.integer(n),sortgrid=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
