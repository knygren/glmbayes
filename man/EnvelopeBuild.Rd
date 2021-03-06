% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnvBuild.R
\name{EnvelopeBuild}
\alias{EnvelopeBuild}
\title{Builds Envelope function for simulation}
\usage{
EnvelopeBuild(
  bStar,
  A,
  y,
  x,
  mu,
  P,
  alpha,
  wt,
  family = "binomial",
  link = "logit",
  Gridtype = 2L,
  n = 1L,
  sortgrid = FALSE
)
}
\arguments{
\item{bStar}{Point at which envelope should be centered (typically posterior mode}

\item{A}{Diagonal precision matrix for log-likelihood function associated with model in standard form}

\item{y}{a vector of observations of length \code{m}}

\item{x}{a design matrix of dimension \code{m * p}}

\item{mu}{a vector giving the prior means of the variables}

\item{P}{a positive-definite symmetric matrix specifying the prior precision matrix of the variables}

\item{alpha}{offset vector}

\item{wt}{a vector of weights}

\item{family}{Family for the envelope. Can take on values binomial, quasibinomial, poisson, quasibinomial, and Gamma}

\item{link}{Link function for the envelope. Valid links for the binomial and quasibinomial families are logit, probit, and cloglog. The valid value for the poisson, quasipoisson, and the Gamma families is the log}

\item{Gridtype}{an optional argument specifying the method used to determine number of likelihood subgradient 
densities used to construct the enveloping function}

\item{n}{number of draws to generate from posterior density. Used here to help determine the size of the grid (Gridtype 1 and 2 only)}

\item{sortgrid}{an optional Logical argument determining whether the final envelope should be sorted in descending order based on probability for each part of the envelope}
}
\value{
The function returns a list consisting of the following components (
the first six of which are matrices with the number of rows equal to the 
number of components in the Grid and columns equal to the number of parameters):
\item{GridIndex}{A matrix containing information on how each dimension should be 
sampled (1 means left tail of a restricted normal, 2 center, 3 right tail, and 4 the 
entire line)}
\item{thetabars}{A matrix containing the points of tangencies associated with each component of the grid}
\item{cbars}{A matrix containing the gradients for the negative log-likelihood at each tangency}
\item{logU}{A matrix containing the log of the cummulative probability associated with each dimension}
\item{logrt}{A matrix containing the log of the probability associated with the right tail (i.e. 
that to the right of the lower bound)}
\item{loglt}{A matrix containing the log of the probability associated with the left tail (i.e.,
that to the left of the upper bound)}
\item{LLconst}{A vector containing constant for each component of the grid used during the accept-reject procedure}
\item{logP}{A matrix containing log-probabilities related to the components of the grid}
\item{PLSD}{A vector containing the probability of each component in the Grid}
}
\description{
Builds an Enveloping function for simulation using a grid and tangencies for 
the posterior density
}
\details{
To construct an enveloping function, we follow the approach in Nygren and Nygren (2006)
which involves the following steps when a maximally sized grid is constructed 
(if the prior for some dimensions is relatively strong, this may not be needed)


1) For each dimension, a constant \code{omega_i} is found that depends on the 
corresponding diagonal element in the precision matrix.

 
2) Corresponding intervals \code{(thetastar_i-0.5 *omega,thetastar_i-0.5 *omega)} 
 are constructed around the posterior mode thetastar for each dimension
 
 
3) The mode as well as the points \code{thetastar_i-omega_i} and \code{thetastar_i+omega_i}
are selected as the components of the points at which tangencies will be found for each of 
the dimensions.  


4) A Grid is constructed with all possible combinations of points and 
negative log-likelihood and gradient for the negative log-likelihood are evaluated (see the
EnvelopeBuild_c.cpp function source code for details)


5) The \code{\link{Set_Grid}} function is called in order to evaluate the log of the density associated
with each of the resulting restricted multivariate normals by evaluating the differences between the cummulative 
density for each dimension between its lower and upper bound.


6)  The Set_LogP is called in order to help set the probabilities with which each of the components
of the grid should be sampled (see remark 6 in Nygren and Nygren (2006)).


Any constants needed by the sampling are added to a list and returned by the function.
}
\examples{
data(menarche2)

summary(menarche2)
plot(Menarche/Total ~ Age, data=menarche2)

Age2=menarche2$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche2$Menarche/menarche2$Total
wt=menarche2$Total

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

Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),as.vector(wt2),
family="binomial",link="logit",Gridtype=as.integer(3), n=as.integer(n),
sortgrid=TRUE)

## These now seem to match

Env2
}
\keyword{internal}
