\name{glmb_Standardize_Model}
\alias{glmb_Standardize_Model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standardize A Non-Gaussian Model}
\description{
Standardizes a Non-Gaussian Model prior to Envelope Creation
}
\usage{
glmb_Standardize_Model(y, x,P, bstar, A1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
a vector of observations of length m
}
  \item{x}{
a design matrix of dimension m*p
}
  \item{P}{
a positive-definite symmetric matrix specifying the prior precision matrix of the variables.
}
  \item{bstar}{a matrix containing the posterior mode from an optimization step
}
  \item{A1}{a matrix containing the posterior precision matrix at the posterior mode
}


}
\details{
Standardizes a Non-Gaussian Model prior to Envelope Creation
}
\value{
A list with the following components
\item{bstar2}{Standardized Posterior Mode}
\item{A}{Standardized Data Precision Matrix}
\item{x2}{Standardized Design Matrix}
\item{mu2}{Standardized Prior Mean vector}
\item{P2}{Standardized Precision Matrix Added to log-likelihood}
\item{L2Inv}{A matrix used when undoing the first step in standardization}
\item{L3Inv}{A matrix used when undoing the second step in standardization}

}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glmb}}
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
    family="binomial",link="logit",Gridtype=as.integer(3), n=as.integer(n),sortgrid=TRUE)

## These now seem to match

Env2


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
