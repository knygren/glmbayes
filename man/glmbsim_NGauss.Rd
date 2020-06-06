\name{glmbsim_NGauss_cpp}
\alias{glmbsim_NGauss_cpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates From Non-Gaussian Models
}
\description{
Simulates From Non-Gaussian Models
}
\usage{
glmbsim_NGauss_cpp(n,y,x,mu,P,offset2,wt,dispersion,famfunc,f1,f2,f3,start,
family,link,Gridtype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of draws to generate from posterior density.
}
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
  \item{offset2}{
this can be used to specify an a priori known component to be included in the linear predictor during
 fitting. This should be NULL or a numeric vector of length equal to the number of cases.
}
  \item{wt}{
a vector of weights to be used in the fitting process. 
}
  \item{dispersion}{
dispersion parameter
}

  \item{famfunc}{
family function used
}
  \item{f1}{
f1 function (negative log-likelihood)}
  \item{f2}{
f2 function (negative log-posterior)
}
  \item{f3}{
f3 function (gradient vector)
}
  \item{start}{
an optional argument providing starting values for the posterior mode optimization.}
  \item{family}{
family
}
  \item{link}{
link function
}
  \item{Gridtype}{
an optional argument specifying the method used to determine number of likelihood subgradient 
densities used to construct the enveloping function.
}


}
\details{
Simulates from Non-Gaussian regression functions
}
\value{
A list with the following components
\item{coefficients}{an n by length(mu) matrix with one sample in each row}
\item{PostMode}{a vector of length(mu) with the estimated posterior mode coefficients}
\item{Prior}{A list with two components. The first being the prior mean vector and the second the prior precision matrix}
\item{iters}{an n by 1 matrix giving the number of candidates generated before acceptance for each sample.}
\item{famfunc}{an object of class famfunc}
\item{Envelope}{an object of class envelope}
\item{dispersion}{the dispersion parameter used in the model}
\item{loglike}{an n by 1 matrix containing the negative loglikelihood for each sample.}


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

## Appears that the type for some of these arguments are important/problematic

outlist<-glmbsim_NGauss_cpp(n=as.integer(n),y=as.vector(y),
x=as.matrix(x),mu=as.vector(mu),P=as.matrix(P),
offset2=as.vector(offset2),wt=as.vector(wt),dispersion=as.numeric(dispersion2),
famfunc=famfunc,f1=f1,f2=f2,f3=f3,start=as.vector(start),family="binomial",link="logit",
Gridtype=as.integer(3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
