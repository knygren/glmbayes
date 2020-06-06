\name{glmbfamfunc}
\alias{glmbfamfunc}
%\alias{glmbfamfunc.default}
\alias{glmbfamfuncset}
\alias{print.glmbfamfunc}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns family functions used during simulation and post processing
}
\description{
This function takes as input a \code{\link{family}} object and returns a set of functions that are used to simulate and summarize models using the \code{\link{glmb}}, and \code{\link{rglmb}} functions.
}
\usage{
glmbfamfunc(family)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{family}{
an object of class \code{\link{family}}
}
}
\details{
This function takes as input a family and returns a set of functions related to the family.
}
\value{
A list with the following components
  \item{f1 }{Negative log-likelihood function}
  \item{f2}{Negative log-posterior function}
  \item{f3}{Gradient function for negative log-posterior function}
  \item{f4}{Deviance function}
  \item{f5}{Vectorized version of Negative log-posterior function}
  \item{f6}{Vectorized version of Gradient for log-posterior function}
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
\code{\link{family}}, \code{\link{rglmb}}
}
\examples{
data(menarche)

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


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
