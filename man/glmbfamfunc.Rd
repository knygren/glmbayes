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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (family) 
UseMethod("glmbfamfunc")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
