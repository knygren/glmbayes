\name{f1_gaussian}
\alias{f1_gaussian}
\alias{f1_poisson}
\alias{f1_binomial_logit}
\alias{f1_binomial_probit}
\alias{f1_binomial_cloglog}
\alias{f1_gamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Provides a vectorized evaluation of the Negative Log-likelihood
}
\description{
Provides a vectorized evaluation of the Negative Log-likelihood
}
\usage{
f1_gaussian(b, y, x, alpha, wt)
f1_poisson(b, y, x, alpha, wt)
f1_binomial_logit(b, y, x, alpha, wt)
f1_binomial_probit(b, y, x, alpha, wt)
f1_binomial_cloglog(b, y, x, alpha, wt)
f1_gamma(b, y, x, alpha, wt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{
A matrix where each column represents a vector of coefficients at which the negative log-likelihood is to be evaluated
}
  \item{y}{
A vector of observations of length m
}
  \item{x}{
A design matrix of dimension m*p
}
  \item{alpha}{
A vector with offsets
}
  \item{wt}{
A vector of weights
}
}
\details{
Provides a vectorized evaluation of the Negative Log-likelihood
}
\value{
A vector with the negative log-likelihood for each vector of coefficients
}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}


\seealso{
\code{\link{glmbenvelope}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (b, y, x, alpha, wt) 
{
    .Call("glmbayes2_f1_gaussian", PACKAGE = "glmbayes2", b, 
        y, x, alpha, wt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
