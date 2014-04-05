\name{f3_gaussian}
\alias{f3_gaussian}
\alias{f3_poisson}
\alias{f3_binomial_logit}
\alias{f3_binomial_probit}
\alias{f3_binomial_cloglog}
\alias{f3_gamma}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Negative Gradient function for posterior density
}
\description{
Negative Gradient function for posterior density
}
\usage{
f3_gaussian(b, y, x, mu, P, alpha, wt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{
A matrix where each column represents a vector of coefficients at which the negative log-posterior is to be evaluated
}
  \item{y}{
A vector of observations of length m
}
  \item{x}{
A design matrix of dimension m*p
}
  \item{mu}{
A prior mean vector of dimension p}
  \item{P}{
A prior precicion matrix of dimension p*p
}
  \item{alpha}{
A vector with offsets
}
  \item{wt}{
A vector of weights
}
}
\details{
Negative Gradient function for posterior density
}
\value{
A matrix where individual subvectors represents the negative gradient for the log-posterior function
}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.

}
\author{
Kjell Nygren
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glmbenvelope}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (b, y, x, mu, P, alpha, wt) 
{
    .Call("glmbayes2_f3_gaussian", PACKAGE = "glmbayes2", b, 
        y, x, mu, P, alpha, wt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
