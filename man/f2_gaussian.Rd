\name{f2_gaussian}
\alias{f2_gaussian}
\alias{f2_poisson}
\alias{f2_binomial_logit}
\alias{f2_binomial_probit}
\alias{f2_binomial_cloglog}
\alias{f2_gamma}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Negative Log-Posterior Functions
}
\description{
Negative Log-Posterior Functions
}
\usage{
f2_gaussian(b, y, x, mu, P, alpha, wt)
f2_poisson(b, y, x, mu, P, alpha, wt,progbar)
f2_binomial_logit(b, y, x, mu, P, alpha, wt,progbar)
f2_binomial_probit(b, y, x, mu, P, alpha, wt,progbar)
f2_binomial_cloglog(b, y, x, mu, P, alpha, wt,progbar)
f2_gamma(b, y, x, mu, P, alpha, wt,progbar)
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
A prior precicion matrix of dimension p*p}
  \item{alpha}{
A vector with offsets
}
  \item{wt}{
A vector of weights
}

\item{progbar}{
A dummy for indicating whether a progress-bar should be displayed during the call
}

}
\details{
Negative Log-Posterior Functions
}
\value{
A vector where each element represents the negative log-posterior for
the corresponding vector of coefficients
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
    .Call("glmbayes2_f2_gaussian", PACKAGE = "glmbayes2", b, 
        y, x, mu, P, alpha, wt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword perline
