\name{ctpnorm}
\alias{ctpnorm}
\alias{ctrnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Central Normal Distribution
}
\description{Distribution function and random generation for the center (between a lower and an upper bound) of the normal distribution with mean equal to mu and standard deviation equal to sigma.
}
\usage{
ctpnorm(a = -Inf, b = Inf, mu = 0, sigma = 1, log.p = TRUE, Diff = FALSE)
ctrnorm(n,lgrt,lglt,mu=0,sigma=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
lower bound}
  \item{b}{
upper bound}
  \item{n}{
number of draws to generate}
  \item{lgrt}{
log of the distribution function between the lower bound and infinity}
  \item{lglt}{
log of the distribution function between negative infinity and the upper bound}
  \item{mu}{
mean parameter}
  \item{sigma}{
standard deviation
}
  \item{log.p}{
Logical argument. If True, the log probability is provided}
  \item{Diff}{
If True the second parameter is the difference between the lower and upper bound
}
}
\details{

The distribution function ctpnorm finds the probability of the center of a normal density (the probability of the area between a lower bound \code{a} and an upper bound \code{b}) while the random number generator ctrnorm samples from a restricted normal density where lgrt is the log of the distribution between the lower bound and infinity and lglt is the log of the distribution function between negative infinity and the upper bound. The sum of the exponentiated values for the two (exp(lgrt)+exp(lglt)) must sum to more than 1.

}
\value{
A vector of length equal to length of a}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}


\seealso{
\code{\link{glmbsim_cpp}}
}
\examples{

ctpnorm(0.2,0.4)
exp(ctpnorm(0.2,0.4))
ctpnorm(0.2,0.4,log.p=FALSE)
log(ctpnorm(0.2,0.4,log.p=FALSE))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
