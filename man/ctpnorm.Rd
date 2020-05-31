\name{ctpnorm}
\alias{ctpnorm}
\alias{ctrnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Finds the probability of the center of a normal density (probability of area between lower and upper bound)}
\description{
Finds the probability of the center of a normal density (probability of area between lower and upper bound)}
\usage{
ctpnorm(a = -Inf, b = Inf, mu = 0, sigma = 1, log.p = TRUE, Diff = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
lower bound}
  \item{b}{
upper bound}
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
Finds the probability of the center of a normal density (probability of area between lower and upper bound)
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
