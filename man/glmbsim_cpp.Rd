\name{glmbsim_cpp}
\alias{glmbsim_cpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates from the posterior density of a model in standard form
}
\description{
Simulates from the posterior density of a model in standard form
}
\usage{
glmbsim_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of draws to generate from posterior density. If length(n) > 1, the length is taken 
to be the number required.
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
  \item{alpha}{
an Offset vector}
  \item{wt}{
an vector of weights to be used in the fitting process. }
  \item{f2}{
fucntion used to calculate the negative of the log-posterior function
}
  \item{Envelope}{
an object of type \code{\link{glmbenvelope}}.
}
  \item{family}{
family used for simulation
}
  \item{link}{
link function used for simulation}
}
\details{
Simulates from the posterior density of a model in standard form
}
\value{
A list consisting of the following:
\item{out}{A matrix with simulated draws from a model in standard form. Each row represents on draw from the density}
\item{draws}{A vector with the number of candidates required before acceptance for each draw}
}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glmb}}, \code{\link{rglmb}}, \code{\link{glmbenvelope}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, y, x, mu, P, alpha, wt, f2, Envelope, family, link) 
{
    .Call("glmbayes2_glmbsim_cpp", PACKAGE = "glmbayes2", n, 
        y, x, mu, P, alpha, wt, f2, Envelope, family, link)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
