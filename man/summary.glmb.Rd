\name{summary.glmb}
\alias{summary.glmb}
\alias{print.summary.glmb}
\title{
Summarizes the output from the \code{\link{glmb}} function
}
\description{
The \code{summary.glmb} function summarizes the output from the \code{\link{glmb}} function. It takes an object of class \code{glmb} as an input. The output to a large extent mirrors the output from the \code{\link{summary.glm}} function. This is particularly true for the print output from the function (i.e. the output of the function \code{\link{print.summary.glmb}}). 
}
\usage{
\method{summary}{glmb}(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{rglmb} for which a summary is desired
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
The \code{summary.glmb} function summarizes the output from the \code{\link{glmb}} function. Key output includes mean residuals, information related to the prior, mean coefficients with associated stats, percentiles for the coefficients, as well as the effective number of parameters and the DIC statistic.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A list with the following elements:
 \item{n }{number of draws generated}
  \item{residuals}{vector of mean deviance residuals}
  \item{coefficients1}{Matrix with the prior mean and maximum likelihood coefficients with associated standard deviations}
  \item{coefficients}{Matrix with columns for the posterior mode, posterior mean, posterior standard deviation, monte carlo error, and tail probabilities (posterior probability of observing a value for the coefficient as extreme as the prior mean)}
  \item{Percentiles}{Matrix with estimated percentiles associated with the posterior density}
  \item{pD}{Estimated effective number of parameters}
  \item{deviance}{Vector with draws for the deviance}
  \item{DIC}{Estimated DIC statistic}
  \item{iters}{Average number of candidates per generated draws}
%% ...
}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
Spiegelhalter, et al,
}
\author{
Kjell Nygren
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glmb}}, \code{\link{rglmb}}, \code{\link{summary.rglmb}}}
\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

n<-1000
mysd<-1

mu<-matrix(0,5)
V0<-((mysd)^2)*diag(5)

glmb.D93<-glmb(n=n,counts ~ outcome + treatment, family = poisson(),mu=mu,Sigma=V0,Gridtype=1)

summary.glmb(glmb.D93)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
