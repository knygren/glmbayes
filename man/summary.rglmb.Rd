\name{summary.rglmb}
\alias{summary.rglmb}
\alias{summary.rglmbdisp}
\alias{print.summary.rglmb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes the output from the \code{\link{rglmb}} function
}
\description{
The \code{summary.rglmb} function summarizes the output from the \code{\link{rglmb}} function. It takes an object of class \code{rglmb} as an input. The output to a large extent mirrors the output from the \code{\link{summary.glm}} function. This is particularly true for the print output from the function (i.e. the output of the function \code{\link{print.summary.rglmb}}).
}
\usage{
\method{summary}{rglmb}(object, ...)
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
The \code{summary.rglmb} function summarizes the output from the \code{\link{rglmb}} function. Key output includes mean residuals, information related to the prior, mean coefficients with associated stats, percentiles for the coefficients, as well as the effective number of parameters and the DIC statistic.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{n}{number of draws generated}
  \item{coefficients1}{Matrix with the prior mean and approximate weight for the prior relative to the data}
  \item{coefficients}{Matrix with columns for the posterior mode, posterior mean, posterior standard deviation, monte carlo error, and tail probabilities (posterior probability of observing a value for the coefficient as extreme as the prior mean)}
  \item{Percentiles}{Matrix with estimated percentiles associated with the posterior density}
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
\code{\link{rglmb}}, \code{\link{glmb}}, \code{\link{summary.glmb}}
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

out<-rglmb(n = 1000, y=y, x=x, mu=mu, P=solve(V1), wt = wt, 
family = binomial(logit), Gridtype = 3) 
summary(out)

mean(out$iters)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
