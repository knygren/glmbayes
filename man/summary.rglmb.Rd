\name{summary.rglmb}
\alias{summary.rglmb}
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

data(carinsca)
carinsca$Merit <- ordered(carinsca$Merit)
carinsca$Class <- factor(carinsca$Class)
options(contrasts=c("contr.treatment","contr.treatment"))
attach(carinsca)


out <- glm(Cost/Claims~Merit+Class,family=Gamma(link="log"),weights=Claims,x=TRUE)
summary(out)

V_Out=summary(out)$cov.scaled  
P_Data=solve(V_Out)

Like_std=summary(out)$coefficients[,2]

y1<-out$y
x1<-out$x
b1<-out$coefficients
wt1<-out$prior.weights
scale<-0.1275 # SAS estimate

dispersion<-1/scale
wt2<-wt1/dispersion
alpha1<-rep(0,length(y1))

mu<-matrix(0,8)
P<-0.1*solve(diag(Like_std*Like_std))
mu[1,1]=-1.1

out2<-rglmb(n=1000,out$y,out$x,mu=mu,P=P,wt=out$prior.weights,dispersion=
dispersion,family=Gamma(link="log"),offset2=rep(0,20),start=out$coefficients,
Gridtype = 3) 

summary.rglmb(out2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
