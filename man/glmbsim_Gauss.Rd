\name{glmbsim_Gauss_cpp}
\alias{glmbsim_Gauss_cpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates From Gaussian Models with dispersion
}
\description{
Simulates From Gaussian Models with dispersion
}
\usage{
glmbsim_Gauss_cpp(n,y,x,mu,P,offset2,wt,dispersion,famfunc,f1,f2,f3,start,
family,link, Gridtype)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  number of draws to generate from posterior density.
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
  \item{offset2}{
this can be used to specify an a priori known component to be included in the linear predictor during
 fitting. This should be NULL or a numeric vector of length equal to the number of cases.
}
  \item{wt}{
vector of weights to be used in the fitting process. 
}
  \item{dispersion}{
dispersion parameter
}

  \item{famfunc}{
family function
}
  \item{f1}{
f1 function (negative log-likelihood)
}
  \item{f2}{
f2 function (negative log-posterior)
}
  \item{f3}{
f3 function (gradient vector)
}
  \item{start}{
an optional argument providing starting values for the posterior mode optimization.}

\item{family}{
family
}

\item{link}{
link function
}

\item{Gridtype}{
an optional argument specifying the method used to determine number of likelihood subgradient 
densities used to construct the enveloping function.
}


}
\details{
Simulates from the Gaussian regression function
}
\value{
A list with the following components
\item{coefficients}{an n by length(mu) matrix with one sample in each row}
\item{PostMode}{a vector of length(mu) with the estimated posterior mode coefficients}
\item{Prior}{A list with two components. The first being the prior mean vector and the second the prior precision matrix}
\item{iters}{an n by 1 matrix giving the number of candidates generated before acceptance for each sample.}
\item{famfunc}{an object of class famfunc}
\item{Envelope}{an object of class envelope}
\item{dispersion}{the dispersion parameter used in the model}
\item{loglike}{an n by 1 matrix containing the negative loglikelihood for each sample.}


}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glmb}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (coefficients, y, x, f1, f4, wt = 1) 
{
    l1 <- length(coefficients[1, ])
    l2 <- length(coefficients[, 1])
    D <- matrix(0, nrow = l2, ncol = 1)
    D2 <- matrix(0, nrow = l2, ncol = 1)
    for (i in 1:l2) {
        b <- as.vector(coefficients[i, ])
        D[i, 1] <- f4(b = b, y = y, x = x, alpha = 0, wt = wt)
        D2[i, 1] <- 2 * f1(b = b, y = y, x = x, alpha = 0, wt = wt)
    }
    Dbar <- mean(D2)
    b <- colMeans(coefficients)
    Dthetabar <- 2 * f1(b = b, y = y, x = x, alpha = 0, wt = wt)
    pD <- Dbar - Dthetabar
    DIC <- pD + Dbar
    list(Deviance = D, Dbar = Dbar, Dthetabar = Dthetabar, pD = pD, 
        DIC = DIC)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
