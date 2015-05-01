\name{glmbdic}
\alias{glmbdic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates DIC and Deviance Information
}
\description{
Calculates DIC and Deviance Information

}
\usage{
glmbdic(coefficients, y, x, alpha=0,f1, f4, wt = 1,dispersion=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coefficients}{
A matrix with coefficients from the rglmb function
}
  \item{y}{
a vector of observations of length m
}
  \item{x}{
a design matrix of dimension m*p
}
  \item{alpha}{
an offset parameter or vector
}
  \item{f1}{
A function used to calculate the negative log-likelihood
}
  \item{f4}{
A function used to calculate the deviance
}
  \item{wt}{
a vector of weights
}
  \item{dispersion}{
dispersion parameter
}

}
\details{
Calculates DIC and Deviance Information
}
\value{
A list with the following components
\item{Deviance}{A n*1 matrix with the deviance for each draw}
\item{Dbar}{Mean for negative 2 times negative log-likelihood}
\item{Dthetabar}{Negative 2 times log-likelihood evaluated at mean parameters}
\item{pD}{Effective number of parameters}
\item{DIC}{DIC statistic}

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
