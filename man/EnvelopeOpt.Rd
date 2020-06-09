\name{EnvelopeOpt}
\alias{EnvelopeOpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimizes Envelope function for simulation
}
\description{
Optimizes the grid for the envelope
}
\usage{
EnvelopeOpt(a1, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a1}{Diagonal elements of data precision matrix for a model in standard form
}
  \item{n}{Number of draws to generate}
}
\details{
This function attempts to find a computationally optimal gridsize by using information on the strenght of the prior and the number of iterations desired. Generally, more data (i.e., larger values for the diagonal elements of the precision matrix) will require a larger grid. The same also holds when the number of desired draws is higher (as the setup costs associated with the larger grid is offset by the savings in the number of candidates per sample).
}
\value{
A vector containing information on how many component each dimension should be split into.
}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}

\seealso{
\code{\link{rglmb}}, \code{\link{EnvelopeBuild_c}}, \code{\link{expand.grid}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (a1, n) 
{
    a1rank <- rank(1/(1 + a1))
    l1 <- length(a1)
    dimcount <- matrix(0, (l1 + 1), l1)
    scaleest <- matrix(0, (l1 + 1), l1)
    intest <- c(1:(l1 + 1))
    slopeest <- c(1:(l1 + 1))
    dimcount[1, ] <- diag(diag(l1))
    scaleest[1, ] <- 1 + a1
    slopeest[1] <- prod(scaleest[1, ])
    for (i in 2:(l1 + 1)) {
        dimcount[i, ] <- dimcount[i - 1, ]
        scaleest[i, ] <- scaleest[i - 1, ]
        for (j in 1:l1) {
            if (a1rank[j] == i - 1) {
                dimcount[i, j] <- 3
                scaleest[i, j] <- 2/sqrt(pi)
            }
        }
        intest[i] <- 3^(i - 1)
        slopeest[i] <- prod(scaleest[i, ])
    }
    evalest <- intest + n * slopeest
    minindex <- 0
    for (j in 1:(l1 + 1)) {
        if (evalest[j] == min(evalest)) {
            minindex <- j
        }
    }
    dimcount[minindex, ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
