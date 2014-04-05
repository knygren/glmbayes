\name{setlogP}
\alias{setlogP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set Probabilities of each part of the grid from which enveloping function is built and also sets a constant used during simulation.
}
\description{
Set Probabilities of each part of the grid from which enveloping function is built and also sets a constant used during simulation.}
\usage{
setlogP(logP, NegLL, cbars, G3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logP}{
A matrix with log probabilities associated with the components of the grid}
  \item{NegLL}{
A vector where each element represents the Negative log-likelihood associated with each component of the grid.
}
  \item{cbars}{
A matrix where each row consist of the gradient for the negative log-likelihood associated with the respective component of the grid
}
  \item{G3}{
A matrix containing the points of tangency for each component of the grid.
}
}
\details{
Set Probabilities of each part of the grid from which enveloping function is built and also sets a constant used during simulation.
}
\value{
A list with two components
\item{logP}{A matrix containing log-probabilities related to the components of the grid}
\item{LLconst}{A vector containing constant for each component of the grid used during the accept-reject procedure}
}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rglmb}}, \code{\link{expand.grid}}, \code{\link{optgrid}},  \code{\link{Set_Grid}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (logP, NegLL, cbars, G3) 
{
    .Call("glmbayes2_setlogP", PACKAGE = "glmbayes2", logP, NegLL, 
        cbars, G3)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
