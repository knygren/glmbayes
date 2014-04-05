\name{Set_Grid}
\alias{Set_Grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setup key constants in the grid used during simulation
}
\description{
Setup key constants in the grid used during simulation
}
\usage{
Set_Grid(GIndex, cbars, Lint)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GIndex}{A matrix containing information on how each dimension should be sampled (1 means left tail of a restricted normal, 2 center, 3 right tail, and 4 the entire line)
  }
  \item{cbars}{
A matrix containing the gradients for the negative log-likelihood at each tangency}
  \item{Lint}{A matrix containing information on the lower and upper bound for each dimension
}
}
\details{
This function sets constants used during the simulation.
}
\value{

\item{Down}{A matrix containing a vector or lower bounds for each component of the grid}
\item{Up}{A matrix containing vector of upper bounds for each components of the grid}
\item{lglt}{A matrix containing a vector of log left-tails for each component of the grid}
\item{lgrt}{A matrix containing a vector of log right-tails for each component of the grid}
\item{lgct}{A matrix containing a vector of log centers for each component of the grid}
\item{lgU}{A matrix containing a vector of log-probabilies for each component of the grid}
\item{lgP}{A matrix containing log probabilities for the log-probabilities for each component of the grid}
}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}

\seealso{
\code{\link{rglmb}}, \code{\link{expand.grid}}, \code{\link{optgrid}},  \code{\link{setlogP}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (GIndex, cbars, Lint) 
{
    .Call("glmbayes2_Set_Grid", PACKAGE = "glmbayes2", GIndex, 
        cbars, Lint)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
