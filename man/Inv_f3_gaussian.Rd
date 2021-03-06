% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Inv_f3_gaussian.R
\name{Inv_f3_gaussian}
\alias{Inv_f3_gaussian}
\title{Derives the inverse of the gradient vector.}
\usage{
Inv_f3_gaussian(cbars, y, x, mu, P, alpha, wt)
}
\arguments{
\item{cbars}{Gradient vectors desired for new thetabars. Typically the output of an initial call to 
Envelopebuild.}

\item{y}{For \code{glm}:
    logical values indicating whether the response vector and model
    matrix used in the fitting process should be returned as components
    of the returned value.

    For \code{glm.fit}: \code{x} is a design matrix of dimension
    \code{n * p}, and \code{y} is a vector of observations of length
    \code{n}.
  }

\item{x}{For \code{glm}:
    logical values indicating whether the response vector and model
    matrix used in the fitting process should be returned as components
    of the returned value.

    For \code{glm.fit}: \code{x} is a design matrix of dimension
    \code{n * p}, and \code{y} is a vector of observations of length
    \code{n}.
  }

\item{mu}{Prior mean}

\item{P}{Prior Precision matrix}

\item{alpha}{offset vector}

\item{wt}{weighting vector}
}
\value{
Refer to Nygren and Nygren (2006) for details. The first set of items refers to Example 2 in section 
3.1. All except the last item in this list of returned items has a number of rows equaling the number of 
components of the grid and a number of columns equaling the number of coefficients in the model. All quantities 
refer to the respective coefficient for each of the components of the grid.
\item{Down}{The lower bounds for the interval to be evaluated. Either negative infinity or a real number.}
\item{Up}{The upper bounds for the interval to be evaluated. Either positive infinity or a real number.}
\item{lglt}{The log of the density between negative infinity and the upper bound}
\item{lgrt}{The log of the density between the lower bound and infinity}
\item{lgct}{The log of the density between the lower and upper bounds}
\item{logU}{The one of the 3 above that is relevant for the component of the grid}
\item{logP}{A two column matrix, the first of which holds sum of logU across the components. The second column is 
0 and is later populated by the Set_logP function}
}
\description{
Computes the inverse of the gradient vector for the gaussian model. Typically done to find the set 
of tangency points that yield the same gradient as an initial set of gradients used for an envelope 
positioned at the posterior mode for a specific dispersion.
}
\examples{
## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb extractAIC-------------------------------------------------------------
extractAIC(glmb.D93)
}
