% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confint.glmb.R
\name{confint.glmb}
\alias{confint.glmb}
\title{Credible Intervals for Model Parameters}
\usage{
\method{confint}{glmb}(object, parm, level = 0.95, ...)
}
\arguments{
\item{object}{a fitted model object of class \code{"glmb"}. Typically the result of a call to \link{glmb}.}

\item{parm}{a specification (not yet implemented) of which parameters are to be given credible sets,
either a vector of numbers or a vector of names. If missing, all parameters are considered.}

\item{level}{the credible interval required.}

\item{\ldots}{additional argument(s) for methods.}
}
\value{
A matrix (or vector) with columns giving lower and
upper credible limits for each parameter. These will be labeled
(1-level)/2 and 1-(1-level)/2 in \% (by default 2.5\% and 97.5\%).
}
\description{
Computes credible intervals for Model Parameters
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
## ----glmb confint-------------------------------------------------------------
confint(glmb.D93)
}
