% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnvNeg_logLik.R,
%   R/rIndependent_Normal_Gamma_reg.R
\name{Neg_logLik}
\alias{Neg_logLik}
\alias{Neg_logLik2}
\title{Negative Log-Likelihood for a Generalized Linear Model}
\usage{
Neg_logLik(b, y, x, alpha, wt, family)

Neg_logLik2(b, y, x, alpha, wt, family)
}
\arguments{
\item{b}{a matrix where each row represents a point at which Negative Log-Likelihood is to be calculated}

\item{y}{A vector}

\item{x}{A matrix}

\item{alpha}{A vector}

\item{wt}{A vector}

\item{family}{A family}
}
\value{
The sum of \code{x} and \code{y}
}
\description{
Negative Log-Likelihood for a Generalized Linear Model
}
\examples{
1+1
10+1
}
\keyword{internal}
