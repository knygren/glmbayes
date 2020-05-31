\name{rmultireg}
\alias{rmultireg}
\concept{bayes}
\concept{multivariate regression}
\concept{simulation}


\title{Draw from the Posterior of a Multivariate Regression}

\description{
\code{ rmultireg} draws from the posterior of a Multivariate Regression model with a natural conjugate prior.
}

\usage{rmultireg(n,Y, X, Bbar, A, nu, V)}

\arguments{
  \item{n}{number of draws to generate from the density. }
  \item{Y    }{ \eqn{n x m} matrix of observations on m dep vars }
  \item{X    }{ \eqn{n x k} matrix of observations on indep vars (supply intercept) }
  \item{Bbar }{ \eqn{k x m} matrix of prior mean of regression coefficients }
  \item{A    }{ \eqn{k x k} Prior precision matrix }
  \item{nu   }{ d.f. parameter for Sigma }
  \item{V    }{ \eqn{m x m} pdf location parameter for prior on Sigma }
}

\details{
  Model: \cr
  \eqn{Y = XB + U} with \eqn{cov(u_i) = \Sigma} \cr
  \eqn{B} is \eqn{k x m} matrix of coefficients; \eqn{\Sigma} is \eqn{m x m} covariance matrix.

  Priors: \cr
  \eqn{\beta} | \eqn{\Sigma} \eqn{\sim}{~} \eqn{N(betabar, \Sigma(x) A^{-1})} \cr
  \eqn{betabar = vec(Bbar)};  \eqn{\beta = vec(B)} \cr
  \eqn{\Sigma} \eqn{\sim}{~} IW(nu, V) 
}

\value{
  A list of the components of a draw from the posterior
  \item{B }{ draw of regression coefficient matrix }
  \item{Sigma }{ draw of Sigma }
}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 2, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. \cr \url{http://www.perossi.org/home/bsm-1}}

\examples{


R=2000
set.seed(66)
n=200
m=1

X=matrix(cbind(rep(1,n),runif(n)),ncol=2) 
k=ncol(X)
B=matrix(c(1,2),ncol=m)


Sigma=matrix(c(1),ncol=m) 
RSigma=chol(Sigma)  

Y= matrix(X \%*\% B) + matrix(rnorm(m*n),ncol=m) \%*\% RSigma 


betabar=rep(0,k*m)
Bbar=matrix(betabar,ncol=m) # 2 x 1 
A=diag(rep(.01,k))
nu=3
V=nu*diag(m)

betadraw=matrix(double(R*k*m),ncol=k*m) 
Sigmadraw=matrix(double(R*m*m),ncol=m*m)

out=rmultireg(n,Y,X,Bbar,A,nu,V)

summary(out$B)
summary(out$Sigma)


}

\keyword{regression}