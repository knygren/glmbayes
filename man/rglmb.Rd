\name{rglmb}
\alias{rglmb}
%\alias{rglmb.default}
\alias{print.rglmb}
%\alias{rnorm_gamma_reg}
\alias{rnnorm_reg_cpp}
%\alias{rnorm_reg_cpp}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Generalized Linear Model Distributions
}
\description{
\code{rglmb} is used to generate iid samples for Bayesian generalized linear models, specified by providing a data vector, a design matrix, a description of the error distribution and a multivariate normal prior
}
\usage{
rglmb(n = 1, y, x, mu, P, wt = 1, dispersion=NULL,nu=NULL,
V=NULL,family = gaussian(), offset2 = rep(0, nobs), start = NULL, Gridtype = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of draws to generate from posterior density. If length(n) > 1, the length is taken 
to be the number required.
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
  \item{wt}{
an optional vector of weights to be used in the fitting process. Should be NULL or a numeric 
vector.
}
  \item{dispersion}{
the dispersion parameter for the family used. Either a single numerical value or NULL (the default). Must be provided for all
families other than the Poisson and Binomial.
}
\item{nu}{ Optional prior information for dispersion parameter (gaussian model only)}
\item{V}{ Optional prior information for dispersion parameter (gaussian model only)}

  \item{family}{
a description of the error distribution and link function to be used in the model. This can be a 
character string naming a family function, a family function or the result of a call to a family 
function. (See family for details of family functions.)

}
  \item{offset2}{
this can be used to specify an a priori known component to be included in the linear predictor during
 fitting. This should be NULL or a numeric vector of length equal to the number of cases.
}
  \item{start}{
an optional argument providing starting values for the posterior mode optimization.
}
  \item{Gridtype}{
an optional argument specifying the method used to determine number of likelihood subgradient 
densities used to construct the enveloping function.
}
}
\details{
The \code{rglmb} function produces iid samples for Bayesian generalized linear models. It has a more minimialistic interface than than the \code{\link{glmb}} function. Core required inputs for the function include the data vector, the design matrix and a prior specification. In addition, the dispersion parameter must currently be provided for the gaussian, Gamma, quasipoisson, and quasibinomial families (future implementations may incorporate a prior for these into the \code{rglmb} function).

For the gaussian family, iid samples from the posterior density is genererated using standard simulation procedures for multivariate normal densities. For all other families, the samples are generated using accept-reject procedures by leveraging the likelihood subgradient approach of Nygren and Nygren (2006). Current implemented models include the gaussian family (identity link function), the poisson and quasipoisson families (log link function), the gamma family (log link function), as well as the binomial and quasibinomial families (logit, probit, and cloglog link functions).  

}
\value{
\code{rglmb} returns a object of class "rglmb".  The function \code{summary} (i.e., \code{\link{summary.rglmb}}) can be used to obtain or print a 
summary of the results.  The generic accessor functions \code{\link{coefficients}}, \code{\link{fitted.values}}, \code{\link{residuals}}, and \code{\link{extractAIC}} can be used to extract various useful features of the value returned by \code{\link{rglmb}}.
An object of class "rglmb" is a list containing at least the following components:
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
Written by Kjell Nygren.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glm}}, \code{\link{glmb}}, \code{\link{rglmbdisp}}, \code{\link{summary.rglmb}}
}
\examples{

data(menarche)

summary(menarche)
plot(Menarche/Total ~ Age, data=menarche)

Age2=menarche$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche$Menarche/menarche$Total
wt=menarche$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 

out<-rglmb(n = 1000, y=y, x=x, mu=mu, P=solve(V1), wt = wt, 
family = binomial(logit), Gridtype = 3) 
summary(out)

mean(out$iters)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
