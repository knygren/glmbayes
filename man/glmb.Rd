\name{glmb}
\alias{glmb}
\alias{glmb.default}
\alias{print.glmb}
\alias{extractAIC.glmb}
\alias{logLik.glmb}
\alias{residuals.glmb}
\title{Producing iid samples for Bayesian Generalized Linear Models}
\description{
\code{glmb} is used to produce iid samples for Bayesian generalized linear models, specified by giving a symbolic description of the linear predictor, 
a description of the error distribution, and a multivariate normal prior.
}
\usage{
glmb(n, formula, family = binomial, dispersion=NULL, mu, Sigma, nu=NULL,V=NULL, Gridtype = 1, 
data, weights, subset2, na.action, start = NULL, etastart, mustart, offset, 
control = list(...), model = TRUE, method = "glm.fit", 
x = FALSE, y = TRUE, contrasts = NULL, ...)
}
\arguments{
  \item{n}{
number of draws to generate from posterior density. If length(n) > 1, the length is taken 
to be the number required.}
  \item{formula}{
an object of class formula (or one that can be coerced to that class): a symbolic 
description of the model to be fitted. The details of model specification are given under 
Details.}
  \item{family}{
a description of the error distribution and link function to be used in the model. This can be a 
character string naming a family function, a family function or the result of a call to a family 
function. (See family for details of family functions.)
}
  \item{dispersion}{
the dispersion parameter for the family used. Either a single numerical value or NULL (the default). Must be provided for all
families other than the Poisson and Binomial.
}
  \item{mu}{
a vector giving the prior means of the variables}
  \item{Sigma}{
a positive-definite symmetric matrix specifying the prior covariance matrix of the variables.}
\item{nu}{ Optional prior information for dispersion parameter (gaussian model only)}
\item{V}{ Optional prior information for dispersion parameter (gaussian model only)}
  \item{Gridtype}{
an optional argument specifying the method used to determine number of likelihood subgradient 
densities used to construct the enveloping function.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by as.data.frame to a data frame)
containing the variables in the model. If not found in data, the variables are taken from environment
(formula), typically the environment from which glm is called.}
  \item{weights}{
an optional vector of prior weights to be used in the fitting process. Should be NULL or a numeric 
vector.}
  \item{subset2}{
an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{na.action}{
a function which indicates what should happen when the data contain NAs. The default is set 
by the na.action setting of options, and is na.fail if that is unset. The factory-fresh 
default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.}
  \item{start}{
starting values for the parameters in the linear predictor.
}
  \item{etastart}{
starting values for the linear predictor.}
  \item{mustart}{
starting values for the vector of means.
}
  \item{offset}{
this can be used to specify an a priori known component to be included in the linear predictor during
 fitting. This should be NULL or a numeric vector of length equal to the number of cases. One or more 
offset terms can be included in the formula instead or as well, and if more than one is specified their 
sum is used. See model.offset.}
  \item{control}{
a list of parameters for controlling the fitting process.
}
  \item{model}{
a logical value indicating whether model frame should be included as a component of the returned value.
}
  \item{method}{
the method to be used in fitting the model. The default method glm.fit uses iteratively reweighted least 
squares (IWLS): the alternative model.frame returns the model frame and does no fitting.
User-supplied fitting functions can be supplied either as a function or a character string naming a function, 
with a function which takes the same arguments as glm.fit. If specified as a character string it is looked up from 
within the stats namespace. 
}
  \item{x,y}{
logical values indicating whether the response vector and model matrix used in the fitting process should be 
returned as components of the returned value. 
 }
  \item{contrasts}{
an optional list. See the contrasts.arg of model.matrix.default
}
  \item{\dots}{
arguments to be used to form the default control argument if it is not supplied directly.
}
}
\details{
The function \code{glmb} is a Bayesian version of the classical \code{\link{glm}} function.  Setup of the models mirrors those for the \code{\link{glm}} function but add additional required arguments \code{mu} and \code{Sigma}, representing a multivariate normal prior. In addition, the dispersion parameter must currently be provided for the gaussian, Gamma, quasipoisson, and quasibinomial families (future implementations may incoporate priors for these into the \code{glmb} function).  The function generates random iid samples from the posterior density associated with the model (assuming a fixed dispersion parameter). 

For the gaussian family, iid samples from the posterior density is genererated using standard simulation procedures for multivariate normal densities. For all other families, the samples are generated using accept-reject procedures by leveraging the likelihood subgradient approach of Nygren and Nygren (2006). Current implemented models include the gaussian family (identity link function), the poisson and quasipoisson families (log link function), the gamma family (log link function), as well as the binomial and quasibinomial families (logit, probit, and cloglog link functions).  

}
\value{
\code{glmb} returns an object of class "glmb". The function \code{summary} (i.e., \code{\link{summary.glmb}}) can be used to obtain or print a 
summary of the results.  The generic accessor functions \code{\link{coefficients}}, \code{\link{fitted.values}}, \code{\link{residuals}}, and \code{\link{extractAIC}} can be used to
extract various useful features of the value returned by code{\link{glmb}}.
An object of class "glmb" is a list containing at least the following components:

\item{glm}{a object of class glm containing the output from a call to the function glm}
\item{coefficients}{an n by length(mu) matrix with one sample in each row}
\item{coef.means}{a vector of length(mu) with the estimated posterior mean coefficients}
\item{coef.mode}{a vector of length(mu) with the estimated posterior mode coefficients}
\item{fitted.values}{an n by length(y) matrix with one sample for the mean values in each row}
\item{linear.predictors}{an n by length(y) matrix with one sample for the linear fit on the link scale in each row}
\item{deviance}{an n by 1 matrix with one sample for the deviance in each row}
\item{pD}{Estimate for effective number of parameters}
\item{Dbar}{Expected value for minus twice the log-likelihood function}
\item{Dthetabar}{Value of minus twice the log-likelihood function evaluated at mean value for coefficients}
\item{DIC}{Estimated Deviance Information criterion}
\item{famfunc}{Family functions used during estimation process}
\item{iters}{an n by 1 matrix giving the number of candidates generated before acceptance for each sample.}
}


\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Written by Kjell Nygren.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{glm}} for the classical version of these models, and \code{\link{rglmb}} for an alternate implementation with less overhead. See the function \code{\link{summary.glmb}} for details on the summary function accompanying this function.
}
\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

n<-1000
mysd<-1

mu<-matrix(0,5)
V0<-((mysd)^2)*diag(5)

glmb.D93<-glmb(n=n,counts ~ outcome + treatment, family = poisson(),mu=mu,Sigma=V0,Gridtype=1)

summary(glmb.D93)

# Menarche Binomial Data Example 

library(MASS)
data(menarche)

summary(menarche)
plot(Menarche/Total ~ Age, data=menarche)

n<-1000
mu<-matrix(0,nrow=2,ncol=1)
V1<-1*diag(2)
V1[1,1]<-100
V1[2,2]<-1

glm.out1<-glm(cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit), data=menarche)
glmb.out1<-glmb(n=n,cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit),mu=mu,Sigma=V1,Gridtype=1, data=menarche)

summary(glm.out1)
summary(glmb.out1)

glm.out2<-glm(cbind(Menarche,Total-Menarche)~Age,family=binomial(probit), data=menarche)
glmb.out2<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age,family=binomial(probit),mu=mu,Sigma=V1,Gridtype=1,data=menarche)

summary(glm.out2)
summary(glmb.out2)

glm.out3 <-glm(cbind(Menarche, Total-Menarche) ~ Age,family=binomial(cloglog), data=menarche)
glmb.out3<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age,family=binomial(cloglog),mu=mu,Sigma=V1,Gridtype=1,data=menarche)

summary(glm.out3)
summary(glmb.out3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
