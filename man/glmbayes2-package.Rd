\name{glmbayes2-package}
\alias{glmbayes2-package}
\alias{glmbayes2}
\docType{package}
\title{glmbayes2}
\description{
Produces iid sampling for Bayesian Generalized Linear models}
\details{
\tabular{ll}{
Package: \tab glmbayes\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2014-03-23\cr
License: \tab GPL-3\cr
}
The \code{glmbayes} package produces iid samples for Bayesian Genereralized Linear Models and is intended as a Bayesian version of the \code{\link{glm}} function for classical models. Estimation can be performed using three main functions. For models with fixed dispersion parameters, the \code{\link{rglmb}} function is the workhorse function and comes with a minimialistic interface for the input and output. It is also suitable for use as part of block Gibbs sampling procedures. The \code{\link{glmb}} function is essentially a wrapper function for the rglmb function that provides an interface closer to that of  \code{\link{glm}} function. The \code{\link{rglmbdisp}} function can be leveraged in order to produce samples for the dispersion parameters associated with the gaussian and Gamma link functions. Most methods defined for the output of the \code{\link{glm}} function are also defined for the \code{\link{rglmb}} and \code{\link{glmb}}, and \code{\link{rglmbdisp}} functions (see their respective documentation for details).

For the regression parameters, multivariate normal priors are assumed. Simulation for the gaussian family with the identify link function is performed using standard procedures for multivariate normal densities.  For all other families and link functions, simulation is performed using the likelihood subgradient approach of Nygren and Nygren (2006). This approach involves the construction of an enveloping function for the full posterior density followed by accept-reject sampling. For models that are approximately multivariate normal, the expected number of draws required per acceptance are bounded from above as noted in Nygren and Nygren (2006). 

Currently implemented models include the gaussian (identity link), poisson/quasipoisson (log link), binomial/quasibinomial (logit, probit, and cloglog links), and Gamma (log link) families. These models all have log-concave likelihood functions that allow us to leverage the likelihood-subgradient approach for the iid sampling. Models that fail to have log-concave likelihood functions are not implemented.  Our demos (viewable by entering the \code{demo()} command) provides examples of each of these families and links. 


The current implementation requires separate use of the \code{\link{rglmbdisp}} function in order to generate samples for dispersion parameters (gaussian, Gammma, quasipoisson, quasibinomaial families).  Our demos include examples of the joint use of the \code{\link{rglmb}} and \code{\link{rglmbdisp}} to produce samples for both regression and dispersion parameters using two-block Gibbs samplers. As these two-block Gibbs samplers likely are geometrically ergodic, future implementations may incorporate these two-block Gibbs samplers into the \code{\link{rglmb}} and \code{\link{glmb}} functions by leveraging theoretical bounds om convergence rates derived using Rosenthal (1996) type drift and minorization conditions.

The \code{\link{rglmb}} function can also be used in Block-Gibbs sampling implementations for Hierarchical Bayesian models. The demos associated with this package contains examples of such models.
}
\author{
Kjell Nygren

Maintainer: Kjell Nygren <kjell.a.nygren@gmail.com>
}
\references{
Nygren and Nygren (2006) Likelihood subgradient densities.
}
\keyword{ package }
\seealso{\code{\link{glm}}

}
\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(),x=TRUE)

n<-1000
mysd<-1

y<-glm.D93$y
x<-glm.D93$x
mu<-matrix(0,5)
#P<-matrix(0.1/((mysd)^2),nrow=5,ncol=5)+1/((mysd)^2)*diag(5)
P<-1/((mysd)^2)*diag(5)

glmtest<-glmb(counts ~ outcome + treatment, family = poisson()
,mu=mu,Sigma=solve(P),n=n,Gridtype=1)

glm.D93
glmtest

summary(glm.D93)
summary(glmtest)

fitted(glm.D93)
colMeans(fitted(glmtest))

residuals(glm.D93)
colMeans(residuals(glmtest))

extractAIC(glm.D93)
extractAIC(glmtest)

deviance(glm.D93)
mean(deviance(glmtest))

logLik(glm.D93)
mean(logLik(glmtest))

predict(glm.D93)
#predict(glmtest)  # Need to add function to produce predictions

colMeans(glmtest$linear.predictors)


}
