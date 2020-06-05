\name{rglmbdisp}
\alias{rglmbdisp}
%\alias{rglmbdisp.default}
\alias{print.rglmbdisp}

\title{
Produces iid samples for the conditional density associated with the posterior dispersion.
}
\description{
\code{rglmbdisp} is used to produce iid samples for the dispersion parameter associated with Bayesian generalized linear models, specified by providing a data vector, a design matrix, a vector for the regression parameters, a description of the error distribution, and an inverse gamma prior.
}
\usage{
rglmbdisp(n, y, x, b, alpha = 0,wt=1, shape, rate,family=gaussian())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of draws to generate from posterior density. If length(n) > 1, the length is taken to be the number required.

}
  \item{y}{
a vector of observations of length m
}
  \item{x}{
a design matrix of dimension m*p
}
  \item{b}{
a vector of regression parameters of length p}
  \item{alpha}{
an optional vector of length m containing offsets
}
  \item{wt}{
an optional vector of weights to be used in the fitting process. Should be NULL or a numeric 
vector.
}
  \item{shape}{
shape parameter for the prior inverse gamma density
}
  \item{rate}{
rate parameter for the prior inverse gamma density
}
  \item{family}{
a description of the error distribution and link function to be used in the model. This can be a 
character string naming a family function, a family function or the result of a call to a family 
function. (See family for details of family functions.)
}
}
\details{
For the gaussian and gamma regression models, this function assumes an inverse gamma prior distribution for the dispersion parameter. In the case of the gaussian regression model, this is a conjugate prior for the conditional density and sampling is performed using standard procedures for the inverse gamme distribution. For the gamma regression model, a likelihood subgradient based approach (see Nygren and Nygren 2006) is used to generate samples from the required density. Future implementations may add an ability to simulate from the dispersion parameters associated with the quasipoisson and quasibinomial models (although these models don't have proper posterior densites).  
}
\value{
A n*1 matrix containing draws from the posterior density
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glmb}}, \code{\link{rglmb}}
}
\examples{

data(carinsca)
carinsca$Merit <- ordered(carinsca$Merit)
carinsca$Class <- factor(carinsca$Class)
options(contrasts=c("contr.treatment","contr.treatment"))
attach(carinsca)
out <- glm(Claims/Insured~Merit+Class,family="poisson")
summary(out,cor=FALSE)

wt1<-out$prior.weights
scale<-0.1275 # SAS estimate

dispersion<-1/scale

out <- glm(Cost/Claims~Merit+Class,family=Gamma(link="log"),weights=Claims,x=TRUE)
summary(out)


###  Dispersion Model  ######################

a0<-0.01
m0<-0.01
b<-out$coefficients

out3<-rglmbdisp(n=10000,y=out$y,x=out$x,b=b,alpha= rep(0, length(out$y))
,wt=out$prior.weights,shape=m0,rate=a0,family=Gamma(link=log))
summary(out3)

### Development Code may have attempts at calculating the DIC (not sure about rglmbdisp)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
