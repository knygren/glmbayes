\name{Prior_Likelihood_Check}
\alias{Prior_Likelihood_Check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparison of Prior Specification and Maximum Likelihood Estimates
}
\description{
Calculates the difference between the prior and the maximum likelihood estimates and standardizes to understand how many standard deviations apart the two are (using both prior and maximum likelihood standard deviations). Issues warning messages if the two are more than 3 standard deviations apart based on either of the two measures.
}
\usage{
Prior_Likelihood_Check(Prior_mean, Prior_std, Like_est, Like_std)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Prior_mean}{
a vector with the intended prior mean
}
  \item{Prior_std}{
a vector with the intended prior standard deviation
}
  \item{Like_est}{
a vector of maximum likelihood estimates typically generated from a call to glm
}
  \item{Like_std}{
an vector of standard deviations associated with the maximum likelihood estimate that typical is generated from a call to the summary function for glm model output
}

}
\details{
This function returns a list of four vectors containing information on the differences between the prior and the maximum likelihood estimates. The differences are provided both as signed differences as well as in absolute value terms. The first column is defined as (Like_est-Prior_mean)/Prior_std, while the second column is defined as (Prior_mean-Like_est)/Like_std. The 3rd and 4th columns provides the same values in absolute value terms. 
}
\value{
A matrix defined as follows
\item{Error_Checks}{A n*4 matrix where the first two columns are signed differences and the last two columns are provided as absolute values.  
}

}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glmb}}
}
\examples{

## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))


glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(),x=TRUE)

summary(glm.D93)

n<-1000
mu<-matrix(0,5)
X<-glm.D93$x
Xmu=glm.D93$x\%*\%glm.D93$coefficients
explambda=diag(as.vector(exp(Xmu)))

# Use approximately conjugate prior for mu with 10% weight on prior

wt_0<-0.1
m_0=exp(log(wt_0/(1-wt_0)))
V0<-solve(m_0*t(X)\%*\%explambda\%*\%X)

diag(sqrt(diag(V0)))
V0
sqrt(diag(V0))


Like_std=summary(glm.D93)$coefficients[,2]
D93_Prior_Error_Checks=Prior_Likelihood_Check(mu,
sqrt(diag(V0)),glm.D93$coefficients,Like_std)

D93_Prior_Error_Checks


mu[1,1]=0+4*sqrt(diag(V0)[1])
mu[1,1]=0+15*Like_std[1]
mu[1,1]=log(mean(counts))
mu

D93_Prior_Error_Checks=Prior_Likelihood_Check(mu,
sqrt(diag(V0)),glm.D93$coefficients,Like_std)

#glmb.D93<-glmb(n=n,counts ~ outcome + treatment, family = poisson(),mu=mu,Sigma=V0,Gridtype=3)
#summary(glmb.D93)

  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
