\name{rnnorm_reg_std_cpp}
\alias{rnnorm_reg_std_cpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Bayesian Generalized Linear Model Distribution in Standard Form
}
\description{
\code{rnnorm_reg_std_cpp} is used to generate iid samplers from Non-Gaussian 
Generalized Linear Models in standard form. The function should only
be called after standardization of a Generalized Linear Model.
}
\usage{
rnnorm_reg_std_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link,progbar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
}
  \item{y}{
a vector of observations of length \code{m}.
}
  \item{x}{
a design matrix of dimension \code{m * p}.
}
  \item{mu}{
a vector of length \code{p} giving the prior means of the variables in the design matrix.
}
  \item{P}{
a positive-definite symmetric matrix of dimension \code{p * p} specifying the prior precision matrix of the variable.}
  \item{alpha}{
this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.}
  \item{wt}{
an optional vector of \sQuote{prior weights} to be used in the fitting process. Should be NULL or a numeric vector. }
  \item{f2}{
function used to calculate the negative of the log-posterior function
}
  \item{Envelope}{
an object of type \code{glmbenvelope}.
}
  \item{family}{
family used for simulation. Used that this is different from the family used in other functions.
}
  \item{link}{
link function used for simulation.
}
  \item{progbar}{dummy for flagging if a progressbar should be produced during the call}

}
\details{
Simulates from the posterior density of a model in standard form
}
\value{
A list consisting of the following:
\item{out}{A matrix with simulated draws from a model in standard form. Each row represents one draw from the density}
\item{draws}{A vector with the number of candidates required before acceptance for each draw}
}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glmb}}, \code{\link{rglmb}}, \code{\link{EnvelopeBuild_c}}
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
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 1 times as large as point estimate 

out<-rglmb(n = 1000, y=y, x=x, mu=mu, P=solve(V1), wt = wt, 
           family = binomial(logit), Gridtype = 3) 
summary(out)

famfunc<-glmbfamfunc(binomial(logit))

f1<-famfunc$f1
f2<-famfunc$f2  # Used in optim and glmbsim_cpp
f3<-famfunc$f3  # Used in optim
f5<-famfunc$f5
f6<-famfunc$f6

dispersion2<-as.numeric(1.0)
start <- mu
offset2=rep(as.numeric(0.0),length(y))
P=solve(V1)
n=1000

## Appears that the type for some of these arguments are important/problematic

outlist<-rnnorm_reg_cpp(n=as.integer(n),y=as.vector(y),
x=as.matrix(x),mu=as.vector(mu),P=as.matrix(P),
offset2=as.vector(offset2),wt=as.vector(wt),dispersion=as.numeric(dispersion2),
famfunc=famfunc,f1=f1,f2=f2,f3=f3,start=as.vector(start),family="binomial",
link="logit",Gridtype=as.integer(3))

### This allows use of the rglmb summary function 
### add interface for glmbsim_NGauss_cpp later

#outlist$call<-match.call()
colnames(outlist$coefficients)<-colnames(x)
class(outlist)<-c(outlist$class,"rglmb")
#summary(outlist)
Env1=outlist$Envelope


###### Adjust weight for dispersion

wt2=wt/dispersion2

######################### Shift mean vector to offset so that adjusted model has 0 mean

alpha=x\%*\%as.vector(mu)+offset2
mu2=0*as.vector(mu)
P2=P
x2=x


#####  Optimization step to find posterior mode and associated Precision

parin=start-mu

opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x),mu=as.vector(mu2),
              P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
              method="BFGS",hessian=TRUE
)

bstar=opt_out$par  ## Posterior mode for adjusted model
bstar
bstar+as.vector(mu)  # mode for actual model
A1=opt_out$hessian # Approximate Precision at mode

## Standardize Model

Standard_Mod=glmb_Standardize_Model(y=as.vector(y), x=as.matrix(x),P=as.matrix(P),
                                    bstar=as.matrix(bstar,ncol=1), A1=as.matrix(A1))

bstar2=Standard_Mod$bstar2  
A=Standard_Mod$A
x2=Standard_Mod$x2
mu2=Standard_Mod$mu2
P2=Standard_Mod$P2
L2Inv=Standard_Mod$L2Inv
L3Inv=Standard_Mod$L3Inv

Env2=EnvelopeBuild_c(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),as.vector(wt2),
family="binomial",link="logit",Gridtype=as.integer(3), 
n=as.integer(n),sortgrid=TRUE)

## These now seem to match

Env1
Env2

#int n, NumericVector y, NumericMatrix x, NumericMatrix mu, NumericMatrix P, 
#NumericVector alpha, NumericVector wt, Function f2,
#Rcpp::List Envelope, Rcpp::CharacterVector family, 
#Rcpp::CharacterVector link, int progbar

### Note: getting the types correct here is important but potentially difficult for users
### May be better to call an R function wrapper that checks (and converts when possible) 
##  to correct types

sim=rnnorm_reg_std_cpp(n=as.integer(n),y=as.vector(y),x=as.matrix(x2),mu=as.matrix(mu2,ncol=1),
P=as.matrix(P2),alpha=as.vector(alpha),wt=as.vector(wt2),
f2=f2,Envelope=Env2,family="binomial",link="logit")

out=L2Inv\%*\%L3Inv\%*\%t(sim$out)

for(i in 1:n){
  out[,i]=out[,i]+mu
}

summary(t(out))
mean(sim$draws)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
