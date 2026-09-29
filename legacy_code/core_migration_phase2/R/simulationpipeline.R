#' Low-Level Simulation Pipeline for Bayesian GLMs
#'
#' @description
#' A detailed overview of the low-level simulation pipeline used by
#' \code{rglmb()} and related functions. These routines implement the
#' optimization -> standardization -> envelope sizing -> envelope construction ->
#' sampling -> back-transformation workflow described in \insertCite{Nygren2006}{glmbayes}.
#'
#' @details
#' (summaries of each step)
#'
#' @references
#' \insertAllCited{}
#'
#' @example inst/examples/Ex_rNormalGLM_std.R
#'
#' @name SimulationPipeline
NULL




#' Return family functions used during simulation and post processing
#'
#' This function takes as input a \code{\link{family}} object and returns a 
#' set of functions that are used during simulation and summarization of models 
#' using the \code{\link{glmb}}, and \code{\link{rglmb}} functions.
#' 
#' @name glmbfamfunc
#' @aliases
#' glmbfamfunc
#' print.glmbfamfunc
#' @param family an object of class \code{\link{family}}
#' @param lik_shape Known shape parameter of the Gamma likelihood; used only for the
#'   \code{Gamma(link = "identity")} branch where the regression coefficient is the
#'   Gamma \emph{rate}.  Ignored for all other family/link combinations.
#'   Defaults to \code{1} (i.e.\ exponential likelihood).
#' @param x an object of class \code{"glmbfamfunc"} for which a printed output is desired.
#' @param \ldots additional optional arguments
#' @return A list (class \code{"glmbfamfunc"}) whose first four components are **always**
#'   present for every supported \code{family} and \code{link}. The names \code{f1}--\code{f4}
#'   are stable: they mean the same roles across families (only the internal formulas change).
#'   \describe{
#'   \item{\code{f1}}{Negative log-likelihood as a function of coefficients \code{b}
#'     (arguments typically \code{b}, \code{y}, \code{x}, optional \code{alpha}, \code{wt}).}
#'   \item{\code{f2}}{Negative log-posterior (likelihood plus Normal prior quadratic form in
#'     \code{b} with precision \code{P} and mean \code{mu}).}
#'   \item{\code{f3}}{Gradient of \code{f2} with respect to \code{b} (same argument pattern as \code{f2}).}
#'   \item{\code{f4}}{Deviance-related quantity (twice negative log-likelihood contrast vs.\
#'     saturated model, with a \code{dispersion} argument for quasi-families); used in DIC-style summaries.}
#'   \item{\code{f7}}{Family-specific matrix: weighted sum of outer products of predictor rows,
#'     i.e.\ a curvature / expected negative Hessian of the log-likelihood w.r.t.\ \code{b}
#'     at the supplied \code{b} (used e.g.\ by \code{\link{directional_tail}} for \code{glmb} fits).}
#'   }
#'   Slots \code{f5} and \code{f6} are **not** returned: they were reserved for alternate or
#'   C++-aligned likelihood/posterior routines and remain commented out in the implementation
#'   (only \code{f1}, \code{f2}, \code{f3}, \code{f4}, and \code{f7} are assigned in the returned list).
#' @details
#'   \code{glmbfamfunc} is the canonical R closure bundle for likelihood, posterior, gradient, and
#'   deviance quantities across all supported family/link combinations.
#'
#'   **Registration requirement.**  A branch must exist inside \code{glmbfamfunc} for every
#'   family/link combination that is used in the package.  If a combination is missing, the
#'   closures \code{f1}--\code{f4} will never be assigned and any downstream code that accesses
#'   them will fail with \emph{"object 'f1' not found"}.  New family/link combinations must
#'   therefore be explicitly added here before they can produce valid DIC, log-likelihood, or
#'   directional-tail output.
#'
#'   **Currently implemented family/link combinations:**
#'
#'   \tabular{ll}{
#'     \strong{Family}                    \tab \strong{Link} \cr
#'     \code{gaussian}                    \tab \code{identity} \cr
#'     \code{poisson}, \code{quasipoisson} \tab \code{log} \cr
#'     \code{poisson}, \code{quasipoisson} \tab \code{identity} \cr
#'     \code{binomial}, \code{quasibinomial} \tab \code{logit} \cr
#'     \code{binomial}, \code{quasibinomial} \tab \code{probit} \cr
#'     \code{binomial}, \code{quasibinomial} \tab \code{cloglog} \cr
#'     \code{binomial}, \code{quasibinomial} \tab \code{identity} \cr
#'     \code{Gamma}                       \tab \code{log} \cr
#'     \code{Gamma}                       \tab \code{identity} \cr
#'   }
#'
#'   Any family/link not in this table will fall through all branches silently and produce
#'   the error above at the point of first use.  For \code{Gamma(link = "identity")}, pass the
#'   known Gamma likelihood shape \code{lik_shape} (default \code{1}) so that \code{f1}--\code{f4}
#'   and \code{f7} use the correct parameterization (coefficient = Gamma rate \eqn{\beta}).
#'
#'   **Relationship to C++ simulation paths.**  Many simulation procedures in the package have
#'   been fully or partially migrated to \code{*.cpp} routines, which receive their own
#'   objective functions directly.  For those paths \code{glmbfamfunc} may not be called at all
#'   during sampling.  However, the R-side post-processing functions
#'   (\code{\link{logLik.glmb}}, \code{\link{summary.rglmb}}, \code{\link{directional_tail}})
#'   always use \code{famfunc\$f1}, \code{famfunc\$f4}, and \code{famfunc\$f7} respectively, so a
#'   registered branch is still required for those outputs even when the sampler itself has moved
#'   to C++.
#' @example inst/examples/Ex_glmbfamfunc.R
#' @export
#' @rdname glmbfamfunc
#' @order 1

glmbfamfunc<-function(family, lik_shape = 1){
  
  # need to add handling for offsets 
  
  # f1-(negative)Log-Likelihood function
  # f2-(negative)Log-Posterior density
  # f3-(negative) Gradient for Log-Posterior density 
  # f4-Deviance function (note multiplies the difference in two times log-likelihood by dispersion)
  
  if(family$family=="gaussian" && family$link=="identity")
  {
    f1<-function(b,y,x,alpha=0,wt=1){
      Xb<-alpha+x%*%b
      -sum(dnorm(y, mean=Xb,sd=sqrt(1/wt),log=TRUE))
    }
    f2<-function(b,y,x,mu,P,alpha=0,wt=1){
      Xb<-alpha+x%*%b
      -sum(dnorm(y, mean=Xb,sd=sqrt(1/wt),log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
    }
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      l2<-length(y)
      ltemp<-length(wt)
      yxb2<-NULL
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      xb<-alpha+x%*%b-y
      yXb2<-P%*%(b-mu)+t(x)%*%(Ptemp%*%xb)
      yXb2
    }
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      
      ## adds 2*log-likelihood at mean=y which measn this is 2* Negative log-likelihood difference from saturated
      ## Multiplication by dispersion likely is/was incorrect
      ## As long as single value for the parameters and dispersion, can do the below
      
      #dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dnorm(y, mean=y,sd=sqrt(1/wt),log=TRUE)))
      2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dnorm(y, mean=y,sd=sqrt(dispersion/wt),log=TRUE))
      
    }
    
    # Edit This
    #	f5<-f2_gaussian
    #	f6<-f3_gaussian
    
    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
      l2<-length(y)
      ltemp<-length(wt)
      yxb2<-NULL
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      Pout<-t(x)%*%(Ptemp)%*%x
      Pout
    }
    
    
    
    
    
    
  }
  
  # Check if Poisson weight should be outside or inside function
  
  
  if((family$family=="poisson"||family$family=="quasipoisson") && family$link=="log")
  {
    f1<-function(b,y,x,alpha=0,wt=1){
      lambda<-t(exp(alpha+x%*%b))
      #-sum(dpois(y, lambda,log=TRUE)*wt)
      -sum(dpois2(y, lambda,log=TRUE)*wt)
      
      
      
    }
    f2<-function(b,y,x,mu,P,alpha=0,wt=1){
      lambda<-t(exp(alpha+x%*%b))
      #-sum(dpois(y, lambda,log=TRUE)*wt)+0.5*t((b-mu))%*%P%*%(b-mu)
      -sum(dpois2(y, lambda,log=TRUE)*wt)+0.5*t((b-mu))%*%P%*%(b-mu)
    }
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      -t(x)%*%((y-exp(alpha+x%*%b))*wt)+P%*%(b-mu)
      
    }
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      #dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dpois(y, y,log=TRUE)*wt))
      (2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dpois2(y, y,log=TRUE)*(wt/dispersion)))
    }
    
    f7<-function(b, y, x, mu, P, alpha = 0, wt = 1) {
      l2 <- length(y)
      l1 <- length(b)
      ltemp <- length(wt)
      
      # Construct diagonal weight matrix
      if (ltemp == 1) {
        Ptemp <- wt * diag(l2)
      } else {
        Ptemp <- diag(wt)
      }
      
      # Compute expected counts
      mu_i <- exp(alpha + x %*% b)
      
      # Initialize output matrix
      Pout <- matrix(0, nrow = l1, ncol = l1)
      
      # Loop over observations
      for (i in 1:l2) {
        xi <- x[i, , drop = FALSE]  # row as matrix
        weight_i <- Ptemp[i, i] * mu_i[i]
        Pout <- Pout + weight_i * (t(xi) %*% xi)
      }
      
      Pout
    }
    
    
    
    
  }

  ## Poisson / quasi-Poisson with identity link: lambda = Xb (no exp transformation).
  ## Used by dGamma(Inv_Dispersion=FALSE) + Poisson(identity) for DIC and logLik post-processing.
  ## Closures differ from the log-link versions in three ways:
  ##   f1/f2: lambda = alpha + x %*% b  (linear, not exp)
  ##   f3:    gradient = X'(wt*(1 - y/lambda)) + P(b-mu)
  ##   f7:    Fisher info = X' diag(wt/lambda) X  (inverted vs log-link wt*lambda)
  if((family$family=="poisson"||family$family=="quasipoisson") && family$link=="identity")
  {
    f1<-function(b,y,x,alpha=0,wt=1){
      lambda <- as.vector(alpha + x %*% b)
      -sum(dpois2(y, lambda, log=TRUE) * wt)
    }
    f2<-function(b,y,x,mu,P,alpha=0,wt=1){
      lambda <- as.vector(alpha + x %*% b)
      -sum(dpois2(y, lambda, log=TRUE) * wt) + 0.5 * t((b-mu)) %*% P %*% (b-mu)
    }
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      lambda <- as.vector(alpha + x %*% b)
      t(x) %*% ((1 - y/lambda) * wt) + P %*% (b-mu)
    }
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      (2*f1(b,y,x,alpha,wt/dispersion) + 2*sum(dpois2(y, y, log=TRUE) * (wt/dispersion)))
    }
    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
      lambda <- as.vector(alpha + x %*% b)
      l1 <- length(b)
      l2 <- length(y)
      ltemp <- length(wt)
      wt_vec <- if(ltemp == 1L) rep(wt, l2) else as.vector(wt)
      Pout <- matrix(0, nrow=l1, ncol=l1)
      for(i in seq_len(l2)){
        xi <- x[i, , drop=FALSE]
        Pout <- Pout + (wt_vec[i] / lambda[i]) * (t(xi) %*% xi)
      }
      Pout
    }
  }
  
  if(family$family %in%  c("binomial","quasibinomial") && family$link=="logit")
  {
    f1<-function(b,y,x,alpha=0,wt=1){
      lambda<-t(exp(alpha+x%*%b))
      p<-lambda/(1+lambda)
      -sum(dbinom(round(wt*y),round(wt), p,log=TRUE))
    }
    f2<-function(b,y,x,mu,P,alpha=0,wt=1){
      lambda<-t(exp(alpha+x%*%b))
      p<-lambda/(1+lambda)
      -sum(dbinom(round(wt*y),round(wt),p,log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
    }
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      p<-1/(1+t(exp(-alpha-x%*%b)))
      t(x)%*%((t(p)-y)*wt)+P%*%(b-mu)
    }
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      #dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dbinom(round(wt*y),round(wt),y,log=TRUE)))
      (2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dbinom(round((wt/dispersion)*y),round(wt/dispersion),y,log=TRUE)))
    }
    
    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
      
      l2<-length(y)
      l1<-length(b)
      ltemp<-length(wt)
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      
      p<-1/(1+t(exp(-alpha-x%*%b)))
      
      
      Pout<-matrix(0,nrow=l1,ncol=l1)
      
      for(i in 1:l2){
        Pout<-Pout+Ptemp[i,i]*p[i]*(1-p[i])*x[i,]%*%t(x[i,])
      }
      
      Pout
    }
    
    
  }
  
  if(family$family %in%  c("binomial","quasibinomial") && family$link=="probit")
  {
    f1<-function(b,y,x,alpha=0,wt=1){
      p<-pnorm(alpha+x%*%b)
      -sum(dbinom(round(wt*y),round(wt), p,log=TRUE))
    }
    f2<-function(b,y,x,mu,P,alpha=0,wt=1){
      p<-pnorm(alpha+x%*%b)
      -sum(dbinom(round(wt*y),round(wt),p,log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
    }
    
    
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      p1<-pnorm(alpha+x%*%b)
      p2<-pnorm(-alpha-x%*%b)
      -t(x)%*%as.matrix(((y*dnorm(alpha+x%*%b)/p1)-(1-y)*dnorm(alpha+x%*%b)/p2)*wt)+P%*%(b-mu)
    }
    
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      #dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dbinom(round(wt*y),round(wt),y,log=TRUE)))
      (2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dbinom(round((wt/dispersion)*y),round(wt/dispersion),y,log=TRUE)))
    }
    
    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
      
      l2<-length(y)
      l1<-length(b)
      ltemp<-length(wt)
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      
      p<-1/(1+t(exp(-alpha-x%*%b)))
      
      
      Pout<-matrix(0,nrow=l1,ncol=l1)
      
      for(i in 1:l2){
        Pout<-Pout+Ptemp[i,i]*p[i]*(1-p[i])*x[i,]%*%t(x[i,])
      }
      
      Pout
    }
    
    
    
    
    
  }
  if(family$family %in%  c("binomial","quasibinomial") && family$link=="cloglog")
  {
    f1<-function(b,y,x,alpha=0,wt=1){
      Xb<-alpha+x%*%b
      p<-1-exp(-exp(Xb))
      -sum(dbinom(round(wt*y),round(wt), p,log=TRUE))
    }
    f2<-function(b,y,x,mu,P,alpha=0,wt=1){
      Xb<-alpha+x%*%b
      p<-1-exp(-exp(Xb))
      -sum(dbinom(round(wt*y),round(wt),p,log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
    }
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      Xb<-alpha+x%*%b
      p1<-1-exp(-exp(Xb))
      p2<-exp(-exp(Xb))
      atemp<-exp(Xb-exp(Xb))
      -t(x)%*%as.matrix(((y*atemp/p1)-(1-y)*atemp/p2)*wt)+P%*%(b-mu)
    }
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dbinom(round(wt*y),round(wt),y,log=TRUE)))
      (2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dbinom(round((wt/dispersion)*y),round(wt/dispersion),y,log=TRUE)))
    }
    
    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
      
      l2<-length(y)
      l1<-length(b)
      ltemp<-length(wt)
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      
      p<-1/(1+t(exp(-alpha-x%*%b)))
      
      
      Pout<-matrix(0,nrow=l1,ncol=l1)
      
      for(i in 1:l2){
        Pout<-Pout+Ptemp[i,i]*p[i]*(1-p[i])*x[i,]%*%t(x[i,])
      }
      
      Pout
    }
    
    
    
  }
  ## Binomial / quasi-Binomial with identity link: theta = Xb in (0,1).
  ## Coefficient is the probability directly (not log-odds).
  if (family$family %in% c("binomial", "quasibinomial") && family$link == "identity") {

    f1 <- function(b, y, x, alpha = 0, wt = 1) {
      theta <- as.vector(alpha + x %*% b)
      -sum(stats::dbinom(round(wt * y), round(wt), theta, log = TRUE))
    }

    f2 <- function(b, y, x, mu, P, alpha = 0, wt = 1) {
      theta <- as.vector(alpha + x %*% b)
      -sum(stats::dbinom(round(wt * y), round(wt), theta, log = TRUE)) +
        0.5 * t(b - mu) %*% P %*% (b - mu)
    }

    ## Gradient of f2 w.r.t. b.
    ## ∂(-log L)/∂b = X' diag(wt) (θ − y) / (θ(1−θ)),  with θ = Xb.
    f3 <- function(b, y, x, mu, P, alpha = 0, wt = 1) {
      theta <- as.vector(alpha + x %*% b)
      t(x) %*% ((theta - y) / (theta * (1 - theta)) * wt) + P %*% (b - mu)
    }

    ## Deviance: 2 × (fitted NLL − saturated NLL).
    f4 <- function(b, y, x, alpha = 0, wt = 1, dispersion = 1) {
      (2 * f1(b, y, x, alpha, wt / dispersion) +
         2 * sum(stats::dbinom(round((wt / dispersion) * y),
                               round(wt / dispersion), y, log = TRUE)))
    }

    ## Fisher information: I(b) = X' diag(wt / (θ(1−θ))) X.
    f7 <- function(b, y, x, mu, P, alpha = 0, wt = 1) {
      theta  <- as.vector(alpha + x %*% b)
      l1     <- length(b)
      l2     <- length(y)
      wt_vec <- if (length(wt) == 1L) rep(wt, l2) else as.vector(wt)
      Pout   <- matrix(0, nrow = l1, ncol = l1)
      for (i in seq_len(l2)) {
        xi   <- x[i, , drop = FALSE]
        Pout <- Pout + (wt_vec[i] / (theta[i] * (1 - theta[i]))) * (t(xi) %*% xi)
      }
      Pout
    }

  }

  if(family$family=="Gamma" && family$link=="log")
  {
    
    f1<-function(b,y,x,alpha=0,wt=1){
      mu<-t(exp(alpha+x%*%b))
      disp2<-1/wt
      -sum(dgamma(y,shape=1/disp2,scale=mu*disp2,log=TRUE))
    }
    
    #  f1<-f1_gamma
    f2 <- function(b, y, x, mu, P, alpha = 0, wt = 1) {
      
      eta    <- alpha + x %*% b
      scale2 <- t(exp(eta - log(wt)))
      disp2  <- 1 / wt
      
      -sum(dgamma(y, shape = 1/disp2, scale = scale2, log = TRUE)) +
        0.5 * t((b - mu)) %*% P %*% (b - mu)
    }
    
    
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      mu2<-t(exp(alpha+x%*%b))
      t(x)%*%(t(1-y/mu2)*wt)+P%*%(b-mu)
    }
    
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      disp2<-1/(wt/dispersion)
      
      (2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dgamma(y,shape=1/disp2,scale=y*disp2,log=TRUE)))
      
    }
    
    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
      
      l2<-length(y)
      l1<-length(b)
      ltemp<-length(wt)
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      
      p<-1/(1+t(exp(-alpha-x%*%b)))
      
      
      Pout<-matrix(0,nrow=l1,ncol=l1)
      
      for(i in 1:l2){
        Pout<-Pout+Ptemp[i,i]*p[i]*(1-p[i])*x[i,]%*%t(x[i,])
      }
      
      Pout
    }
    
    
    
  }	

  ## Gamma with identity link: coefficient IS the Gamma rate β.
  ## Y_i | β ~ Gamma(shape = k, rate = β),  mean = k/β,  Var = k/β².
  ## lik_shape k is the known Gamma shape parameter (default 1 = exponential).
  if (family$family == "Gamma" && family$link == "identity") {

    k <- lik_shape   # capture in closures

    f1 <- function(b, y, x, alpha = 0, wt = 1) {
      beta <- as.vector(alpha + x %*% b)   # rate = Xb
      -sum(stats::dgamma(y, shape = k, rate = beta, log = TRUE) * wt)
    }

    f2 <- function(b, y, x, mu, P, alpha = 0, wt = 1) {
      beta <- as.vector(alpha + x %*% b)
      -sum(stats::dgamma(y, shape = k, rate = beta, log = TRUE) * wt) +
        0.5 * t(b - mu) %*% P %*% (b - mu)
    }

    ## Gradient of f2 w.r.t. b.
    ## ∂(-log L)/∂b  =  X' diag(wt) (y - k/β),  with β = Xb.
    f3 <- function(b, y, x, mu, P, alpha = 0, wt = 1) {
      beta <- as.vector(alpha + x %*% b)
      t(x) %*% ((y - k / beta) * wt) + P %*% (b - mu)
    }

    ## Deviance: 2 × (saturated NLL − fitted NLL).
    ## Saturated rate for obs i: β_sat,i = k / y_i  →  dgamma(y_i, k, k/y_i).
    f4 <- function(b, y, x, alpha = 0, wt = 1, dispersion = 1) {
      wt_eff <- wt / dispersion
      2 * f1(b, y, x, alpha, wt_eff) +
        2 * sum(stats::dgamma(y, shape = k, rate = k / y, log = TRUE) * wt_eff)
    }

    ## Fisher information: I(b) = X' diag(wt × k/β²) X.
    f7 <- function(b, y, x, mu, P, alpha = 0, wt = 1) {
      beta   <- as.vector(alpha + x %*% b)
      l1     <- length(b)
      l2     <- length(y)
      wt_vec <- if (length(wt) == 1L) rep(wt, l2) else as.vector(wt)
      Pout   <- matrix(0, nrow = l1, ncol = l1)
      for (i in seq_len(l2)) {
        xi   <- x[i, , drop = FALSE]
        Pout <- Pout + (wt_vec[i] * k / beta[i]^2) * (t(xi) %*% xi)
      }
      Pout
    }

  }

  out=list(f1=f1,f2=f2,f3=f3,f4=f4,
           #f5=f5,
           #f6=f6,
           f7=f7)
  
  out$call<-match.call()
  
  class(out)<-"glmbfamfunc"
  out
  
  
}


## print.* S3 method provided by glmbayesCore (Imports); local registration commented out.
# #' @rdname glmbfamfunc
# #' @order 2
# #' @method print glmbfamfunc
# #' @export
# 
# 
# print.glmbfamfunc<-function(x,...)
# {
#   cat("Call:\n")
#   print(x$call)
#   cat("\nNegative Log-Likelihood Function:\nf1<-")
#   print(x$f1)
#   cat("\n(Negative Log-Posterior Function:\nf2<-")
#   print(x$f2)
#   cat("\nNegative Log-Posterior Gradient Function:\nf3<-")
#   print(x$f3)
#   cat("\nDeviance Function:\nf4<-")
#   print(x$f4)
#   
#   
# }


#' Standardize A Non-Gaussian Model
#'
#' Standardizes a Non-Gaussian Model prior to Envelope Creation
#' @param y a vector of observations of length m
#' @param x a design matrix of dimension m*p
#' @param P a positive-definite symmetric matrix specifying the prior precision matrix of the variables.
#' @param bstar a matrix containing the posterior mode from an optimization step
#' @param A1 a matrix containing the posterior precision matrix at the posterior mode
#' @return A list with the following components
#' \item{bstar2}{Standardized Posterior Mode}
#' \item{A}{Standardized Data Precision Matrix}
#' \item{x2}{Standardized Design Matrix}
#' \item{mu2}{Standardized Prior Mean vector}
#' \item{P2}{Standardized Precision Matrix Added to log-likelihood}
#' \item{L2Inv}{A matrix used when undoing the first step in standardization described below}
#' \item{L3Inv}{A matrix used when undoing the second step in standardization described below}
#' @details This functions starts with basic information about the model in the argument list and then
#' uses the following steps to further standardize the model (the model is already assumed to have a 0 prior mean vector
#' when this step is applied).
#' 
#' 1) An eigenvalue composition is applied to the posterior precision matrix, and the model is (as an interim step)
#' standardized to have a posterior precision matrix equal to the identity matrix. Please note that this means
#' that the prior precision matrix after this step is \code{"smaller"} than the identity matrix.
#' 
#' 
#' 2) A diagonal matrix epsilon is pulled out from the standardized prior precision matrix so that the remaining
#' part of the prior precision matrix still is positive definite. That part is then treated as part of the likelihood
#' for the rest of the standardization and simulation and only the part connected to epsilon is treated as part of the prior. 
#' Note that the exact epsilon chosen seems not to matter. Hence there are many possible ways of doing this 
#' standardization and future versions of this package may tweak the current approach 
#' if it helps improve numerical accuracy or acceptance rates.
#' 
#' 
#' 3) The model is next standardized (using a second eigenvalue decomposition) so that the prior (i.e., the portion connected to epsilon) is the identity 
#' matrix. The standardized model then simutaneously has the feature that the prior precision matrix is the 
#' identity matrix and that the data precision A (at the posterior mode) is a diagonal matrix. Hence the variables
#' in the standardized model are approximately independent at the posterior mode.
#' 
#' The steps here are based on the procedure described in \insertCite{Nygren2006}{glmbayes}.
#'
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_glmb_Standardize_Model.R
#' @export


glmb_Standardize_Model<-function(y, x, P, bstar, A1){
  
  return(.glmb_Standardize_Model_cpp(y, x, P, bstar, A1))
  
}


#' Envelope Sizing and Optimization
#'
#' \code{EnvelopeSize()} is the high-level entry point that constructs
#' per-dimension grids and expected draw counts, while \code{EnvelopeOpt()}
#' performs the adaptive optimization used when \code{Gridtype = 2}.
#'
#' These functions implement the grid sizing logic used in envelope construction
#' for rejection sampling. They make use of the theory described in
#' \insertCite{Nygren2006}{glmbayes} and the general implementation outlined in
#' \insertCite{glmbayesSimmethods}{glmbayes}. 
#' 
#' @param a Numeric vector of diagonal precisions for the log-likelihood
#'   (posterior precision is \eqn{1 + a_i}).
#' @param G1 Numeric matrix of candidate grid points (3 * l1).
#' @param Gridtype Integer code controlling grid sizing logic:
#'   \itemize{
#'     \item 1 = static threshold test
#'     \item 2 = adaptive optimization via \code{EnvelopeOpt()}
#'     \item 3 = always three-point grid
#'     \item 4 = always single-point grid
#'   }
#' @param n Integer; number of posterior draws to generate (used for grid sizing).
#' @param n_envopt Integer; effective sample size passed to \code{EnvelopeOpt}.
#'   Defaults to -1, which means "use n".
#' @param use_opencl Logical; if \code{TRUE}, attempt GPU acceleration.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#'
#' @param a1 Numeric vector of diagonal elements of the data precision matrix
#'   (used by \code{EnvelopeOpt}).
#' @param core_cnt Integer; number of OpenCL cores or parallel workers available
#'   (default 1). When >1, envelope build cost is scaled down to reflect
#'   parallel construction.
#'
#' @section Gridtype Logic and Candidates per Draw:
#'
#' The envelope sizing logic follows the analysis of \insertCite{Nygren2006}{glmbayes}.
#'
#' \describe{
#'
#'   \item{Gridtype 1: Static Threshold}{
#'     For each dimension \eqn{i}, if
#'     \eqn{\sqrt{1 + a_i} \leq 2/\sqrt{\pi} \approx 1.128379},
#'     then a single tangent at the posterior mode suffices.
#'     Expected candidates per draw in that dimension:
#'     \eqn{\sqrt{1 + a_i}}. Otherwise, a symmetric three-point envelope is used at
#'     \eqn{(\theta^\star_i - \omega_i, \theta^\star_i, \theta^\star_i + \omega_i)},
#'     with expected candidates per draw bounded above by
#'     \eqn{2/\sqrt{\pi}}.
#'   }
#'
#'   \item{Gridtype 2: Adaptive Optimization}{
#'     Each dimension is assigned either a single-point or three-point envelope
#'     by minimizing
#'     \deqn{T_\mathrm{total}(g_i) = T_\mathrm{build}(g_i) + T_\mathrm{sample}(n, acc_i(g_i)).}
#'     The optimizer balances build cost (grows with number of tangents) against
#'     sampling cost (decreases as acceptance improves).
#'     Expected candidates per draw:
#'     \eqn{\prod_j \mathrm{scaleest}_{i,j}}, where each factor is either
#'     \eqn{\sqrt{1+a_j}} (single-point) or \eqn{2/\sqrt{\pi}} (three-point),
#'     depending on the optimization outcome.
#'   }
#'
#'   \item{Gridtype 3: Always Three-Point}{
#'     Every dimension uses a symmetric three-point envelope.  
#'     Expected candidates per draw:
#'     \deqn{\left(\tfrac{2}{\sqrt{\pi}}\right)^k}
#'     for \eqn{k} dimensions, as shown in Theorem 3 of
#'     \insertCite{Nygren2006}{glmbayes}.
#'   }
#'
#'   \item{Gridtype 4: Always Single-Point}{
#'     Every dimension uses a single tangent at the posterior mode.  
#'     Expected candidates per draw:
#'     \deqn{\prod_{i=1}^k \sqrt{1 + a_i}}
#'     (Example 1 in \insertCite{Nygren2006}{glmbayes}).
#'   }
#'
#' }
#'
#' @details
#' \code{EnvelopeSize()} returns the constructed grid (\code{G2}),
#' index vectors (\code{GIndex1}), expected draw count (\code{E_draws}),
#' and the per-dimension grid index.  
#'
#' \code{EnvelopeOpt()} implements the adaptive optimization used in
#' \code{Gridtype = 2}, ranking dimensions by posterior variance and
#' promoting them to three-point tangents when the tradeoff is favorable.
#'
#' @return
#' \describe{
#'   \item{\code{EnvelopeSize()}}{A list with components \code{G2}, \code{GIndex1},
#'   \code{E_draws}, and \code{gridindex}.}
#'   \item{\code{EnvelopeOpt()}}{An integer vector of length \eqn{l1} with entries
#'   1 (single-point) or 3 (three-point).}
#' }
#'
#' @seealso \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeEval}}, \code{\link{EnvelopeSort}};
#' \code{\link{rNormal_reg}}, \code{\link{rglmb}} for user-facing sampling that uses these grids.
#' Vignettes: \insertCite{glmbayesSimmethods,glmbayesChapterA08}{glmbayes}.
#'
#' @references
#' \insertAllCited{}
#' @example inst/examples/Ex_EnvelopeOpt.R
#' @export
#' @usage EnvelopeSize(a, G1, Gridtype = 2L, n = 1000L, n_envopt = -1,
#'                     use_opencl = FALSE, verbose = FALSE)
#' @rdname EnvelopeSize
#' @export
EnvelopeSize <- function(a,
                         G1,
                         Gridtype   = 2L,
                         n          = 1000L,   # <-- updated default
                         n_envopt   = -1,
                         use_opencl = FALSE,
                         verbose    = FALSE) {
  .EnvelopeSize_cpp(a, G1, Gridtype, n, n_envopt, use_opencl, verbose)
}


#' @usage EnvelopeOpt(a1,n,core_cnt=1L)
#' @rdname EnvelopeSize
#' @export


EnvelopeOpt<-function(a1,n,core_cnt=1L){
  
  core_cnt <- as.integer(core_cnt)
  if (is.na(core_cnt) || core_cnt < 1L) core_cnt <- 1L
  
  
  a1rank<-rank(1/(1+a1))
  l1<-length(a1)
  
  dimcount<-matrix(0,(l1+1),l1)
  scaleest<-matrix(0,(l1+1),l1)
  intest<-c(1:(l1+1))
  slopeest<-c(1:(l1+1))
  
  dimcount[1,]<-diag(diag(l1))
  scaleest[1,]<-sqrt(1+a1)
  slopeest[1]<-prod(scaleest[1,])
  
  for(i in 2:(l1+1)){
    dimcount[i,]<-dimcount[i-1,]
    scaleest[i,]<-scaleest[i-1,]
    for(j in 1:l1){
      if(a1rank[j]==i-1){ 
        dimcount[i,j]<-3
        scaleest[i,j]<-2/sqrt(pi) 
      }
    }
    ##    intest[i]<-3^(i-1)
    intest[i]<-(3^(i-1))
    slopeest[i]<-prod(scaleest[i,])
  }
  evalest<-(intest/core_cnt)+n*slopeest
  minindex<-0
  for(j in 1:(l1+1)){if(evalest[j]==min(evalest)){minindex<-j}}
  
  ##  message("Estimated draws per Acceptance: ", slopeest[minindex])
  
  dimcount[minindex,]
  
}



#' GPU-Accelerated Envelope Construction for Posterior Simulation
#'
#' @name EnvelopeBuild
#'
#' @details
#' Constructs an enveloping function for posterior simulation using a grid of
#' tangency points. The envelope is used in accept-reject sampling to guarantee
#' iid draws from the posterior distribution. The implementation follows
#' \insertCite{Nygren2006}{glmbayes}, with extensions for GPU acceleration
#' (via OpenCL), dynamic grid optimization, and parallelized evaluation.
#'
#' The envelope is typically built around the posterior mode \eqn{\theta^\star} for a model in standard 
#' form (which in this context means a model with a diagonal posterior precision matrix
#' and prior identity precision matrix - \code{glmb_Standardize_Model}). It uses dimension-specific width parameters 
#' \eqn{\omega_i} derived from the precision matrix. Tangency points are selected per dimension, and the full grid is
#' formed via Cartesian expansion. Negative log-likelihood and gradient values
#' are computed at each grid point, either on CPU or GPU depending on the
#' \code{use_opencl} flag. These values are used to construct a piecewise
#' envelope function that dominates the posterior density.
#'
#' @section Models in standard form:
#'
#' The standard-form restriction and its closed-form truncated-normal integrals
#' follow \insertCite{Nygren2006}{glmbayes}. See
#' \insertCite{glmbayesChapterA08}{glmbayes} for the full theoretical details (standard form restriction,
#' closed-form truncated-normal integrals, and the resulting log-scale tractability).
#'
#' In the implementation, these standard-form quantities determine the grid-based
#' tangency shifts and the precomputed region constants (e.g., the log-CDF pieces)
#' that drive the mixture weights used by the envelope sampler.

#' @section Construction of restricted subgradient densities:
#'
#' For the full restricted density construction and the resulting envelope
#' constants, see \insertCite{glmbayesChapterA08}{glmbayes}.
#'
#' In the implementation, these theory objects become the precomputed
#' region log-constants (via closed-form CDF pieces) that are used to normalize
#' the envelope mixture weights for the accept-reject sampler.
#' 
#' @section Mixture construction and tractable probabilities:
#'
#' The mixture construction and its tractable region probabilities are derived
#' in \insertCite{glmbayesChapterA08}{glmbayes}. In the implementation, these theory objects become the
#' precomputed region log-constants and the mixture weights (`PLSD`) used by
#' the envelope-based accept-reject sampler.
#' @section Log-scale properties of the envelope function:
#'
#' The log-scale form of the envelope factor and the subgradient inequality that
#' imply envelope dominance are given in \insertCite{glmbayesChapterA08}{glmbayes}. In the implementation, these
#' properties allow pointwise evaluation in the log-domain and provide the
#' theoretical basis for the rejection test inside the sampler.
#' @section Use of the envelope during sampling:
#'
#' The standardized sampler called through \code{.rNormalGLM_std_cpp()}
#' uses the envelope to generate posterior samples via rejection sampling. Although not exported,
#' this routine is called internally by \code{.rNormalGLM_cpp()}, which in turn is invoked by
#' the user-facing function \code{rNormal_reg()}. Together, these routines implement
#' envelope-based sampling for generalized linear models with log-concave likelihood functions
#' and multivariate normal priors.
#'
#' The envelope provides a mixture of restricted likelihood-subgradient densities,
#' each defined over a region \eqn{A_i}, with associated mixture weights
#' \eqn{\tilde{p}_i} stored in \code{PLSD}. The sampling proceeds as follows:
#'
#' \enumerate{
#'   \item A region index \eqn{J(i)} is drawn from the discrete distribution
#'         defined by \code{PLSD}.
#'   \item A candidate \eqn{\theta_i} is drawn from the restricted density
#'         \eqn{q^{\bar{\theta}_{J(i)}}_{A_{J(i)}}}, using the normal CDF bounds
#'         \code{loglt} and \code{logrt}, and subgradient vector \code{cbars}.
#'         Simulation for each dimension uses the internal C++ function
#'         \code{ctrnorm_cpp()}, which explicitly uses these inputs.
#'   \item The log-likelihood \eqn{\log f(y \mid \theta_i)} is computed and
#'         stored in \code{testll[0]} using the appropriate likelihood function
#'         \code{f2}.
#' }
#'
#' The acceptance test is performed using the inequality
#' \deqn{
#'   \log(U_2) \le \mathrm{LLconst}[J(i)] + \mathrm{cbars}[J(i), ]^{T} \theta_i
#'   + \log f(y \mid \theta_i),
#' }
#' which is equivalent to
#' \deqn{
#'   \log(U_2) \le \log f(y \mid \theta_i) - \left( \log f(y \mid \bar{\theta}_{J(i)}) - c(\bar{\theta}_{J(i)})^{T}(\theta_i - \bar{\theta}_{J(i)}) \right),
#' }
#' 
#' where:
#' \itemize{
#'   \item \code{LLconst[J(i)]} stores the precomputed quantity
#'         \eqn{-\log f(y \mid \bar{\theta}_{J(i)}) - c(\bar{\theta}_{J(i)})^{T} \bar{\theta}_{J(i)}},
#'         computed during envelope construction via \code{EnvelopeSet_LogP_C2()}.
#'   \item \code{cbars[J(i), ]} is the precomputed subgradient vector
#'         \eqn{c(\bar{\theta}_{J(i)})}, extracted via \code{cbars(J(i), _)}.
#'         It defines the exponential tilt direction used to evaluate the envelope.
#'   \item \code{testll[0]} is the log-likelihood at the candidate draw
#'         \eqn{\theta_i}, evaluated using the model specified by
#'         \code{family} and \code{link}.
#'   \item \eqn{-\log(U_2)} is the threshold from a uniform draw
#'         \eqn{U_2 \sim \mathrm{Unif}(0,1)}.
#' }
#'
#' The right-hand side of this inequality is always non-positive, and equals zero
#' when \eqn{\theta_i = \bar{\theta}_{J(i)}}. This reflects the fact that the envelope
#' is tangent to the log-likelihood at each \eqn{\bar{\theta}_j}, and lies above it elsewhere.
#'
#' This procedure guarantees that accepted samples are drawn from the posterior
#' \eqn{\pi(\theta \mid y)}. The envelope ensures bounded rejection probability,
#' and the mixture structure allows efficient sampling across regions. The output
#' \code{out} contains accepted draws, and \code{draws} records the number of
#' attempts per sample.
#'
#' The components returned by \code{EnvelopeBuild()} are used in specific steps of the
#' sampling procedure as follows:
#' \itemize{
#'   \item \code{PLSD} is used to randomly select a region index \eqn{J(i)} from the envelope mixture.
#'   \item \code{loglt} and \code{logrt} define the truncated normal bounds for each dimension,
#'         used together with \code{cbars} to generate candidate values \eqn{\theta_i}.
#'   \item \code{cbars} provides the subgradient vectors \eqn{c(\bar{\theta}_j)} used both for
#'         candidate generation and for computing the acceptance test.
#'   \item \code{LLconst} stores precomputed constants used in the acceptance inequality,
#'         avoiding recomputation of posterior terms at tangency points.
#'   \item \code{logU} stores the per-dimension log-density contributions for each region,
#'         computed during envelope setup. These values are summed to produce \code{logP},
#'         which determines the mixture weights \code{PLSD}.
#'   \item \code{logP} contains the total log-probabilities for each grid component,
#'         which are normalized to form the mixture weights \code{PLSD}.
#'   \item \code{thetabars} stores the tangency points \eqn{\bar{\theta}_j} used to define
#'         subgradients and region-specific densities.
#'   \item \code{GridIndex} encodes the sampling type (tail, center, line) used for each
#'         dimension and region, guiding how each coordinate is simulated.
#' }
#' @section Algorithmic steps (linked to theory):
#'
#' The implementation of \code{EnvelopeBuild} follows the envelope construction
#' in \insertCite{Nygren2006}{glmbayes} for models in standard form (see Section 3--3.3 there). 
#' Each computational step corresponds to a theoretical guarantee:
#'
#' 1. **Compute width parameters \eqn{\omega_i} from the diagonal precision matrix.**  In particular,
#' let \eqn{\theta^{\ast}} denote the unique posterior mode. For each dimension \eqn{i},
#'  define
#' \deqn{
#'   \omega_{i} :=
#'   \frac{\sqrt{2} - \exp\!\big(-1.20491 - 0.7321\,\sqrt{0.5 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2}}\big)}
#'        {\sqrt{1 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2}}}.
#' }
#' 
#'    As seen from the above, the widths \eqn{\omega_i} are derived from the local curvature of the 
#'    log-likelihood at the posterior mode. This ensures that the three-interval construction per
#'    dimension below yields an envelope whose efficiency does not deteriorate with sample size.
#'
#' 2. **Use the width parameters to construct intervals around the posterior mode \eqn{\theta^\star}.**  Specifically,
#' we set
#' \deqn{
#'   \ell_{i,1} = \theta^{\ast}_{i} - 0.5\,\omega_{i}, \quad
#'   \ell_{i,2} = \theta^{\ast}_{i} + 0.5\,\omega_{i},
#' }
#' and construct three intervals per dimension:
#' \deqn{
#'   A_{i,1} = (-\infty,\ell_{i,1}), \quad
#'   A_{i,2} = [\ell_{i,1},\ell_{i,2}], \quad
#'   A_{i,3} = (\ell_{i,2},\infty).
#' }
#'
#' For each dimension \eqn{i}, let \eqn{J_{i} = \{1,2,3\}} and define
#' \eqn{J = \prod_{i=1}^{p} J_{i}}, which has \eqn{3^{p}} elements. Each
#' \eqn{j \in J} is a vector \eqn{(j_{1},\ldots,j_{p})}, and we define
#' \deqn{
#'   A^{\ast}_{j} = \prod_{i=1}^{p} A_{i,j_{i}}.
#' }
#' The collection \eqn{A^{\ast} = \{A^{\ast}_{j} : j \in J\}} forms a partition of \eqn{\Theta}.
#'
#'    
#' 3. **For each member of the partition, select tangency points \eqn{\theta^\star \pm \omega_i}.**
#' 
#'  For each \eqn{j \in J}, define index sets
#' \deqn{
#'   C_{j1} = \{i : j_{i} = 1\}, \quad
#'   C_{j2} = \{i : j_{i} = 2\}, \quad
#'   C_{j3} = \{i : j_{i} = 3\}.
#' }
#' 
#' The tangency points \eqn{\bar{\theta}_{j}} are then defined componentwise by
#' \deqn{
#'   \bar{\theta}_{j,i} =
#'   \begin{cases}
#'     \theta^{\ast}_{i} - \omega_{i}, & i \in C_{j1}, \\
#'     \theta^{\ast}_{i},              & i \in C_{j2}, \\
#'     \theta^{\ast}_{i} + \omega_{i}, & i \in C_{j3}.
#'   \end{cases}
#' }
#'   
#'    The tangency points are hence chosen so that the envelope touches the log-likelihood at
#'    representative points in each interval, guaranteeing dominance and tightness.
#'
#' 4. **Build the full grid of tangency points (Cartesian product across dimensions).**  
#'    
#'    The Cartesian product of per-dimension partitions yields the \eqn{3^p} restricted
#'    densities described in the paper, ensuring coverage of the full parameter space.
#'
#' 5. **Evaluate negative log-likelihood and gradients at each grid point to construct the 
#'    likelihood subgradient densities and to facilitate accept rejection sampling**  
#'      
#'    The subgradients \eqn{c(\bar{\theta})} enter the likelihood-subgradient density construction
#'    (\insertCite{Nygren2006}{glmbayes}; see also \insertCite{glmbayesChapterA08}{glmbayes}), and both subgradients and
#'    negative log-likelihoods (through \eqn{h_{\bar{\theta}}(\cdot)}) are used in the accept-reject procedure. 
#'    CPU and GPU routines compute these values efficiently.
#'    - On CPU: via \code{f2_f3_non_opencl}.
#'    - On GPU: via \code{f2_f3_opencl}, which computes these in parallel across faces
#'    
#' 6. **Call \code{EnvelopeSet_Grid_C2_pointwise} to evaluate restricted multivariate normal log-densities.**  
#'    Each restricted density corresponds to a subset of the partition, normalized
#'    as in Remark 5 of \insertCite{Nygren2006}{glmbayes}.
#'
#' 7. **Call \code{EnvelopeSet_LogP_C2} to compute component log-probabilities and constants.**  
#'    The constants \eqn{\tilde{a}} and mixture weights \eqn{\tilde{p}_i} are computed
#'    explicitly as in Remark 6 of the paper, ensuring that the mixture envelope is properly normalized.
#'
#' 8. **Normalize probabilities (\code{PLSD}) and optionally sort grid components.**  
#'    Normalization implements Claim 2 of the paper so the mixture forms a valid
#'    dominating density for the posterior. Sorting is an implementation detail to
#'    improve sampling efficiency.
#'
#' @section Theory reference (JASA paper and vignette):
#' Definitions, claims, theorems, remarks, and examples through Remark 16
#' (including standard form, the \eqn{3^p} partition, and sampling remarks) are in
#' \insertCite{Nygren2006}{glmbayes}. An expanded narrative is in
#' \code{vignette("Chapter-A08", package = "glmbayes")}.
#'
#' @section Subgradient density formulation:
#' Each grid component corresponds to a tilted multivariate normal density,
#' normalized using the moment-generating function (MGF). In the single-point
#' case, centered at the posterior mode \eqn{\theta^\star}, the density is:
#' \deqn{
#' f(\theta) = \frac{1}{(2\pi)^{p/2} |A|^{-1/2} \cdot \text{MGF}_A(c)} \exp\left( -\frac{1}{2} (\theta - \mu)^T A (\theta - \mu) + c^T (\theta - \theta^\star) \right)
#' }
#' where:
#' - \eqn{A} is the precision matrix,
#' - \eqn{\mu} is the prior mean vector,
#' - \eqn{c} is the gradient of the log-likelihood at \eqn{\theta^\star},
#' - \eqn{\text{MGF}_A(c)} is the moment-generating function:
#' \deqn{
#' \text{MGF}_A(c) = \exp\left( \frac{1}{2} c^T A^{-1} c \right)
#' }
#'
#' This closed-form density dominates the posterior locally and is used when
#' \code{Gridtype = 1}. For richer envelopes, multiple such components are
#' constructed at tangency points \eqn{\theta_j}, each with its own gradient
#' \eqn{c_j}, and combined into a mixture:
#' \deqn{
#' f_{\text{env}}(\theta) = \sum_{j=1}^{K} p_j f_j(\theta)
#' }
#' where the weights \eqn{p_j} are computed using log-CDF differences and constants:
#' \deqn{
#' \log p_j = \log \Phi(U_j) - \log \Phi(L_j) - \text{NegLL}_j + \text{LLconst}_j
#' }
#' @section Gridtype logic:
#' The \code{Gridtype} argument controls how many tangency points are used per dimension:
#' - 1: Threshold rule. If \eqn{1 + a_i \le 2/\sqrt{\pi}}, use a single-point envelope at the mode;
#'      otherwise use three points.
#' - 2: Dynamic optimization via \code{EnvelopeOpt}, which balances grid build cost and
#'      expected acceptance rate. Grid size is scaled by \code{n} and the number of
#'      OpenCL cores when GPU is enabled.
#' - 3: Always use three points per dimension.
#' - 4: Always use a single point (mode only).
#' @section Supported families and links:
#' The following families and link functions are supported:
#' - Binomial: logit, probit, cloglog
#' - Quasibinomial: logit, probit
#' - Poisson: log
#' - Quasipoisson: log
#' - Gamma: log
#' - Gaussian: identity
#'
#' GPU acceleration (\code{use_opencl = TRUE}) is available for all of the above
#' except Gaussian, which is always evaluated on CPU.
#' @section GPU acceleration:
#' When \code{use_opencl = TRUE}, likelihood and gradient evaluations are
#' offloaded to the GPU using OpenCL. This can substantially reduce runtime for
#' high-dimensional models or large grids. Results are mathematically equivalent
#' to the CPU version, but small numerical differences may occur due to
#' floating-point arithmetic. If reproducibility across hardware is critical,
#' prefer the CPU path.
#'
#' If OpenCL support was not detected at compile time, the flag is ignored and
#' the CPU implementation is used. Diagnostic messages are printed when
#' \code{verbose = TRUE}.
#' @section Verbose output:
#' When \code{verbose = TRUE}, the function prints:
#' - Grid type, number of draws, OpenCL usage, and detected core count.
#' - Grid size after expansion.
#' - Time-stamped messages when entering the grid loop, starting likelihood
#'   evaluations, starting gradient evaluations, and invoking GPU kernels.
#' - Messages when setting grid values, computing log-probabilities, and sorting.
#' 
#' Any constants needed by the sampling are added to a list and returned.
#'
#' @param bStar     Point at which envelope should be centered (typically posterior mode).
#' @param A         Diagonal precision matrix for the log-likelihood in standard form.
#' @param y         A vector of observations of length \code{m}.
#' @param x A design matrix of dimension \code{m * p}.
#' @param mu        A vector giving the prior means of the variables.
#' @param P         Prior precision matrix of the variables (positive-definite).
#' @param alpha     Offset vector.
#' @param wt        A vector of weights.
#' @param family    Family for the envelope: \code{binomial}, \code{quasibinomial}, \code{poisson}, \code{quasipoisson}, or \code{Gamma}.
#' @param link      Link function ("logit", "probit", "cloglog" for binomial; "log" for Poisson/Gamma).
#' @param Gridtype  Method to determine the number of subgradient densities in the grid.
#' @param n         Number of draws from the posterior (used for grid sizing).
#' @param n_envopt Effective sample size passed to EnvelopeOpt for grid construction.
#'   Defaults to match `n`. Larger values encourage tighter envelopes.
#' @param sortgrid  Logical; if \code{TRUE}, sort the envelope descending by component probability.
#' @param use_opencl Logical; if \code{TRUE}, use OpenCL for gradient evaluations.
#' @param verbose   Logical; if \code{TRUE}, print progress messages.
#'
#' @param GridIndex A matrix indicating, for each grid component, whether the component
#'   lies in the left tail, center, or right tail of the density. Rows correspond to
#'   grid components; columns correspond to standardized variables.
#' @param cbars     A matrix containing the subgradient of the (adjusted) negative log-likelihood
#'   at each grid component.
#' @param Lint      A matrix storing the lower and upper bounds for each grid component,
#'   depending on whether sampling is from the left, center, or right.
#'
#' @param logP      A matrix (typically two columns) with information for each grid component.
#'   The first column usually holds the output from \code{EnvelopeSetGrid()}, corresponding to
#'   the restricted normal density.
#' @param NegLL     A vector of negative log-likelihood evaluations at each grid component.
#' @param G3        A matrix of tangency points used in the grid.
#'    
#' @return
#' \describe{
#'
#'   \item{\code{EnvelopeBuild()}}{A list of envelope components used for accept-reject sampling:
#'     \describe{
#'       \item{\code{GridIndex}}{Integer matrix encoding sampling type (tail, center, line) per dimension and region.}
#'       \item{\code{thetabars}}{Matrix of tangency points \eqn{\bar{\theta}_j} for each grid region.}
#'       \item{\code{cbars}}{Matrix of subgradients \eqn{c(\bar{\theta}_j)} of the negative log-likelihood at tangency.}
#'       \item{\code{loglt}}{Matrix of log left-tail probabilities per dimension and region.}
#'       \item{\code{logrt}}{Matrix of log right-tail probabilities per dimension and region.}
#'       \item{\code{logU}}{Matrix of selected per-dimension log-density contributions (tail/center) for each region.}
#'       \item{\code{logP}}{Matrix of total log-probabilities per region (first column); used to derive mixture weights.}
#'       \item{\code{PLSD}}{Vector of normalized mixture weights over grid regions used to draw region indices.}
#'       \item{\code{LLconst}}{Vector of acceptance-test constants per region used in the inequality for rejection sampling.}
#'     }
#'   }
#'
#'   \item{\code{EnvelopeSetGrid()}}{A list of matrices computed for grid-based log-density evaluation:
#'     \describe{
#'       \item{\code{Down}}{Lower bounds for truncated-normal evaluation per dimension and region.}
#'       \item{\code{Up}}{Upper bounds for truncated-normal evaluation per dimension and region.}
#'       \item{\code{lglt}}{Log left-tail probabilities (from \eqn{(-\infty, \mathrm{Up}]}) per dimension and region.}
#'       \item{\code{lgrt}}{Log right-tail probabilities (from \eqn{[\mathrm{Down}, \infty)}) per dimension and region.}
#'       \item{\code{lgct}}{Log central-interval probabilities (from \eqn{[\mathrm{Down}, \mathrm{Up}]}) per dimension and region.}
#'       \item{\code{logU}}{Selected log-probability per grid cell based on \code{GridIndex} (tail or center).}
#'       \item{\code{logP}}{Matrix with row-wise sums of \code{logU} (first column) used to form mixture weights.}
#'     }
#'   }
#'
#'   \item{\code{EnvelopeSetLogP()}}{A list with updated mixture-weight and acceptance constants:
#'     \describe{
#'       \item{\code{logP}}{Input \code{logP} with its second column populated by the log of unnormalized visit probabilities per region (mixture denominators).}
#'       \item{\code{LLconst}}{Vector of acceptance constants \eqn{-\log f(y \mid \bar{\theta}_j) - c(\bar{\theta}_j)^{T}\bar{\theta}_j} used in the accept-reject test.}
#'     }
#'   }
#' }
#'
#' @seealso \code{\link{EnvelopeSize}}, \code{\link{EnvelopeEval}}, \code{\link{EnvelopeSort}},
#' \code{\link{glmb_Standardize_Model}}; \code{\link{rNormal_reg}}, \code{\link{rglmb}}, \code{\link{glmb}}.
#' Theory and vignettes: \insertCite{Nygren2006}{glmbayes};
#' \insertCite{glmbayesChapterA08,glmbayesSimmethods,glmbayesChapterA10,glmbayesChapter12}{glmbayes}.
#'
#' @references
#' \insertAllCited{}
#' @example inst/examples/Ex_EnvelopeBuild.R
#' @importFrom Rdpack reprompt
#' @usage EnvelopeBuild(bStar,A,y,x,mu,P,alpha,wt,family = "binomial",link = "logit", 
#' Gridtype = 2L,n = 1L,n_envopt=NULL,sortgrid = FALSE,use_opencl = FALSE,verbose = FALSE)
#' @rdname EnvelopeBuild
#' @export

EnvelopeBuild <- function(
    bStar, A, y, x, mu, P, alpha, wt,
    family     = "binomial",
    link       = "logit",
    Gridtype   = 2L,
    n          = 1L,
    n_envopt   = NULL,       # effective sample size for EnvelopeOpt
    sortgrid   = FALSE,
    use_opencl = FALSE,
    verbose    = FALSE
) {
  # normalize n_envopt: if not supplied, fall back to n
  if (is.null(n_envopt)) {
    n_envopt <- n
  }
  
  if (length(n_envopt) != 1L || is.na(n_envopt) || n_envopt < 1) {
    stop("`n_envopt` must be a positive integer scalar.")
  }
  # coerce safely to integer
  n_envopt <- as.integer(n_envopt)
  
  if (family == "gaussian") {
    return(.EnvelopeBuild_Ind_Normal_Gamma_cpp(
      bStar, A, y, x, mu, P, alpha, wt,
      family = family, link = link,
      Gridtype = Gridtype, n = n,    
      n_envopt  = n_envopt, 
      sortgrid = sortgrid,
      use_opencl = use_opencl,
      verbose    = verbose
    ))
  }
  
  .EnvelopeBuild_cpp(
    bStar, A, y, x, mu, P, alpha, wt,
    family    = family,
    link      = link,
    Gridtype  = Gridtype,
    n         = n,
    n_envopt  = n_envopt,
    sortgrid  = sortgrid,
    use_opencl = use_opencl,
    verbose    = verbose
  )
}



#' @usage EnvelopeSetGrid(GridIndex, cbars, Lint)
#' @rdname EnvelopeBuild
#' @export
EnvelopeSetGrid <- function(GridIndex, cbars, Lint) {
  .EnvelopeSet_Grid_cpp(GridIndex, cbars, Lint)
}



#' @usage EnvelopeSetLogP(logP, NegLL, cbars, G3)
#' @rdname EnvelopeBuild
#' @export
EnvelopeSetLogP <- function(logP, NegLL, cbars, G3) {
  .EnvelopeSet_LogP_cpp(logP, NegLL, cbars, G3)
}



#' Evaluate Negative Log-Likelihood and Gradients
#'
#' \code{EnvelopeEval()} evaluates the negative log-likelihood and gradients
#' at a grid of parameter values, optionally using OpenCL acceleration.
#'
#' The lower-level helpers `f2_f3_non_opencl` and `f2_f3_opencl`
#' are internal C++ kernels used by the CPU and OpenCL backends.
#' The internal routine `run_opencl_pilot` benchmarks OpenCL performance
#' on a pilot subset of the grid to estimate runtime before full evaluation.
#'
#' These functions implement the grid evaluation logic used in envelope
#' construction for rejection sampling. They make use of the theory described
#' in \insertCite{Nygren2006}{glmbayes} and the general implementation outlined
#' in \insertCite{glmbayesSimmethods}{glmbayes}.
#'
#' @param G4 Numeric matrix of parameter values (parameters * grid points).
#' @param y Numeric response vector.
#' @param x Numeric design matrix.
#' @param mu Numeric matrix of offsets or prior means.
#' @param P Numeric matrix representing the portion of the prior precision
#'   shifted into the likelihood.
#' @param alpha Numeric offset vector of length `m`
#' @param wt Numeric vector of weights.
#' @param family Character string; model family (e.g. \code{"gaussian"}).
#' @param link Character string; link function (e.g. \code{"identity"}).
#' @param use_opencl Logical; if \code{TRUE}, attempt OpenCL acceleration.
#' @param verbose Logical; if \code{TRUE}, print diagnostic output.
#' @details
#' The evaluation workflow has several layers:
#' **1. High-level dispatch (`EnvelopeEval`)**
#'
#' * `EnvelopeEval()` is the user-facing entry point. It accepts a grid of
#'   parameter values (`G4`) and the data (`y`, `x`, `mu`, `P`, `alpha`, `wt`).
#' * If the grid is large (>= 14 columns), it first calls
#'   `run_opencl_pilot` to benchmark OpenCL performance and optionally
#'      report estimated runtime.
#' * It then dispatches to either the CPU or GPU backend:
#'   - If `use_opencl = TRUE` and the family is not `"gaussian"`, it calls
#'     `f2_f3_opencl` (an internal C++ kernel).
#'   - Otherwise, it calls `f2_f3_non_opencl` (the CPU kernel).
#'   
#' **2. CPU backend (`f2_f3_non_opencl`)**
#'
#' * This function evaluates the negative log-likelihood and gradients using
#'   standard CPU routines.
#' * It inspects the `family` and `link` arguments and routes to the correct
#'   pair of kernels (`f2_*` for the likelihood, `f3_*` for the gradient).
#' * For example:
#'   - `"binomial"` with `"logit"` calls `f2_binomial_logit()` and
#'     `f3_binomial_logit()`.
#'   - `"poisson"` calls `f2_poisson()` and `f3_poisson()`.
#'   - `"gaussian"` calls `f2_gaussian()` and `f3_gaussian()`.
#' * These kernels ultimately rely on the same C math routines that R itself
#'   uses (from the `nmath`/`rmath` libraries), ensuring numerical consistency
#'   with base R functions like `dnorm`, `dpois`, etc.
#'
#' **3. GPU backend (`f2_f3_opencl`)**
#'
#' * This function mirrors the CPU backend but executes the likelihood and
#'   gradient calculations on an OpenCL device (GPU or CPU).
#' * It flattens the input matrices/vectors and allocates output buffers.
#' * It then constructs a full OpenCL program by concatenating:
#'   - a generic OpenCL support header (`OPENCL.CL`),
#'   - OpenCL ports of R's `rmath`, `nmath`, and `dpq` libraries,
#'   - and the family/link-specific kernel source (e.g.
#'     `f2_f3_binomial_logit.cl`).
#' * The resulting program is compiled and passed to a kernel runner
#'   (`f2_f3_kernel_runner`) which executes the likelihood and gradient
#'   calculations in parallel on the device.
#' * This ensures that the GPU backend produces results consistent with the
#'   CPU backend, but can scale to much larger grids efficiently.
#'
#' **4. Pilot timing (`run_opencl_pilot`)**
#'
#' * This helper runs a small subset of the grid through the OpenCL backend
#'   to estimate runtime.
#' * It is used by `EnvelopeEval()` to inform users (when `verbose = TRUE`)
#'   whether OpenCL acceleration is likely to be beneficial.
#'
#' **5. Returned values**
#'
#' * All backends return a list with:
#'   - `NegLL`: numeric vector of negative log-likelihood values.
#'   - `cbars`: numeric matrix of gradients (parameters * grid points).
#'   
#' **6. Role of likelihood and gradients in sampling**
#'
#' * The outputs of `EnvelopeEval()` - the negative log-likelihood values
#'   (`NegLL`) and the gradient matrix (`cbars`) - are not endpoints in
#'   themselves. They form the *envelope* used in the rejection sampler
#'   implemented by internal functions such as
#'   `.rNormalGLM_std_cpp()`.
#' * This routine is called by `.rNormalGLM_cpp()`, which underlies the
#'   user-facing function `rNormal_reg()`. Together they implement
#'   envelope-based posterior sampling for GLMs with log-concave likelihoods
#'   and multivariate normal priors.
#'
#' **7. Simulation execution (accept/reject procedure)**
#'
#' The acceptance test is performed using
#' \deqn{
#'   \log(U_2) \leq
#'   \log f(y \mid \theta_i) -
#'   \Big(\log f(y \mid \bar{\theta}_{J(i)}) -
#'        c(\bar{\theta}_{J(i)})^T(\theta_i - \bar{\theta}_{J(i)})\Big) \leq 0
#' }
#'
#' **Connections between code and notation:**
#' * The arguments `G4` (in `EnvelopeEval`) and `b` (in `f2_f3_*`) both
#'   represent the grid of tangency points \eqn{\bar{\theta}_j}.
#' * The output `NegLL` corresponds to
#'   \eqn{-\log f(y \mid \bar{\theta}_{J(i)})}, i.e. the negative
#'   log-likelihood evaluated at each tangency point.
#' * The output `cbars` corresponds to the subgradient vectors
#'   \eqn{c(\bar{\theta}_{J(i)})}, which define the tangent hyperplanes
#'   used in the envelope construction.
#'
#' **Precomputation for efficiency:**
#' * Both `NegLL` and `cbars` are computed once during envelope construction,
#'   prior to the simulation stage.
#' * This means the sampler does not need to recompute likelihoods or
#'   gradients at every candidate draw - it simply reuses the stored values
#'   (`NegLL`, `cbars`, and `LLconst`) in the acceptance inequality.
#'
#' This design ensures that the envelope is tangent to the log-likelihood at
#' each \eqn{\bar{\theta}_j}, lies above it elsewhere, and that the
#' accept-reject procedure can run efficiently while still producing samples
#' from the true posterior \eqn{\pi(\theta \mid y)}.
#'
#' @seealso \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeSize}}, \code{\link{EnvelopeSort}};
#' \code{\link{rNormal_reg}}, \code{\link{rglmb}}. Vignettes:
#' \insertCite{glmbayesSimmethods,glmbayesChapterA08,glmbayesChapterA10,glmbayesChapter12}{glmbayes}.
#'
#' @references
#' \insertAllCited{}
#' @return
#' \describe{
#'   \item{EnvelopeEval}{List with components \code{NegLL} (numeric vector of
#'   negative log-likelihood values) and \code{cbars} (numeric matrix of gradients).}
#'   \item{f2_f3_non_opencl}{List with components \code{qf} (negative log-likelihood)
#'   and \code{grad} (gradients) from the CPU kernel.}
#'   \item{f2_f3_opencl}{List with components \code{qf} and \code{grad} from the
#'   OpenCL kernel.}
#'   \item{run_opencl_pilot}{Numeric scalar giving estimated runtime (seconds)
#'   for OpenCL evaluation on a pilot subset of the grid.}
#' }
#'
#' @example inst/examples/Ex_EnvelopeEval.R
#'
#' @rdname EnvelopeEval
#' @export
#' @usage EnvelopeEval(G4, y, x, mu, P, alpha, wt,
#'                     family, link,
#'                     use_opencl = FALSE, verbose = FALSE)
EnvelopeEval <- function(G4, y, x, mu, P, alpha, wt,
                         family, link,
                         use_opencl = FALSE,
                         verbose = FALSE) {
  if (!is.matrix(G4)) stop("G4 must be a numeric matrix")
  if (!is.numeric(y)) stop("y must be numeric")
  if (!is.matrix(x)) stop("x must be a numeric matrix")
  if (!is.matrix(mu)) stop("mu must be a numeric matrix")
  if (!is.matrix(P)) stop("P must be a numeric matrix")
  if (!is.numeric(alpha)) stop("alpha must be numeric")
  if (!is.numeric(wt)) stop("wt must be numeric")
  if (!is.character(family) || length(family) != 1L) stop("family must be a string")
  if (!is.character(link) || length(link) != 1L) stop("link must be a string")
  
  
  .EnvelopeEval_cpp(G4, y, x, mu, P, alpha, wt,
               family, link,
               use_opencl, verbose)
}


#' Builds Dispersion-Aware Envelope for Simulation
#'
#' @description
#' Constructs a dispersion-aware envelope for simulation in Gaussian models with uncertain variance.
#' This function extrapolates the coefficient envelope across a high-probability interval for the
#' dispersion parameter \code{sigma^2}, and builds a global upper bound for the log-posterior remainder.
#' It also computes mixture weights for envelope faces and adjusts the Gamma proposal for precision.
#'
#' The envelope is constructed using the slopes of the face constants with respect to dispersion,
#' evaluated at an anchor point. The resulting structure supports exact i.i.d. sampling via accept-reject correction.
#'
#' The procedure follows these steps:
#' \enumerate{
#'
#'   \item \strong{Posterior precision (Gamma) parameters.}
#'   Using the prior and posterior-predictive RSS, set
#'
#'   \deqn{ \mathrm{shape2} = \mathrm{Shape} + n_{\mathrm{obs}}/2, \quad
#'          \mathrm{rate3}  = \mathrm{Rate} + \mathrm{RSS}_{\mathrm{post}}/2 }
#'
#'   These parameterize the posterior precision
#'   \eqn{v = 1/\sigma^2 \sim \mathrm{Gamma}(\mathrm{shape2}, \mathrm{rate3})}.
#'
#'   \item \strong{Central credible interval for dispersion (low, upp).}
#'   Choose a central mass level \code{max_disp_perc} (e.g., 0.99) for precision, then invert the corresponding
#'   Gamma quantiles to dispersion:
#'
#'   \deqn{ \mathrm{low} = 1 / Q_{\Gamma}(\mathrm{max\_disp\_perc}; \mathrm{shape2}, \mathrm{rate3}), \quad
#'          \mathrm{upp} = 1 / Q_{\Gamma}(1 - \mathrm{max\_disp\_perc}; \mathrm{shape2}, \mathrm{rate3}) }
#'
#'   The interval \eqn{[\mathrm{low}, \mathrm{upp}]} is the domain over which all envelopes must dominate.
#'
#'   \item \strong{Face slopes at an anchor (dispstar).}
#'
#'   \deqn{ \mathrm{dispstar} = \mathrm{rate3} / (\mathrm{shape2} - 1) }
#'
#'   (posterior mean of \eqn{\sigma^2}). Compute \eqn{\mathrm{New\_LL\_Slope}_j} for each face \eqn{j}.
#'
#'   \item \strong{Linear extrapolation of face constants.}
#'
#'   \deqn{ \theta^{\mathrm{low}}_{\bar{j}} = \theta^{\mathrm{base}}_{\bar{j}} + (\mathrm{low} - \mathrm{dispstar}) \cdot \mathrm{New\_LL\_Slope}_j }
#'
#'   \deqn{ \theta^{\mathrm{upp}}_{\bar{j}} = \theta^{\mathrm{base}}_{\bar{j}} + (\mathrm{upp} - \mathrm{dispstar}) \cdot \mathrm{New\_LL\_Slope}_j }
#'
#'   \item \strong{Global upper line and endpoint maxima.}
#'
#'   \deqn{ \mathrm{max\_low} = \max_j \theta^{\mathrm{low}}_{\bar{j}}, \quad
#'          \mathrm{max\_upp} = \max_j \theta^{\mathrm{upp}}_{\bar{j}} }
#'
#'   \deqn{ \mathrm{new\_slope} = (\mathrm{max\_upp} - \mathrm{max\_low}) / (\mathrm{upp} - \mathrm{low}), \quad
#'          \mathrm{new\_int}   = \mathrm{max\_low} - \mathrm{new\_slope} \cdot \mathrm{low} }
#'
#'   \item \strong{Face slack and mixture weights.}
#'
#'   \deqn{ \mathrm{lg\_prob\_factor}_j =
#'          \max\big(\theta^{\mathrm{upp}}_{\bar{j}} - \mathrm{max\_upp},\;
#'                    \theta^{\mathrm{low}}_{\bar{j}} - \mathrm{max\_low}\big) }
#'
#'   Combine with \eqn{\mathrm{New\_logP2}_j = \mathrm{logP}_j + \tfrac{1}{2}\|\bar{c}_j\|^2}
#'   to form mixture weights
#'   \eqn{\mathrm{PLSD}_j \propto \exp\!\bigl(\mathrm{New\_logP2}_j + \mathrm{lg\_prob\_factor}_j\bigr)}.
#'
#'   \item \strong{Gamma tilt and dispersion-axis envelope.}
#'
#'   \deqn{ \mathrm{dispstar} = (\mathrm{upp} - \mathrm{low}) / \log(\mathrm{upp}/\mathrm{low}) }
#'
#'   \deqn{ \mathrm{lm\_log2} = \mathrm{new\_slope} \cdot \mathrm{dispstar}, \quad
#'          \mathrm{lm\_log1} = \mathrm{new\_int} + \mathrm{new\_slope} \cdot \mathrm{dispstar}
#'                             - \mathrm{new\_slope} \cdot \log(\mathrm{dispstar}) }
#'
#'   Tilt the Gamma proposal via \eqn{\mathrm{shape3} = \mathrm{shape2} - \mathrm{lm\_log2}}.
#' }
#'
#' @param Env        Envelope object from \code{\link{EnvelopeBuild}}, containing tangency points and gradients
#' @param Shape      Prior shape parameter for precision \code{v = 1 / sigma^2}
#' @param Rate       Prior rate parameter for precision
#' @param mu         Prior mean parameter
#' @param P          Prior precision matrix for coefficients
#' @param y          Numeric response vector of length `m`
#' @param x          a design matrix of dimension \code{m * p}
#' @param alpha      Numeric offset vector of length `m`
#' @param n_obs      Number of observations
#' @param RSS_post   Expected posterior weighted residual sum of squares
#'   (i.e., \eqn{\mathbb{E}[\mathrm{RSS}(\beta)\mid y,\phi]} under the Normal
#'   posterior for \eqn{\beta} at fixed dispersion \eqn{\phi = d}). This value
#'   is used for dispersion anchoring / Gamma updates.
#' @param RSS_ML   Residual sum of squares associated with MLE estimate
#' @param max_disp_perc Truncation level for dispersion (default 0.99)
#' @param disp_lower lower bound truncation for dispersion 
#' @param disp_upper upper bound truncation for dispersion
#' @param verbose Option to have verbose output
#' @param wt          weight vector
#' @param use_parallel Logical. Whether to use parallel processing.
#'
#' @return
#' \describe{
#' \item{\code{EnvelopeDispersionBuild()}}{A list containing:
#' \describe{
#' \item{\code{Env_out}}{Envelope object with updated mixture weights (\code{PLSD})}
#' \item{\code{gamma_list}}{Posterior Gamma tilt parameters
#' \describe{
#' \item{\code{shape3}}{Adjusted shape parameter after slope correction}
#' \item{\code{rate2}}{Posterior rate parameter, defined as \code{Rate + rss_min_global/2}}
#' \item{\code{disp_upper}}{Upper bound of the dispersion interval \eqn{\sigma^2}}
#' \item{\code{disp_lower}}{Lower bound of the dispersion interval \eqn{\sigma^2}}
#' }}
#' \item{\code{UB_list}}{Upper-bound diagnostics
#' \describe{
#' \item{\code{RSS_ML}}{Residual sum of squares at the maximum-likelihood estimate}
#' \item{\code{RSS_Min}}{Minimum residual sum of squares across envelope faces}
#' \item{\code{max_New_LL_UB}}{Maximum extrapolated face constant at the upper dispersion bound}
#' \item{\code{max_LL_log_disp}}{Log-posterior upper bound evaluated at \code{disp_upper}}
#' \item{\code{lm_log1}}{Intercept term of the global upper line approximation}
#' \item{\code{lm_log2}}{Slope term of the global upper line approximation}
#' \item{\code{lg_prob_factor}}{Per-face slack factors used in mixture weighting}
#' \item{\code{lmc1}}{Linear extrapolation constant (intercept)}
#' \item{\code{lmc2}}{Linear extrapolation constant (slope)}
#' \item{\code{UB2min}}{Minimum UB2 value across faces, used for diagnostics}
#' }}
#' \item{\code{diagnostics}}{Internal diagnostic values
#' \describe{
#' \item{\code{dispstar}}{Anchor dispersion value (posterior mean or geometric mean)}
#' \item{\code{New_LL_Slope}}{Vector of slopes of face constants at \code{dispstar}}
#' \item{\code{shape2}}{Posterior shape parameter before tilt correction}
#' \item{\code{rate3}}{Posterior rate parameter before tilt correction}
#' \item{\code{shape3}}{Adjusted shape parameter (same as in \code{gamma_list})}
#' \item{\code{max_low}}{Maximum extrapolated face constant at the lower dispersion bound}
#' \item{\code{max_upp}}{Maximum extrapolated face constant at the upper dispersion bound}
#' \item{\code{new_slope}}{Slope of the global upper line across dispersion bounds}
#' \item{\code{new_int}}{Intercept of the global upper line across dispersion bounds}
#' \item{\code{prob_factor}}{Normalized mixture weights across faces}
#' \item{\code{UB2min}}{Minimum UB2 diagnostic value (duplicated for consistency)}
#' }}
#' }}
#' \item{\code{EnvBuildLinBound()}}{Numeric vector of slopes of face constants with respect to dispersion, evaluated at the anchor \code{dispstar}}
#' \item{\code{thetabar_const()}}{Numeric vector of base face constants computed from tangency points and gradient vectors under prior precision \code{P}}
#' \item{\code{Inv_f3_with_disp()}}{Numeric matrix of inverse function evaluations at a given dispersion and face subset, returned by the C++ routine \code{_glmbayes_Inv_f3_with_disp}}
#' \item{\code{UB2()}}{Numeric scalar representing the UB2 upper-bound criterion for a given dispersion and face, defined as \eqn{(1/dispersion) * (RSS - rss\_min\_global)}}
#' \item{\code{rss_face_at_disp()}}{Numeric scalar giving the residual sum of squares for a specified face at a given dispersion, computed from cached matrices and the inverse function evaluation}
#' }
#'    
#'
#' @details
#' This function is designed to complement \code{\link{EnvelopeBuild}} for Gaussian models
#' with Normal-Gamma priors. It enables exact sampling of both coefficients and dispersion
#' by constructing a joint envelope that respects posterior curvature in both dimensions.
#'
#' The dispersion anchor point is chosen as the log-scale center of the credible interval,
#' and the Gamma proposal is tilted to match the envelope slope at this point.
#' Theory and narrative: \insertCite{Nygren2006}{glmbayes}; vignettes
#' \code{Chapter-A07}, \code{Chapter-A11}; \insertCite{glmbayesChapterA08,glmbayesIndNormGammaVignette}{glmbayes}.
#' @section Use in accept/reject procedure:
#'
#' The accept/reject sampler relies on a decomposition of the log-posterior into
#' a test statistic and several bounding terms. Each component is constructed so
#' that its sign is controlled, ensuring the validity of the accept/reject step.
#'
#' \describe{
#'
#'   \item{test1 (log-likelihood bound)}{
#'     Placeholder: explain how test1 is formed and why it is non-positive.
#'   }
#'
#'   \item{UB1 (if applicable)}{
#'     Placeholder: describe UB1's role and why it is non-negative.
#'   }
#'
#'   \item{UB2 (residual sum of squares bound)}{
#'     Placeholder: explain how UB2 is constructed from RSS differences and why it is non-negative.
#'   }
#'
#'   \item{UB3A (face-wise quadratic/linear envelope surplus)}{
#'     Placeholder: explain how lg_prob_factor, lmc1, and lmc2 are derived and why UB3A >= 0.
#'   }
#'
#'   \item{UB3B (dispersion-axis envelope surplus)}{
#'     Placeholder: explain how lm_log1, lm_log2, and max_New_LL_UB are used and why UB3B >= 0.
#'   }
#'
#' }
#'
#' Together, these components define
#' \deqn{ test = test1 - UB2 - UB3A - UB3B, }
#' with \eqn{test1 \le 0} and each UB term \eqn{\ge 0}, ensuring the accept/reject
#' procedure is valid and unbiased.
#' @seealso \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeOrchestrator}},
#'   \code{\link{EnvelopeCentering}} (for obtaining \code{RSS_post} and anchored dispersion),
#'   \code{\link{rindepNormalGamma_reg}}, \code{\link{rlmb}};
#'   \code{\link{glmb}}, \code{\link{glmbfamfunc}}.
#' @references
#' \insertAllCited{}
#' @example inst/examples/Ex_EnvelopeDispersionBuild.R
#' @usage EnvelopeDispersionBuild(
#'   Env, Shape, Rate, P, y, x, alpha, n_obs, RSS_post, RSS_ML,
#'   mu, wt, max_disp_perc = 0.99,
#'   disp_lower = NULL, disp_upper = NULL,
#'   verbose = FALSE, use_parallel = TRUE
#' )
#' @export 
#' @order 1


EnvelopeDispersionBuild <- function(Env, Shape, Rate, P, y, x, alpha, n_obs, RSS_post, RSS_ML,
                                    mu, wt, max_disp_perc = 0.99, disp_lower = NULL, disp_upper = NULL,
                                    verbose = FALSE, use_parallel = TRUE) {
  .EnvelopeDispersionBuild_cpp(
        Env, Shape, Rate, P, y, x, alpha, n_obs,
        RSS_post, RSS_ML, mu, wt,
        max_disp_perc, disp_lower, disp_upper,
        verbose, use_parallel)
}


#' Sorts Envelope function for simulation
#'
#' Sorts Enveloping function for simulation.
#' how frequently each component of the resulting grid should be sampled
#' during simulation.
#' @param l1 dimension for model (number of independent variables in X matrix)
#' @param l2 dimension for Envelope (number of components)
#' @param GIndex matrix containing information on how each dimension should be sampled (1 means left tail of a restricted normal, 2 center, 3 right tail, and 4 the entire line)
#' @param G3 A matrix containing the points of tangencies associated with each component of the grid
#' @param cbars A matrix containing the gradients for the negative log-likelihood at each tangency
#' @param logU A matrix containing the log of the cummulative probability associated with each dimension
#' @param logrt A matrix containing the log of the probability associated with the right tail (i.e. that to the right of the lower bound)
#' @param loglt A matrix containing the log of the probability associated with the left tail (i.e., that to the left of the upper bound)
#' @param logP A matrix containing log-probabilities related to the components of the grid
#' @param LLconst A vector containing constant for each component of the grid used during the accept-reject procedure
#' @param PLSD A vector containing the probability of each component in the Grid
#' @param a1 A vector containing the diagonal of the standardized precision matrix
#' @param E_draws Bound on Expected number of candidates per accepted draw
#' @param lg_prob_factor vector of lg_prob_factors used for the Envelope connected to the independent normal gamma prior
#' @param UB2min Vector containing min for UB2 for each component (relevant for EnvelopeDispersionBuild)
#' @return The function returns a list consisting of the following components (the first six of which are matrics with number of rows equal to the number of components in the Grid and columns equal to the number of parameters):
#' \item{GridIndex}{A matrix containing information on how each dimension should be sampled (1 means left tail of a restricted normal, 2 center, 3 right tail, and 4 the entire line)}
#' \item{thetabars}{A matrix containing the points of tangencies associated with each component of the grid}
#' \item{cbars}{A matrix containing the gradients for the negative log-likelihood at each tangency}
#' \item{logU}{A matrix containing the log of the cummulative probability associated with each dimension}
#' \item{logrt}{A matrix containing the log of the probability associated with the right tail (i.e. that to the right of the lower bound)}
#' \item{loglt}{A matrix containing the log of the probability associated with the left tail (i.e., that to the left of the upper bound)}
#' \item{LLconst}{A vector containing constant for each component of the grid used during the accept-reject procedure}
#' \item{logP}{A matrix containing log-probabilities related to the components of the grid}
#' \item{PLSD}{A vector containing the probability of each component in the Grid}
#' \item{E_draws}{A containing a computed theoretical bound on the expected number of draws}
#' \item{sort_ok}{Logical; \code{TRUE} if sort succeeded, \code{FALSE} if memory allocation failed and unsorted envelope was returned (sampler remains valid)}
#' @details This function sorts the envelope in descending order based on the 
#' probability associated with each component in the Grid. Sorting helps 
#' speed up simulation once the envelope is constructed. If memory allocation 
#' fails (e.g. for very large grids), the function returns the unsorted envelope 
#' with \code{sort_ok = FALSE}; the sampler remains valid but may have poorer acceptance.
#'
#' Used after \code{\link{EnvelopeBuild}} and (for Normal--Gamma models)
#' \code{\link{EnvelopeDispersionBuild}}; see \insertCite{Nygren2006,glmbayesChapterA08}{glmbayes}.
#' @seealso \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeOrchestrator}},
#'   \code{\link{EnvelopeDispersionBuild}}, \code{\link{rNormal_reg}}, \code{\link{rglmb}}.
#' @references
#' \insertAllCited{}
#' @example inst/examples/Ex_EnvelopeSort.R
#' @export

EnvelopeSort <- function(l1, l2,
                         GIndex, G3, cbars, logU, logrt, loglt,
                         logP, LLconst, PLSD, a1, E_draws,
                         lg_prob_factor = NULL,
                         UB2min=NULL
                         # ,thetabar_const_base = NULL,
                         # New_LL_Slope=NULL,
                         # shape3_face=NULL
                         ) {   # <-- new optional arg
  # Enforce matrix inputs for logU and LLconst
  if (!is.matrix(logU)) {
    stop("EnvelopeSort: logU must be a matrix, got ", typeof(logU))
  }
  if (!is.matrix(LLconst)) {
    stop("EnvelopeSort: LLconst must be a matrix, got ", typeof(LLconst))
  }

  n_grid <- length(PLSD)
  sort_ok <- TRUE
  tryCatch({
    # Order indices by decreasing PLSD
    ord <- order(PLSD, decreasing = TRUE)
    sel <- ord[seq_len(l2)]  # top l2 rows

    # Reorder inputs once (preserve matrix structure for logU, LLconst)
    GIndex <- GIndex[sel, , drop = FALSE]
    G3     <- G3[sel, , drop = FALSE]
    cbars  <- cbars[sel, , drop = FALSE]
    logU   <- logU[sel, , drop = FALSE]
    logrt  <- logrt[sel, , drop = FALSE]
    loglt  <- loglt[sel, , drop = FALSE]
    logP   <- logP[sel, 1, drop = TRUE]
    LLconst<- LLconst[sel, , drop = FALSE]
    PLSD   <- PLSD[sel]

    if (!is.null(lg_prob_factor)) {
      stopifnot(length(lg_prob_factor) == l2)
      lg_prob_factor <- lg_prob_factor[sel]
    }

    if (!is.null(UB2min)) {
      stopifnot(length(UB2min) == l2)
      UB2min <- UB2min[sel]
    }
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("cannot allocate", msg, ignore.case = TRUE) ||
        grepl("memory", msg, ignore.case = TRUE) ||
        grepl("allocation", msg, ignore.case = TRUE)) {
      warning(
        "EnvelopeSort: memory allocation failed (grid size ", format(n_grid, big.mark = ","),
        "). Caller will use unsorted envelope; sampler remains valid but may have poorer acceptance.",
        call. = FALSE
      )
      return(list(sort_ok = FALSE))
    } else {
      stop(e)
    }
  })
  

  # if (!is.null(thetabar_const_base)) {
  #   stopifnot(length(thetabar_const_base) == l2)
  #   thetabar_const_base <- thetabar_const_base[sel]
  # }
  # 
  # if (!is.null(New_LL_Slope)) {
  #   stopifnot(length(New_LL_Slope) == l2)
  #   New_LL_Slope <- New_LL_Slope[sel]
  # }
  # 
  # if (!is.null(shape3_face)) {
  #   stopifnot(length(shape3_face) == l2)
  #   shape3_face <- shape3_face[sel]
  # }
  # 
  
  
    
  # Build output list
  outlist <- list(
    GridIndex = GIndex,
    thetabars = G3,
    cbars     = if (l1 == 1) matrix(cbars, nrow = l2, ncol = l1) else cbars,
    logU      = logU,
    logrt     = if (l1 == 1) matrix(logrt, nrow = l2, ncol = l1) else logrt,
    loglt     = if (l1 == 1) matrix(loglt, nrow = l2, ncol = l1) else loglt,
    LLconst   = LLconst,
    logP      = logP,
    PLSD      = PLSD,
    a1        = a1,
    E_draws   = E_draws
  )
  if (!is.null(lg_prob_factor)) {
    outlist$lg_prob_factor <- lg_prob_factor
  }
  if (!is.null(UB2min)) {
    outlist$UB2min <- UB2min
  }
  outlist$sort_ok <- sort_ok

  ## Adding to enable face specific dispersion bounds
  # if (!is.null(thetabar_const_base)) {
  #   outlist$thetabar_const_base <- thetabar_const_base
  # }
  # 
  # if (!is.null(New_LL_Slope)) {
  #   outlist$New_LL_Slope <- New_LL_Slope
  # }
  #   
  # if (!is.null(shape3_face)) {
  #   outlist$shape3_face <- shape3_face
  # }
  
  outlist
}


#' The Bayesian Generalized Linear Model Distribution in Standard Form
#'
#' \code{rNormalGLM_std} is used to generate iid samplers from Non-Gaussian Generalized
#' Linear Models in standard form. The function should only be called after standardization 
#' of a Generalized Linear Model.
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x a design matrix of dimension \code{m * p}.
#' @param mu a vector of length \code{p} giving the prior means of the variables in the design matrix.
#' @param P a positive-definite symmetric matrix of dimension \code{p * p} specifying the prior precision
#'  matrix of the variable.
#' @param alpha this can be used to specify an \emph{a priori} known component 
#' to be included in the linear predictor during fitting. This should be 
#' \code{NULL} or a numeric vector of length equal to the number of cases. 
#' One or more offset terms can be included in the formula instead or as well,
#' and if more than one is specified their sum is used. See 
#' \code{\link{model.offset}}.
#' @param wt an optional vector of \sQuote{prior weights} to be used in the fitting process. 
#' Should be NULL or a numeric vector.
#' @param f2 function used to calculate the negative of the log-posterior function
#' @param Envelope an object of type \code{glmbenvelope}.
#' @param family family used for simulation. Used that this is different from the 
#' family used in other functions.
#' @param link link function used for simulation.
#' @param progbar dummy for flagging if a progressbar should be produced during the call
#' @return A list consisting of the following:
#' \item{out}{A matrix with simulated draws from a model in standard form. Each row represents 
#' one draw from the density}
#' \item{draws}{A vector with the number of candidates required before 
#' acceptance for each draw}
#' @details This function uses the information contained in the constructed envelope list in order to sample 
#' from a model in standard form. The simulation proceeds as follows in order to generate each draw in the required
#' number of samples.
#' 
#' 1)  A random number between 0 and 1 is generated and is used together with the information in the PLSD vector 
#' (from the envelope) in order to identify the part of the grid from which a candidate is to be generated.
#'
#' 2) For the part of the grid selected, the dimensions are looped through and a candidate component for each dimension
#' is generated from a restricted normal using information from the Envelope (in particular, the values for logrt, loglt, 
#' and cbars corresponding to that the part of the grid selected and the dimension sampled)
#' 
#' 3) The log-likelihood for the standardized model is evaluated for the generated candidate (note that the 
#' log-likelihood here includes the portion of the prior that was shifted to the log-likelihood)
#' 
#'  4) An additional random number is generated and the log of this random number is compared to a 
#'  log-acceptance rate that is calculated based on the candidate and the LLconst component from the Envelope 
#'  component selected in order to determine if the candidate should be accepted or rejected 
#' 
#' 5) If the candidate was not accepted, the process above is repeated from step 1 until a candidate is accepted
#' 
#' 
#' 
#' @example inst/examples/Ex_rNormalGLM_std.R
#' @export



rNormalGLM_std<-function(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link, progbar = 1L){
  
  return(.rNormalGLM_std_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link, progbar = progbar))
  
}


#' The Bayesian Gaussian Regression with Independent Normal-Gamma Prior in Standard Form
#'
#' \code{rIndepNormalGammaReg_std} generates iid samples from a Bayesian Gaussian
#' regression model with an independent Normal-Gamma prior, in standard form. The
#' function should only be called after standardization and envelope construction
#' (e.g., via \code{\link{EnvelopeOrchestrator}}).
#'
#' @param n Number of draws to generate. If \code{length(n) > 1}, the length is
#'   taken to be the number required.
#' @param y A vector of observations of length \code{m}.
#' @param x A design matrix of dimension \code{m * p}.
#' @param mu A matrix of prior means (typically standardized to zero) of dimension
#'   \code{p * 1}.
#' @param P A positive-definite matrix of dimension \code{p * p} specifying the
#'   prior precision component shifted into the log-likelihood.
#' @param alpha A numeric vector of length \code{m} for the offset-adjusted mean
#'   component. See \code{\link{model.offset}}.
#' @param wt An optional vector of prior weights. Should be \code{NULL} or a
#'   numeric vector.
#' @param f2 Function used to calculate the negative of the log-posterior
#'   (kept for signature parity with \code{rNormalGLM_std}).
#' @param Envelope An envelope object containing \code{PLSD}, \code{loglt},
#'   \code{logrt}, \code{cbars}, and related components from
#'   \code{\link{EnvelopeOrchestrator}}.
#' @param gamma_list A list with \code{shape3}, \code{rate2}, \code{disp_lower},
#'   and \code{disp_upper} from the Gamma posterior for the dispersion.
#' @param UB_list A list with \code{lg_prob_factor}, \code{UB2min}, \code{RSS_Min},
#'   and other bounds from \code{\link{EnvelopeOrchestrator}}.
#' @param family Character vector specifying the family (e.g., \code{"gaussian"}).
#' @param link Character vector specifying the link (e.g., \code{"identity"}).
#' @param progbar Logical. Whether to display a progress bar during simulation.
#' @param verbose Logical. Whether to print diagnostic messages.
#'
#' @return A list with components:
#' \item{beta_out}{A matrix of simulated regression coefficients in standardized
#'   space. Each row is one draw.}
#' \item{disp_out}{A vector of dispersion draws for each sample.}
#' \item{iters_out}{A vector of iteration counts (candidates per acceptance) for
#'   each draw.}
#' \item{weight_out}{A vector of weights (typically all ones).}
#'
#' @details
#' This function uses the envelope and dispersion bounds from
#' \code{\link{EnvelopeOrchestrator}} to sample from the joint posterior of
#' coefficients and dispersion via rejection sampling. It is typically called
#' internally by \code{rindepNormalGamma_reg()}, but may be used directly for
#' custom split workflows (e.g., after constructing the envelope separately).
#'
#' @seealso
#' \code{\link{EnvelopeOrchestrator}} for envelope construction,
#' \code{\link{rNormalGLM_std}} for the non-Gaussian standardized sampler,
#' \code{\link{rindepNormalGamma_reg}} for the full simulation routine.
#'
#' @example inst/examples/Ex_rIndepNormalGammaReg_std.R
#'
#' @export
rIndepNormalGammaReg_std <- function(n, y, x, mu, P, alpha, wt, f2, Envelope,
                                     gamma_list, UB_list, family, link,
                                     progbar = TRUE, verbose = FALSE) {
  .rIndepNormalGammaReg_std_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope,
                                gamma_list, UB_list, family, link,
                                progbar, verbose)
}


# Helpers --------------------------------------------------------------------


#' @noRd

dpois2<-function(x,lambda,log=TRUE){
  
  test=max(abs(round(x)-x))
  
  if(test>0){
    warning("Non-Integer Values to Poisson Density - Switching to Gamma Function to Evaluate Factorial")
    return(-lambda+x*log(lambda)-log(gamma(x+1)))
    
  } 
  
  return(dpois(x,lambda,log=TRUE))
}


