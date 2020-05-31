\name{summary.glmb}
\alias{summary.glmb}
\alias{print.summary.glmb}
\title{
Summarizes the output from the \code{\link{glmb}} function
}
\description{
The \code{summary.glmb} function summarizes the output from the \code{\link{glmb}} function. It takes an object of class \code{glmb} as an input. The output to a large extent mirrors the output from the \code{\link{summary.glm}} function. This is particularly true for the print output from the function (i.e. the output of the function \code{\link{print.summary.glmb}}). 
}
\usage{
\method{summary}{glmb}(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{rglmb} for which a summary is desired
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
The \code{summary.glmb} function summarizes the output from the \code{\link{glmb}} function. Key output includes mean residuals, information related to the prior, mean coefficients with associated stats, percentiles for the coefficients, as well as the effective number of parameters and the DIC statistic.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A list with the following elements:
 \item{n }{number of draws generated}
  \item{residuals}{vector of mean deviance residuals}
  \item{coefficients1}{Matrix with the prior mean and maximum likelihood coefficients with associated standard deviations}
  \item{coefficients}{Matrix with columns for the posterior mode, posterior mean, posterior standard deviation, monte carlo error, and tail probabilities (posterior probability of observing a value for the coefficient as extreme as the prior mean)}
  \item{Percentiles}{Matrix with estimated percentiles associated with the posterior density}
  \item{pD}{Estimated effective number of parameters}
  \item{deviance}{Vector with draws for the deviance}
  \item{DIC}{Estimated DIC statistic}
  \item{iters}{Average number of candidates per generated draws}
%% ...
}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
Spiegelhalter, et al,
}
\author{
Kjell Nygren
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glmb}}, \code{\link{rglmb}}, \code{\link{summary.rglmb}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, ...) 
{
    mres <- colMeans(residuals(object))
    l1 <- length(object$coef.means)
    n <- length(object$coefficients[, 1])
    percentiles <- matrix(0, nrow = l1, ncol = 7)
    se <- sqrt(diag(var(object$coefficients)))
    mc <- se/object$n
    mc <- se/n
    priorrank <- matrix(0, nrow = l1, ncol = 1)
    pval1 <- matrix(0, nrow = l1, ncol = 1)
    pval2 <- matrix(0, nrow = l1, ncol = 1)
    for (i in 1:l1) {
        percentiles[i, ] <- quantile(object$coefficients[, i], 
            probs = c(0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99))
        test <- append(object$coefficients[, i], object$Prior$mean[i])
        test2 <- rank(test)
        priorrank[i, 1] <- test2[n + 1]
        priorrank[i, 1] <- test2[n + 1]
        pval1[i, 1] <- priorrank[i, 1]/(n + 1)
        pval2[i, 1] <- min(pval1[i, 1], 1 - pval1[i, 1])
    }
    glmsummary <- summary(object$glm)
    se1 <- sqrt(diag(glmsummary$cov.scaled))
    Tab1 <- cbind(`Prior Mean` = object$Prior$mean, 
    Prior.sd = as.numeric(sqrt(diag(object$Prior$Variance))), 
        `Max Like.` = as.numeric(object$glm$coefficients), Like.sd = se1)
    TAB <- cbind(Post.Mode = as.numeric(object$coef.mode), 
    Post.Mean = object$coef.means, 
        Post.Sd = se, `MC Error` = as.numeric(mc), `Pr(tail)` = as.numeric(pval2))
    TAB2 <- cbind(`1.0\%` = percentiles[, 1], `2.5\%` = percentiles[, 
        2], `5.0\%` = percentiles[, 3], Median = as.numeric(percentiles[, 
        4]), `95.0\%` = percentiles[, 5], `97.5\%` = as.numeric(percentiles[, 
        6]), `99.0\%` = as.numeric(percentiles[, 7]))
    rownames(TAB2) <- rownames(TAB)
    res <- list(call = object$call, n = n, residuals = mres, 
        coefficients1 = Tab1, coefficients = TAB, Percentiles = TAB2, 
        pD = object$pD, deviance = object$deviance, DIC = object$DIC, 
        iters = mean(object$iters))
    class(res) <- "summary.glmb"
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
