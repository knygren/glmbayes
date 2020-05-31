\name{summary.rglmb}
\alias{summary.rglmb}
\alias{print.summary.rglmb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes the output from the \code{\link{rglmb}} function
}
\description{
The \code{summary.rglmb} function summarizes the output from the \code{\link{rglmb}} function. It takes an object of class \code{rglmb} as an input. The output to a large extent mirrors the output from the \code{\link{summary.glm}} function. This is particularly true for the print output from the function (i.e. the output of the function \code{\link{print.summary.rglmb}}).
}
\usage{
\method{summary}{rglmb}(object, ...)
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
The \code{summary.rglmb} function summarizes the output from the \code{\link{rglmb}} function. Key output includes mean residuals, information related to the prior, mean coefficients with associated stats, percentiles for the coefficients, as well as the effective number of parameters and the DIC statistic.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{n}{number of draws generated}
  \item{coefficients1}{Matrix with the prior mean and approximate weight for the prior relative to the data}
  \item{coefficients}{Matrix with columns for the posterior mode, posterior mean, posterior standard deviation, monte carlo error, and tail probabilities (posterior probability of observing a value for the coefficient as extreme as the prior mean)}
  \item{Percentiles}{Matrix with estimated percentiles associated with the posterior density}
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
\code{\link{rglmb}}, \code{\link{glmb}}, \code{\link{summary.glmb}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, ...) 
{
    n <- length(object$coefficients[, 1])
    l1 <- length(object$PostMode)
    percentiles <- matrix(0, nrow = l1, ncol = 7)
    se <- sqrt(diag(var(object$coefficients)))
    mc <- se/n
    Priorwt <- (se/sqrt(diag(solve(object$Prior$Precision))))^2
    priorrank <- matrix(0, nrow = l1, ncol = 1)
    pval1 <- matrix(0, nrow = l1, ncol = 1)
    pval2 <- matrix(0, nrow = l1, ncol = 1)
    for (i in 1:l1) {
        percentiles[i, ] <- quantile(object$coefficients[, i], 
            probs = c(0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99))
        test <- append(object$coefficients[, i], object$Prior$mean[i])
        test2 <- rank(test)
        priorrank[i, 1] <- test2[n + 1]
        pval1[i, 1] <- priorrank[i, 1]/(n + 1)
        pval2[i, 1] <- min(pval1[i, 1], 1 - pval1[i, 1])
    }
    Tab1 <- cbind(`Prior Mean` = object$Prior$mean, Prior.sd = 
    as.numeric(sqrt(diag(solve(object$Prior$Precision)))), 
        Approx.Prior.wt = Priorwt)
    TAB <- cbind(Post.Mode = as.numeric(object$PostMode), Post.Mean = colMeans(coef(object)), 
        Post.Sd = se, `MC Error` = as.numeric(mc), `Pr(tail)` = as.numeric(pval2))
    TAB2 <- cbind(`1.0\%` = percentiles[, 1], `2.5\%` = percentiles[, 
        2], `5.0\%` = percentiles[, 3], Median = as.numeric(percentiles[, 
        4]), `95.0\%` = percentiles[, 5], `97.5\%` = as.numeric(percentiles[, 
        6]), `99.0\%` = as.numeric(percentiles[, 7]))
    rownames(Tab1) <- rownames(TAB)
    rownames(TAB2) <- rownames(TAB)
    res <- list(call = object$call, n = object$n, coefficients1 = Tab1, 
        coefficients = TAB, Percentiles = TAB2, n = object$n)
    class(res) <- "summary.rglmb"
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
