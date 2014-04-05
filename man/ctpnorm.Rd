\name{ctpnorm}
\alias{ctpnorm}
\alias{ctrnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Finds the probability of the center of a normal density (probability of area between lower and upper bound)}
\description{
Finds the probability of the center of a normal density (probability of area between lower and upper bound)}
\usage{
ctpnorm(a = -Inf, b = Inf, mu = 0, sigma = 1, log.p = TRUE, Diff = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
lower bound}
  \item{b}{
upper bound}
  \item{mu}{
mean parameter}
  \item{sigma}{
standard deviation
}
  \item{log.p}{
Logical argument. If True, the log probability is provided}
  \item{Diff}{
If True the second parameter is the difference between the lower and upper bound
}
}
\details{
Finds the probability of the center of a normal density (probability of area between lower and upper bound)
}
\value{
A vector of length equal to length of a}
\references{
Nygren, K.N, and Nygren, L.M, Likelihood Subgradient Densities.  Journal of the American Statistical Association.
}
\author{
Kjell Nygren
}


\seealso{
\code{\link{glmbsim_cpp}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (a = -Inf, b = Inf, mu = 0, sigma = 1, log.p = TRUE, 
    Diff = FALSE) 
{
    if (is.vector(a) == FALSE || is.vector(a) == FALSE) {
        stop("Arguments a and b must be vectors")
    }
    if (is.numeric(a) == FALSE || is.numeric(b) == FALSE || is.numeric(mu) == 
        FALSE || is.numeric(sigma) == FALSE) {
        stop("Non-numeric argument to mathematical function")
    }
    if (is.logical(log.p) == FALSE || is.logical(Diff) == FALSE) {
        stop("Arguments log.P and Diff must be logical")
    }
    if (length(a) != length(b)) {
        stop("a and b must be equal length vectors")
    }
    if (length(mu) > 1 || length(sigma) > 1) {
        stop("mu and sigma must have length 1")
    }
    if (sigma <= 0) {
        stop("sigma must be positive")
    }
    l1 <- length(a)
    output <- 0 * c(1:length(a))
    for (i in 1:length(a)) {
        if (Diff == FALSE) {
            b2 <- b[i]
            if (length(a) > 1) {
                Diff2 <- b[i] - a[i]
            }
            if (length(a) == 1) {
                Diff2 <- b[i] - a
            }
            if (Diff2 < 0) {
                stop("Argument b must be greater than argument a")
            }
        }
        if (Diff == TRUE) {
            b2 <- a[i] + b[i]
            Diff2 <- b[i]
            if (Diff2 <= 0) {
                stop("Argument b Must be positive")
            }
        }
        g1 <- mu - a
        g2 <- b2 - mu
        if (log.p == FALSE) {
            if (g1 >= g2) {
                output[i] <- pnorm(q = b2, mean = mu, sd = sigma) * 
                  (1 - exp(pnorm(q = a, mean = mu, sd = sigma, 
                    log.p = TRUE) - pnorm(q = b2, mean = mu, 
                    sd = sigma, log.p = TRUE)))
            }
            if (g1 < g2) {
                output[i] <- pnorm(q = a, mean = mu, sd = sigma, 
                  lower.tail = FALSE) * (1 - exp(pnorm(q = b2, 
                  mean = mu, sd = sigma, log.p = TRUE, lower.tail = FALSE) - 
                  pnorm(q = a, mean = mu, sd = sigma, log.p = TRUE, 
                    lower.tail = FALSE)))
            }
            if (g1 == Inf && g2 == Inf) {
                output[i] <- 1
            }
            if (output[i] == 0) {
                output[i] <- Diff2 * dnorm(x = 0.5 * (b2 + a), 
                  mean = mu, sd = sigma)
                if (output[i] == 0) {
                  stop("Function ctpnorm evaluates to 0")
                }
            }
        }
        if (log.p == TRUE) {
            output[i] <- NA
            output[i] <- log(pnorm(q = b, mean = mu, sd = sigma) - 
                pnorm(q = a, mean = mu, sd = sigma))
        }
    }
    output
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
