\name{glmbenvelope}
\alias{glmbenvelope}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructs Enveloping function for simulation
}
\description{
Constructs Enveloping function for simulation
}
\usage{
glmbenvelope(bStar, A, f2, f3,f5,f6, y, x, mu, P, alpha, wt = 1, Gridtype = 3, n = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bStar}{Point at which envelope should be centered (typically posterior mode)
}
  \item{A}{Diagonal precision matrix for log-likelihood function associated with model in standard form
}
  \item{f2}{Negative log-posterior function
}
  \item{f3}{Gradient function for Negative log-posterior function
  
  }
  \item{f5}{Vectorized version of Negative log-posterior function
}
  \item{f6}{Vectorized version of Gradient function for Negative log-posterior function}
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
  \item{alpha}{offset vector 
}
  \item{wt}{
a vector of weights
}
  \item{Gridtype}{
an optional argument specifying the method used to determine number of likelihood subgradient 
densities used to construct the enveloping function.
}
  \item{n}{
number of draws to generate from posterior density. Used here to help determine the size of the grid (Gridtype 1 and 2 only)
}
}
\details{
This function identfies tangency points and intervals for the enveloping function, evaluates the log-likelihood function and gradients for the log-likelihood at each point, and calculates additional constants.
}
\value{
The function returns a list consisting of the following components (the first six of which are matrics with number of rows equal to the number of components in the Grid and columns equal to the number of parameters):
 \item{GridIndex}{A matrix containing information on how each dimension should be sampled (1 means left tail of a restricted normal, 2 center, 3 right tail, and 4 the entire line) }
  \item{thetabars }{A matrix containing the points of tangencies associated with each component of the grid}
  \item{cbars}{A matrix containing the gradients for the negative log-likelihood at each tangency}
  \item{logU}{A matrix containing the log of the cummulative probability associated with each dimension}
  \item{logrt}{A matrix containing the log of the probability associated with the right tail (i.e. that to the right of the lower bound)}
  \item{loglt}{A matrix containing the log of the probability associated with the left tail (i.e., that to the left of the upper bound)}
  \item{LLconst}{A vector containing constant for each component of the grid used during the accept-reject procedure}  
  \item{logP}{A matrix containing log-probabilities related to the components of the grid }  
  \item{PLSD}{A vector containing the probability of each component in the Grid}  


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
\code{\link{rglmb}}, \code{\link{expand.grid}}, \code{\link{optgrid}}, \code{\link{Set_Grid}}, \code{\link{setlogP}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (bStar, A, f2, f3, y, x, mu, P, alpha, wt = 1, Gridtype = 3, 
    n = 1) 
{
    a1 <- diag(A)
    l1 <- length(a1)
    omega <- as.matrix((sqrt(2) - exp(-1.20491 - 0.7321 * sqrt(0.5 + 
        a1)))/sqrt(1 + a1))
    G1 <- as.matrix(c(1, 1, 1)) \%*\% t(bStar) + as.matrix(c(-1, 
        0, 1)) \%*\% t(omega)
    Lint <- as.matrix(c(1, 1)) \%*\% t(bStar) + as.matrix(c(-0.5, 
        0.5)) \%*\% t(omega)
    G2 <- list()
    GIndex1 <- list()
    length(G2) <- l1
    length(GIndex1) <- l1
    if (Gridtype == 2) {
        gridindex <- optgrid(a1, n = n)
    }
    for (i in 1:l1) {
        if (Gridtype == 1) {
            if ((1 + a1[i]) <= (2/sqrt(pi))) {
                G2[i] <- list(G1[2, i])
                GIndex1[i] <- 4
            }
            if ((1 + a1[i]) > (2/sqrt(pi))) {
                G2[i] <- list(G1[, i])
                GIndex1[i] <- list(c(1, 2, 3))
            }
        }
        if (Gridtype == 2) {
            if (gridindex[i] == 1) {
                G2[i] <- list(G1[2, i])
                GIndex1[i] <- 4
            }
            if (gridindex[i] == 3) {
                G2[i] <- list(G1[, i])
                GIndex1[i] <- list(c(1, 2, 3))
            }
        }
        if (Gridtype == 3) {
            G2[i] <- list(G1[, i])
            GIndex1[i] <- list(c(1, 2, 3))
        }
        if (Gridtype == 4) {
            G2[i] <- list(G1[2, i])
            GIndex1[i] <- 4
        }
    }
    G3 <- expand.grid(G2)
    GIndex <- expand.grid(GIndex1)
    l2 <- length(GIndex[, 1])
    cbars <- matrix(0, nrow = l2, ncol = l1)
    NegLL <- matrix(0, nrow = l2, ncol = 1)
    Up <- matrix(0, nrow = l2, ncol = l1)
    Down <- matrix(0, nrow = l2, ncol = l1)
    logP <- matrix(0, nrow = l2, ncol = 2)
    logU <- matrix(0, nrow = l2, ncol = l1)
    loglt <- matrix(0, nrow = l2, ncol = l1)
    logrt <- matrix(0, nrow = l2, ncol = l1)
    PLSD <- matrix(0, nrow = l2, ncol = 1)
    LLconst <- matrix(0, nrow = l2, ncol = 1)
    for (j in 1:l2) {
        cbars[j, 1:l1] <- t(f3(b = t(as.matrix(G3[j, 1:l1])), 
            y = y, x = x, mu = mu, P = P, alpha = alpha, wt = wt))
        NegLL[j] <- f2(b = t(as.matrix(G3[j, 1:l1])), y = y, 
            x = x, mu = mu, P = P, alpha = alpha, wt = wt)
        LLconst[j] <- -t(as.matrix(cbars[j, 1:l1])) \%*\% t(as.matrix(G3[j, 
            1:l1])) + NegLL[j]
        for (k in 1:l1) {
            if (GIndex[j, k] == 1) {
                Down[j, k] <- -Inf
                Up[j, k] <- Lint[1, k] + cbars[j, k]
                logrt[j, k] <- 0
                loglt[j, k] <- pnorm(Up[j, k], log.p = TRUE)
                logU[j, k] <- pnorm(Up[j, k], log.p = TRUE)
                if (logU[j, k] == -Inf) {
                  stop("Function pnorm returned -Inf")
                }
                logP[j, 1] <- logP[j, 1] + loglt[j, k]
            }
            if (GIndex[j, k] == 2) {
                Down[j, k] <- Lint[1, k] + cbars[j, k]
                Up[j, k] <- Lint[2, k] + cbars[j, k]
                loglt[j, k] <- pnorm(Up[j, k], log.p = TRUE)
                logrt[j, k] <- pnorm(Down[j, k], lower.tail = FALSE, 
                  log.p = TRUE)
                logU[j, k] <- ctpnorm(a = Down[j, k], b = Up[j, 
                  k], log.p = TRUE)
                logP[j, 1] <- logP[j, 1] + logU[j, k]
            }
            if (GIndex[j, k] == 3) {
                Down[j, k] <- Lint[2, k] + cbars[j, k]
                Up[j, k] <- Inf
                logrt[j, k] <- pnorm(Down[j, k], lower.tail = FALSE, 
                  log.p = TRUE)
                loglt[j, k] <- 0
                logU[j, k] <- pnorm(Down[j, k], lower.tail = FALSE, 
                  log.p = TRUE)
                if (logU[j, k] == -Inf) {
                  stop("Function pnorm returned -Inf")
                }
                logP[j, 1] <- logP[j, 1] + logrt[j, k]
            }
            if (GIndex[j, k] == 4) {
                Down[j, k] <- -Inf
                Up[j, k] <- Inf
                loglt[j, k] <- 0
                logrt[j, k] <- 0
                logU[j, k] <- 0
                logP[j, 1] <- logP[j, 1] + logU[j, k]
            }
        }
        logP[j, 2] <- logP[j, 1] - NegLL[j] + 0.5 * t(as.matrix(cbars[j, 
            1:l1])) \%*\% as.matrix(cbars[j, 1:l1]) + as.matrix(G3[j, 
            1:l1]) \%*\% as.matrix(cbars[j, 1:l1])
    }
    maxlogP <- max(logP[, 2])
    PLSD <- exp(logP[, 2] - maxlogP)
    sumP <- sum(PLSD)
    PLSD <- PLSD/sumP
    Envelope <- data.frame(GIndex = GIndex, G3 = G3, cbars = cbars, 
        logU = logU, logrt = logrt, loglt = loglt, logP = logP[, 
            1], LLconst = LLconst, PLSD = PLSD)
    Envelope = Envelope[order(-Envelope$PLSD), ]
    list(GridIndex = Envelope[1:l2, 1:l1], thetabars = Envelope[1:l2, 
        (l1 + 1):(2 * l1)], cbars = Envelope[1:l2, (2 * l1 + 
        1):(3 * l1)], logU = Envelope[1:l2, (3 * l1 + 1):(4 * 
        l1)], logrt = Envelope[1:l2, (4 * l1 + 1):(5 * l1)], 
        loglt = Envelope[1:l2, (5 * l1 + 1):(6 * l1)], LLconst = Envelope$LLconst, 
        logP = Envelope$logP, PLSD = Envelope$PLSD, a1 = a1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
