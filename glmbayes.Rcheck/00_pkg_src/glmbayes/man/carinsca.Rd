\name{carinsca}
\alias{carinsca}
\docType{data}
\title{
Canadian Automobile Insurance Claims for 1957-1958
}
\description{
The data give the Canadian automobile insurance experience for policy years 1956 and 1957 as of June 30, 1959. The data includes virtually every insurance company operating in Canada and was collated by the Statistical Agency (Canadian Underwriters' Association - Statistical Department) acting under instructions from the Superintendent of Insurance. The data given here is for private passenger automobile liability for non-farmers for all of Canada excluding Saskatchewan. 
}
\usage{data(carinsca)}
\format{
  A data frame with 20 observations on the following 6 variables.
  \describe{
    \item{\code{Merit}}{Merit Rating:\cr
3 - licensed and accident free 3 or more years\cr
2 - licensed and accident free 2 years\cr
1 - licensed and accident free 1 year\cr
0 - all others}
    \item{\code{Class}}{1 - pleasure, no male operator under 25\cr
2 - pleasure, non-principal male operator under 25\cr
3 - business use\cr
4 - unmarried owner or principal operator under 25\cr
5 - married owner or principal operator under 25}
    \item{\code{Insured}}{Earned car years}
    \item{\code{Premium}}{Earned premium in 1000's\cr
(adjusted to what the premium would have been had all cars been written at 01 rates)}
    \item{\code{Claims}}{Number of claims}
    \item{\code{Cost}}{Total cost of the claim in 1000's of dollars}
  }
}
\details{
One could apply Poisson regression to the number of claims and gamma regression to the cost per claim. 
}
\source{
Bailey, R. A., and Simon, LeRoy J. (1960). Two studies in automobile insurance ratemaking. ASTIN Bulletin, 192-217. 
}
\references{
Data downloaded from http://www.statsci.org/data/general/carinsca.html. \cr
That site also contains classical Poisson and Gamma regression examples.
}
\examples{
data(carinsca)
## maybe str(carinsca) ; plot(carinsca) ...
}
\keyword{datasets, Bayesian Poisson Regression, Bayesian Gamma Regression}
