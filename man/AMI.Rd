\name{AMI}
\alias{AMI}
\docType{data}
\title{
amitriptyline overdose data}
\description{
Data with information on 17 overdoses of the drug amitriptyline

}
\usage{data(AMI)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{TOT}}{total TCAD plasma level}
    \item{\code{AMI}}{amount of amitriptyline present in the TCAD plasma level}
    \item{\code{GEN}}{gender (male = 0, female = 1)}
    \item{\code{AMT}}{amount of drug taken at time of overdose}
    \item{\code{PR}}{PR wave measurement}
    \item{\code{DIAP}}{diastolic blood pressure}
    \item{\code{QRS}}{QRS wave measurement}
  }
}
\details{
Details to be added...
}
\references{Data downloaded from https://data.library.virginia.edu/getting-started-with-multivariate-multiple-regression/. \cr

Fox, J and Weisberg, S (2011). An {R} Companion to Applied Regression, Second Edition. Thousand Oaks CA: Sage. URL: http://socserv.socsci.mcmaster.ca/jfox/Books/Companion
Johnson, R and Wichern, D (2007). Applied Multivariate Statistical Analysis, Sixth Edition. Prentice-Hall.
Rudorfer, MV “Cardiovascular Changes and Plasma Drug Levels after Amitriptyline Overdose.” Journal of Toxicology-Clinical Toxicology, 19 (1982), 67-71.

}
\examples{
data(menarche2)
}
\keyword{datasets, Bayesian Binomial Regression}
