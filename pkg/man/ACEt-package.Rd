\name{ACEt-package}
\alias{ACEt-package}
\alias{ACEt}
\docType{package}
\title{
Estimating age modification effect on genetic and environmental variance components in twin models
}
\description{
Novel twin models that are able to estimate the dynamic behaviour of the variance components in the classical twin models with respect to age. 
}
\details{
\tabular{ll}{
Package: \tab ACEt\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2016-01-11\cr
License: \tab GPL \cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Liang He

Maintainer: Liang He <liang.he@duke.edu>

}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{

data(data_ace)

AtCtE(data_ace$mz, data_ace$dz)
AtCtEp(data_ace$mz, data_ace$dz)

}
