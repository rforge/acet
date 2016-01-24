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
require(MASS)

num_m <- 10000
num_d <- 10000

T_m <- rep(1:50, each=100)
T_d <- rep(1:50, each=100)

var_c <- 1
var_a <- 1
var <- 1

beta_c <- -0.20
beta_a <- 2

pheno_m <- matrix(0, num_m/2,2)	
		
for(i in 1:(num_m/2))
{
	 var_ct <- var_c*(T_m[i]^beta_c)
	var_at <- var_a*(3-10*((T_m[i]/50-0.5)^beta_a))
	sigma <- matrix(c(var_at+var_ct+var,var_at+var_ct,var_at+var_ct,var_at+var_ct+var),2,2)

	pheno_m[i,1:2] <- mvrnorm(1, rep(0,2), sigma)
}

pheno_d <- matrix(0, num_d/2,2)	
		
for(i in 1:(num_d/2))
{
	 var_ct <- var_c*(T_d[i]^beta_c)
	var_at <- var_a*(3-10*((T_d[i]/50-0.5)^beta_a))
	sigma <- matrix(c(var_at+var_ct+var,0.5*var_at+var_ct,0.5*var_at+var_ct,var_at+var_ct+var),2,2)

	pheno_d[i,1:2] <- mvrnorm(1, rep(0,2), sigma)
}

#AtCtE(cbind(pheno_m,T_m), cbind(pheno_d,T_d))
#AtCtEp(cbind(pheno_m,T_m), cbind(pheno_d,T_d))

}
