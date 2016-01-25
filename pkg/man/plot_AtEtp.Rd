\name{plot_AtEtp}
\alias{plot_AtEtEp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot variance curves for the AE(t)-p model  
}
\description{
Plot variance curves of the A and E components with respect to age modelled by penalized B-splines for the AE(t)-p model .
}
\usage{
plot_AtEtp(AtEtp_mcmc, data_m, data_d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AtEtp_mcmc}{
%%     ~~Describe \code{pheno_m} here~~
}
  \item{data_m}{
%%     ~~Describe \code{pheno_m} here~~
}
  \item{data_d}{
%%     ~~Describe \code{pheno_d} here~~
}
 
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

data(data_ace)
#result <- AtEtp(data_ace$mz, data_ace$dz)
#result_mc <- AtEtp_mcmc(result)
#plot_AtEtp(result_mc, data_ace$mz, data_ace$dz)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
