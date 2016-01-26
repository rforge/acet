\name{AtCtEp_mcmc}
\alias{AtCtEp_mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute CIs for the ACE(t)-p model
}
\description{
Compute CIs for the ACE(t)-p model using the MCMC methods and the bootstrap method
}
\usage{
AtCtEp_mcmc(AtCtEp, iter_num = 10000, sd = 0.1, burnin = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AtCtEp}{
An object of AtCtEp obtained from the function 'AtCtEp'.
}
  \item{iter_num}{
The number of MCMC iteration
}
  \item{sd}{
The standard error of the normal proposal distribution. The default value is 0.1.
}
  \item{burnin}{
The number of burn-in, which must be smaller than \eqn{iter_num}
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

data(data_ace)

result <- AtCtEp(data_ace$mz, data_ace$dz)
result_mc <- AtCtEp_mcmc(result)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
