\name{plot_AtCtEp}
\alias{plot_AtCtEp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot variance curves for the ACE(t)-p model  
}
\description{
Plot variance curves of the A and C components with respect to age modelled by penalized B-splines for the ACE(t)-p model .
}
\usage{
plot_AtCtEp(AtCtEp_mcmc, data_m, data_d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AtCtEp_mcmc}{
An object of AtCtEp_mcmc obtained from the function 'AtCtEp_mcmc'.
}
  \item{data_m}{
The data for the MZ twins used in the function 'AtCtEp'. 
}
  \item{data_d}{
The data for the DZ twins used in the function 'AtCtEp'.
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
Liang He
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
result <- AtCtEp(data_ace$mz, data_ace$dz)
result_mc <- AtCtEp_mcmc(result)
plot_AtCtEp(result_mc, data_ace$mz, data_ace$dz)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
