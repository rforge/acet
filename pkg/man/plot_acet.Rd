\name{plot_acet}
\alias{plot_acet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot variance curves for various ACE(t) or ACE(t)-p models  
}
\description{
Plot variance curves of the A, C and E components with respect to age modelled by B-splines or P-splines.
}
\usage{
plot_acet(acet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{acet}{
An object obtained from the B-splines or P-splines functions. For the P-splines functions, an object from the MCMC method must be used.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(data_ace)
result <- AtCtEp(data_ace$mz, data_ace$dz)
result_mc <- acetp_mcmc(result, iter_num=5000)
plot_acet(result_mc)

result <- AtCtEt(data_ace$mz, data_ace$dz, model=c('d','d','c'))
plot_acet(result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
