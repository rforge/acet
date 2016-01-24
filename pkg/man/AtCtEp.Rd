\name{AtCtEp}
\alias{AtCtEp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The ACE(t)-p model 
}
\description{
The ACE(t)-p model with the A and C variance components as functions with respect to age modelled by P-splines.
}
\usage{
AtCtEp(data_m, data_d, knot_a = 8, knot_c = 8, eps = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_m}{
%%     ~~Describe \code{ph_m} here~~
}
  \item{data_d}{
%%     ~~Describe \code{ph_d} here~~
}
  \item{knot_a}{
%%     ~~Describe \code{knot_a} here~~
}
  \item{knot_c}{
%%     ~~Describe \code{knot_c} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
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

data_m <- cbind(pheno_m, T_m)
data_d <- cbind(pheno_d, T_d)

#result <- AtCtEp(data_m, data_d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
