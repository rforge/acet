\name{AtCEt}
\alias{AtCEt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The ACE(t) model  
}
\description{
The ACE(t) model with the A and E variance components as functions with respect to age modelled by B-splines.
}
\usage{
AtCEt(data_m, data_d, knot_a = 5, knot_e = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_m}{
An \eqn{N_m} x 3 data matrix for MZ twins. \eqn{N_m} is the number of MZ twin pairs. The first two columns are traits and the third column is age (or other covariates). 
}
  \item{data_d}{
An \eqn{N_d} x 3 data matrix for DZ twins. \eqn{N_d} is the number of DZ twin pairs. The first two columns are traits and the third column is age (or other covariates). 
}
  \item{knot_a}{
The number of knots of the B-spline for the A component. The default value is 5.
}
  \item{knot_e}{
The number of knots of the B-spline for the E component. The default value is 5.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{n_beta_a }{The number of spline coefficients for the A component.}
  \item{n_beta_e }{The number of spline coefficients for the E component.}
  \item{var }{The estimated variance of the C component.}
  \item{beta_a }{The estimated spline coefficients of the A component.}
  \item{beta_e }{The estimated spline coefficients of the E component.}
  \item{hessian }{The approximated Hessian matrix from the minus log-likelihood.}
  \item{con }{The indicator of convergence of the optimization algorithm.}
  \item{lik }{The minus log-likelihood.}
  \item{knots_a }{The values of the knots for the A component.}
  \item{knots_e }{The values of the knots for the E component.}
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
data(data_ace)

result <- AtCEt(data_ace$mz, data_ace$dz)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
