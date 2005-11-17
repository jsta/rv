\name{rvnorm}
\title{Generate Random Variables from a Gaussian (Normal) Sampling Model}
\alias{rvnorm}
\usage{
  rvnorm(n=1, mean=0, sd=1, var=NULL, precision)
}
\description{
  Generates a random vector from a Gaussian sampling model.
}
\arguments{
  \item{n}{integer: number of variables to generate}
  \item{mean}{mean, may be a rv}
  \item{sd}{standard deviation; scalar or vector (constant or rv, not matrix)}
  \item{var}{variance, can be given instead of sd. Scalar, vector, or matrix.}
  \item{precision}{inverse variance or variance matrix, may be given instead of sd or var}
}
\note{
If any of the arguments are random, 
the resulting simulations may have non-normal marginal distributions;
for example, if an inverse-chi-squared scalar rv \code{var} 
and zero \code{mean} is given, the resulting rv will have a t-distribution.
}
\value{
  A random vector (rv object) of length \code{n}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvnorm(mean=1:10, sd=1:10) # A vector of length 10.
  Sigma <- diag(1:10)
  y <- rvnorm(mean=1:10, var=Sigma)
}

\keyword{classes}

