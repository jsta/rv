\name{rvpois}
\title{Generate Random Vectors from a Poisson Sampling Model}
\alias{rvpois}
\usage{
  rvpois(n=1, lambda)
}
\description{
  Generates random variables from a Poisson sampling model.
}
\arguments{
  \item{n}{integer: number of variables to generate}
  \item{lambda}{a vector of (positive) mean parameters; (may be random)}
}
\note{
If any of the arguments are random, 
the resulting simulations may have non-Poisson marginal distributions.
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  x <- rvpois(lambda=10)  # A Poisson rv with mean 10
  lbd <- rvchisq(1,1)     # Some positive rv
  y <- rvpois(lambda=lbd) # Not a Poisson rv, although each simulation is a draw from Poisson.
}

\keyword{classes}


