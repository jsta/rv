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
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvpois(lambda=10)  # A Poisson rv with mean 10
  lbd <- rvchisq(1,1)     # Some positive rv
  y <- rvpois(lambda=lbd) # Not a Poisson rv, although each simulation is a draw from Poisson.
}

\keyword{classes}


