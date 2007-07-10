\name{rvmultinom}
\title{Generate Random Variables from a Multinomial Sampling Model}
\alias{rvmultinom}
\usage{
  rvmultinom(n=1, size=1, prob)
}
\description{
  Generates a random vector from a multinomial sampling model.  
}
\arguments{
  \item{n}{integer, number of random variables to generate}
  \item{size}{integer or integer-valued rv: the number of trials (size of each sample)}
  \item{prob}{vector (of length at least 3) prior probabilities of successes of each trial (may be constant or an rv object)}
}
\details{
  The length of \code{prob} determines the number of bins.

  The vector \code{prob} will be normalized to have sum 1.


  If \code{length(prob)} is two, \code{rvbinom} is called instead.

  NOTE. Case of random \code{n} or \code{size} or \code{prob} --- not yet 
  optimized for speed.
}
\value{
  A random array of dimensions \code{length(prob)} times \code{n}.
}
\references{ 
Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations 
Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  y <- rvmultinom(n=3, size=1, prob=c(0.20, 0.30, 0.50))
}

\keyword{classes}

