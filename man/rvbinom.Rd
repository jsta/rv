\name{rvbinom}
\title{Generate Random Variables from a Binomial Sampling Model}
\alias{rvbinom}
\usage{
  rvbinom(n=1, size, prob)
}
\description{
  Generates a random vector from a binomial sampling model.
}
\arguments{
  \item{n}{integer, number of random variables to generate}
  \item{size}{integer or integer-valued rv: the number of trials (size of each sample)}
  \item{prob}{prior probability of success of each trial (may be constant or an rv object)}
}
\details{
  \code{rvbinom} generates a random vector with given
  length, 
  the distribution for size and the distribution for the probability of success.
}
\value{

}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  s <- 1+rvpois(1,lambda=3)        # A prior distribution on the 'size' parameter.
  rvbinom(1, size=s, prob=0.5)     # The 'size' is random.
  p <- rvbinom(1, 10, prob=0.5)/10 # Prior probability of success.
  rvbinom(1, size=10, prob=p)      # Now the probability is random.
  rvbinom(1, size=s, prob=p)       # Both the size and the probability are random.
}

\keyword{classes}

