\name{rvdiscrete}
\title{Generate Random Vectors from a Discrete Sampling Model}
\alias{rvdiscrete}
\usage{
  rvdiscrete(n=1, x, prob=NULL)
}
\description{
  Generates random variables from a discrete distribution.
}
\arguments{
  \item{n}{integer: number of scalars to generate}
  \item{x}{values of the distribution}
  \item{prob}{probabilities (optional, default: all equal)}
}
\details{
  Computes a random vector of length \code{n}, consisting of identicallly distributed 
  discrete random scalars with the discrete distribution with
  values \code{x} and corresponding probabilities \code{prob}.
  If \code{prob} is not given, all values are considered equally distributed. 
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  # 8 people draw a number each from 1..10 with replacement.
  # What is the probability that the highest number of the eight is "10"?
  u <- rvdiscrete(n=8, x=1:10) # 8 iid variables from the discrete uniform 1:10.
  s <- sort(u) # order distribution
  Pr(s[8]==10)
}

\keyword{classes}

