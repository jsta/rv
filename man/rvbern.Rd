\name{rvbern}
\title{Generate a Random Vector from a Bernoulli Sampling Model}
\alias{rvbern}
\usage{
  rvbern(n=1, prob, logical=FALSE)
}
\description{
  \code{rvbern} generates a random vector where each
  simulation comes from a Bernoulli sampling distribution.
}
\arguments{
  \item{n}{number of random scalars to draw}
  \item{prob}{probability of ``success"; may be a random vector itself}
  \item{logical}{logical; return a logical random variable instead}
}
\details{
  \code{rvbern} is a special case of \code{rvbinom} with the argument size=1.

  If \code{logical} is \code{TRUE}, the function returns a logical
  random variable which has TRUE for 1, FALSE for 0. 
  (The printed summary of this object is slightly different from a regular
  continuous numeric random variable.)

}
\value{
  A random vector (an rv object) of length \code{n}.
}
\note{
  The resulting vector will not be 
  independent and identically distributed Bernoulli unless
  \code{prob} is a fixed number.
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
  rvbern(10, prob=0.5)
  rvbinom(10, size=1, prob=0.5) # Equivalent
  print(rvbern(1, 0.5))
  print(rvbern(1, 0.5, logical=TRUE)) # won't show the quantiles
  print(as.logical(rvbern(1, 0.5))) # equivalent
}

\keyword{classes}

