\name{rvbern}
\title{Generate a Random Vector from a Bernoulli Sampling Model}
\alias{rvbern}
\usage{
  rvbern(n=1, prob)
}
\description{
  \code{rvbern} generates a random vector where each
  simulation comes from a Bernoulli sampling distribution.
}
\arguments{
  \item{n}{number of random scalars to draw}
  \item{prob}{probability of ``success"; may be a random vector itself}
}
\details{
  \code{rvbern} is a special case of \code{rvbinom} with the argument size=1.
}
\value{
  A random vector (an rv object) of length \code{n}.
}
\note{
  The resulting vector will not be 
  independent and identically distributed Bernoulli unless
  \code{prob} is a fixed number.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  rvbern(10, prob=0.5)
  rvbinom(10, size=1, prob=0.5) # Equivalent
}

\keyword{classes}

