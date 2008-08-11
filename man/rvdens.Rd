\name{rvdens}
\title{Sample from an arbitrary density function using grid approximation}
\alias{rvdens}
\usage{
  rvdens(n=1, FUN, range, unitprecision=10, \dots)
}
\description{
  \code{rvdens} generates a random vector where each
  simulation comes from a Bernoulli sampling distribution.
}
\arguments{
  \item{n}{number of random scalars to draw}
  \item{FUN}{density function}
  \item{range}{range to discretize over}
  \item{unitprecision}{how many points per unit length}
  \item{\dots}{other arguments passed on to \code{FUN}}
}
\details{
  %The first argument of \code{FUN} must be the ...
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
  #
}

\keyword{classes}

