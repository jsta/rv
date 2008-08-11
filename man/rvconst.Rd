\name{rvconst}
\title{Random Vector with a Point-Mass Distribution}
\alias{rvconst}
\usage{
rvconst(n=1,x)
}
\description{
  Coerces a given vector of constants into a random vector
  with 1 simulation in each component.
}
\arguments{
  \item{n}{integer: number of variables to generate}
  \item{x}{a vector of constants}
}
\details{
  Coerces a given vector of constants into a random vector
  with 1 simulation in each component.
}
\value{
  A random vector (rv object) of length \code{n}.
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

