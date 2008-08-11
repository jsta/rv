\name{range.rv}
\title{Distribution of the Range of a Random Vector}
\alias{range.rv}
\usage{
  \method{range}{rv}(\dots, na.rm=FALSE, finite=FALSE)
}
\description{
  \code{range.rv}
  returns a 2-component random vector containing 
  the distributions of the minimum and the maximum values of all the given arguments.
}
\arguments{
  \item{\dots}{further arguments passed to or from other methods}
  \item{na.rm}{logical, indicating if \link{NA}s should be omitted}
  \item{finite}{logical, indicating if all non-finite elements should be omitted}
}
\details{
  This is the rv-compatible version of the function \code{\link{range}}.
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
\seealso{
  \code{\link{quantile.rv}}
}
\examples{
  x <- rvnorm(mean=1:10, sd=1)
  print(range(x))
  print(quantile(x, c(0,1)))
}

\keyword{classes}

