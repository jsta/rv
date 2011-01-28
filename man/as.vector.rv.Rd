\name{as.vector.rv}
\title{Coerce an rv object}
\alias{as.vector.rv}
\usage{
  \method{as.vector}{rv}(x, mode="any")
}
\description{
  \code{as.vector.rv} coerces a given \code{rv} object into a vector;
  matrices lose their dimension attributes, but \code{rv} objects
  stay as \code{rv} objects (since they are considered to be ``vectors'').
}
\arguments{
  \item{x}{an object}
  \item{mode}{(currently not used)}
}
\details{
  \code{as.vector.rv} removes the dimension attribute and returns the rv object.
  Needed for compatibility with code that uses \code{as.vector}.
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
  x <- rvmatrix(rvnorm(10), 2, 5)
  as.vector(x)
}

\keyword{classes}

