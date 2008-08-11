\name{as.list.rv}
\title{Coerce a random vector object to a list}
\alias{as.list.rv}
\usage{
  \method{as.list}{rv}(x, \dots)
}
\description{
  \code{as.list.rv} coerces a given \code{rv} object into a list.
}
\arguments{
  \item{x}{an rv object}
  \item{\dots}{arguments passed on to other methods}
}
\details{
  Each component of the argument is extracted into a component of an
  enclosing list, which is returned.
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
  x <- rvnorm(10)
  L <- as.list(x)
}

\keyword{classes}
