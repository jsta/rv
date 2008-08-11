\name{as.rv.bugs}
\title{Coerce an R2WinBUGS object into Random Variable Objects}
\alias{as.rv.bugs}
\alias{as.rvsummary.bugs}
\usage{
  \method{as.rv}{bugs}(x, list.=TRUE, \dots)
  \method{as.rvsummary}{bugs}(x, list.=TRUE, \dots)
}
\description{
  \code{as.rv.bugs} coerces an \code{R2WinBUGS} object to
  a list of \code{rv} objects or to a named rv object (vector).

  \code{as.rvsummary.bugs} works similarly but coerces the resulting
  \code{rv} objects into \code{rvsummary} objects.
  
}
\arguments{
  \item{x}{a bugs (R2WinBUGS) object}
  \item{list.}{logical; return a list of \code{rv} objects instead of a single \code{rv} object (vector)?}
  \item{\dots}{(ignored)}
}
\details{
 
}
\value{
  A named \emph{list} of random vectors or a named random vector.
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

\keyword{classes}

