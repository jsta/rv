\name{outer.rv}
\title{Outer Product of Random Arrays}
\alias{outer.rv}
\usage{
outer.rv(X, Y=NULL, FUN="*", \dots)
}
\description{
  \code{outer.rv}
}
\arguments{
  \item{X}{First argument for function \code{FUN}}
  \item{Y}{Second argument for function \code{FUN}; if missing, \code{X} is used instead}
  \item{FUN}{a function to use on the outer products; a character string or a function}
  \item{\dots}{optional arguments to be passed to \code{FUN}}
}
\details{
  Implements the outer product for random arrays.

  Note. \code{outer} is not a generic function; thus \code{outer(x)} will not
  work if \code{x} is an rv object.
  You must write \code{outer.rv(x)} explicitly.

  See the function \code{outer} for further details.
}
\value{
  A random array.
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

