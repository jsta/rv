\name{Math.rv}
\title{Mathematical functions and Operators for rv Objects}
\alias{Math.rv}
\alias{Ops.rv}
\alias{!.rv}
\alias{Math.rvsim}
\alias{Ops.rvsim}
\alias{cumsum.rv}
\alias{cumprod.rv}
\alias{cummin.rv}
\alias{cummax.rv}
\usage{
  \method{Math}{rv}(x, ...)
  \method{Ops}{rv}(e1, e2)
}
\description{
  Mathematical functions and operators adapted to work with 
  random variable (rv) objects.
}
\arguments{
  \item{x}{object}
  \item{e1}{object}
  \item{e2}{object}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  The operator method preserves the names of the longer vector
  (or those of the first if the lengths match).
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
  -x
  names(x) <- paste("x[", seq_along(x), "]", sep="")
  x + 1:10
  1:2 + x
  cumsum(x)
  cumprod(exp(x))
}

\keyword{classes}
