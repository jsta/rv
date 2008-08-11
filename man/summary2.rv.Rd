\name{summary.rv}
\title{Numerical Summary of a Random Variable}
\alias{summary.rv}
\alias{summary.rvfactor}
\alias{summary.rvmixed}
\alias{summary.rvsummary}
\usage{
  \method{summary}{rv}(object, \dots)
  \method{summary}{rvfactor}(object, all.levels=TRUE, \dots)
%\method{summary}{rvmixed}(object, all.levels=FALSE, \dots)
  \method{summary}{rvsummary}(object, \dots)
  \method{summary}{rvsummary_numeric}(object, \dots)
  \method{summary}{rvsummary_integer}(object, \dots)
  \method{summary}{rvsummary_logical}(object, \dots)
  \method{summary}{rvsummary_rvfactor}(object, all.levels=TRUE, \dots)
}
\description{
  Gives a numerical summary of the random variable in the format of a data frame. 
}
\arguments{
  \item{object}{object to summarize}
  \item{all.levels}{show summary for all levels even if there are too many to display in one line}
  \item{\dots}{rv object}
}
\details{
  The objects are first coerced to \code{rvsummary} objects,
  then passed on to the \code{summary.rvsummary} method, 
  which creates a nicely formatted data frame of the object.
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
  x <- rvarray(rvnorm(6), c(2,3))
  summary(x)
  summary(as.rvsummary(x))
  summary(rvfactor(trunc(x)))
}

\keyword{internal}

