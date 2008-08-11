\name{cbind.rv}
\title{Combine random vectors by columns or rows}
\alias{cbind.rv}
\alias{rbind.rv}
\usage{
 \method{cbind}{rv}(\dots, deparse.level = 1)
 \method{rbind}{rv}(\dots, deparse.level = 1)
}
\description{
  Combines random vectors by columns (\code{cbind.rv}) or rows (\code{rbind.rv}).
}
\arguments{
  \item{\dots}{vectors or matrices, can be rv objects}
  \item{deparse.level}{(passed on to cbind)}
}
\details{
  See \link{cbind} and \link{rbind} for details.
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
  y <- rvnorm(10)
  cbind(x,y)
  rbind(x,y)
}

\keyword{classes}

