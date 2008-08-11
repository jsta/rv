\name{quantile.rv}
\title{Distribution of a Quantile of a Random Vector}
\alias{quantile.rv}
\usage{
  \method{quantile}{rv}(x, \dots)
}
\description{
  \code{quantile.rv} returns 
  the distribution of the quantile of a random vector 
  (as a random variable).
}
\arguments{
  \item{x}{an object}
\item{\dots}{further arguments passed to or from other methods}
}
\value{
  A random vector (rv object) with 
  components giving the distribution of the desired quantiles.
}
\note{
  \code{quantile.rv} does not return the simulated quantiles of
  the quantiles of the argument \code{x}.
  This is done by \code{\link{rvquantile}}.
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
  x <- rvnorm(30)
  quantile(x)
}

\keyword{classes}

