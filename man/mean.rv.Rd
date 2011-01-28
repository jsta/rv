\name{mean.rv}
\title{Distribution of the Arithmetic Mean of a Random Vector}
\alias{mean.rv}
\usage{
  \method{mean}{rv}(x, \dots)
}
\description{
  \code{mean.rv} computes the distribution of the arithmetic average of its
  argument \code{rv} object.
}
\arguments{
  \item{x}{an object}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  \code{mean} gives the distribution (that is, a random variable object)
  of the statistic \eqn{\frac{1}{n}\sum_{i=1}^n x_i}{sum x_i/n}
  (\code{sum(x)/length(x)}).

  In particular, \code{mean(x)} of a random vector \code{x} of length one 
  is equal to \code{x} as it would be in the case of numerical x.

  To find the expectation of a random vector \code{x}
  (that is, the individual means of random components in a vector),
  use \code{rvmean(x)} (same as \code{E(x)} and \code{Pr(x)}).
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\seealso{
  \code{\link{rvmean}}
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}

\examples{
  y <- rvnorm(10, mean=0, sd=1)
  m1 <- mean(y)
  m2 <- rvnorm(1, mean=0, sd=1/sqrt(10))
  print(c(m1, m2)) # should have the same distribution
}
\keyword{classes}

