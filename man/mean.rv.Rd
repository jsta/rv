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
  \code{mean.rv} gives the distribution (that is, a random variable object)
  of the statistic \eqn{\frac{1}{n}\sum_{i=1}^n x_i}{sum x_i/n}
  (sum(x)/length(x)).

  In particular, \code{mean(x)} of a random vector \code{x} of length one 
  is equal to \code{x} as it would be in the case of numerical x.

  To find the expectation of a random vector \code{x}
  (that is, the individual means of random components in a vector),
  use \code{rvmean(x)} (same as \code{E(x)} and \code{Pr(x)}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\seealso{
  \code{\link{rvmean}}
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}

\examples{
  #
}
\keyword{classes}

