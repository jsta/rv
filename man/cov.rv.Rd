\name{Covariance/Variance/Standard deviation Statistic}
\title{Distribution of the Covariance Statistic of a Random Vector and Array}
\alias{cov.rv}
\alias{var.rv}
\alias{sd.rv}
\usage{
  cov.rv(x, y=NULL, \dots)
  var.rv(x, \dots)
  sd.rv(x, \dots)
}
\description{
  \code{cov}...
}
\arguments{
  \item{x}{an object}
  \item{y}{}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  \code{cov} gives the distribution (that is, a random variable object)
  of the covariance statistic, as computed by the 
  standard R function \code{\link{cov}}. 

This is implemented simply by applying the numerical 
\code{cov} function to the rows of the simulation matrices
of \code{x} and \code{y} and forming a new rv object
from the resulting vector of simulations.

Similarly, \code{var} computes the variance.

%%This has been implemented by ``trapping" the original function \code{cov}
%%since \code{cov} is not a generic method.
}
\value{
  A random vector or array.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\seealso{
  See \code{\link[stats]{cov}} for details.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
}

\keyword{internal}

