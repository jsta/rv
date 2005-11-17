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
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvnorm(30)
  quantile(x)
}

\keyword{classes}

