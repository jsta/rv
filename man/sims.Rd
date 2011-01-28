\name{sims}
\title{Retrieve the Simulations of Random Vectors}
\alias{sims}
\alias{sims.rv}
\alias{sims.rvsummary}
\alias{sims.default}
\usage{
 sims(x, \dots)
 \method{sims}{default}(x, \dots)
 \method{sims}{rv}(x, dimensions=FALSE, n.sims=getnsims(), \dots)
 \method{sims}{rvsummary}(x, dimensions=FALSE, \dots)
}
\description{
  Returns the simulation matrix for the random variable object \code{x}.
}
\arguments{
  \item{x}{a random variable  object} 
  \item{n.sims}{(optional) number of simulations}
  \item{dimensions}{logical, try to preserve the dimensions of \code{x}}
  \item{\dots}{arguments passed on}
}
\details{
  \code{sims} returns the matrix of simulations for a given random variable object \code{x}. 

  The first index of the matrix indicates the number of the simulation draw (``simulations are in rows'').
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
  setnsims(n.sims=2500)
  x <- rvnorm(24)
  dim(x) <- c(2,3,4)
  dim(sims(x))                  # 2500x24
  dim(sims(x, dimensions=TRUE)) # 2500x2x3x4
}
\keyword{classes}

