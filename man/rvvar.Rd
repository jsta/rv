\name{rvvar}
\title{Variances of Components of Random Vectors}
\alias{rvvar}
\alias{rvsd}
\alias{rvsd.rv}
\alias{rvsd.rvsummary}
\alias{rvvar.rv}
\alias{rvvar.rvsummary}
\usage{
  rvvar(x)
  \method{rvvar}{rv}(x)
  \method{rvvar}{rvsummary}(x)
  rvsd(x)
  \method{rvsd}{rv}(x)
  \method{rvsd}{rvsummary}(x)
}
\description{
  Computes variances of the simulations of components of a random vector of array.
}
\arguments{
  \item{x}{an object}
}
\details{
  \code{rvvar} computes the means of the simulations 
  of all individual components of a random vector (rv) object.

That is, \code{rvvar} applies the function \code{var} to 
the vector of simulations of each component of \code{x},
thus computing "columnwise" variances of the
matrix of simulations of \code{x}.

\code{rvsd} applies the function \code{sd} to 
the vector of simulations of each component of \code{x},
thus computing "columnwise" standard deviations of the
matrix of simulations of \code{x}.
}
\value{
  A numeric vector or array (of the same dimension as that of \code{x})
}
\seealso{
  \code{\link{rvmin}}, \code{\link{rvmax}}, \code{\link{rvmedian}}, 
  \code{\link{rvsd}}.
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
  x <- rvnorm(mean=0, var=1:10)
  rvvar(x)
  rvsd(x)
}

\keyword{classes}

