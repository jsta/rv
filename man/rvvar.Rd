\name{rvvar}
\title{Variances of Components of Random Vectors}
\alias{rvvar}
\usage{
  rvvar(x, \dots)
}
\description{
  Computes variances of the simulations of components of a random vector of array.
}
\arguments{
  \item{x}{an object}
  \item{\dots}{further arguments passed to \code{var}}
}
\details{
  \code{rvvar} computes the means of the simulations 
  of all individual components of a random vector (rv) object.

That is, \code{rvvar} applies the function \code{var} to 
the vector of simulations of each component of \code{x},
thus computing "columnwise" variances of the
matrix of simulations of \code{x}.
}
\value{
  A numeric vector or array (of the same dimension as that of \code{x})
}
\seealso{
  \code{\link{rvmin}}, \code{\link{rvmax}}, \code{\link{rvmedian}}, 
  \code{\link{rvsd}}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvnorm(mean=0, var=1:10)
  print(rvvar(x))
}

\keyword{classes}

