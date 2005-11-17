\name{rvvar}
\title{Variances of Components of Random Vectors}
\alias{rvvar}
\usage{
rvvar(x, \dots)
}
\description{
  \code{rvvar}
}
\arguments{
  \item{x}{an object}
  \item{\dots}{further arguments passed to \code{var}}
}
\details{
  \code{rvvar} computes the means of the simulations 
  of all individual components of a random vector (rv) object.
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
  x <- rvnorm(mean=0, sd=1:10)
  print(rvvar(x))
}

\keyword{classes}

