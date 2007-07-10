\name{rvsd}
\title{Variances of Components of Random Vectors}
\alias{rvsd}
\usage{
rvsd(x, \dots)
}
\description{
  Computes the means of the simulations 
  of all individual components of a random vector (rv) object.
}
\arguments{
  \item{x}{an object}
  \item{\dots}{further arguments passed to \code{var}}
}
\details{
\code{rvsd} applies the function \code{sd} to 
the vector of simulations of each component of \code{x},
thus computing "columnwise" standard deviations of the
matrix of simulations of \code{x}.
}
\seealso{
  \code{\link{rvmin}}, \code{\link{rvmax}}, \code{\link{rvmedian}}, 
  \code{link{rvvar}}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  rvsd(rvnorm(mean=0, sd=1:10))
}

\keyword{classes}

